
# libraries
library(SummarizedExperiment)
library(tidyverse)
library(survminer)
library(survival)
library(broom)

##
## load data

# exprs data
data_list <- readRDS(
  "data/data_list.rds"
  )

# sig score data
sig_scores <- readRDS(
  "data/sig_scores.rds"
  )

###
### survival analysis
###

## prepare data for survival analysis
##

kidney_list <- lapply(data_list, function(x){
  x %>% colData() %>% as_tibble(rownames = "row_id")
})

kidney_list$tgca <- kidney_list$tgca %>% 
  dplyr::select(
    row_id, days_to_last_follow_up, gender,
    vital_status, age_at_index, days_to_death,
    ajcc_pathologic_stage
  ) %>% 
  mutate(
    status = case_when(
      vital_status == "Alive" ~ 0,
      TRUE ~ 1
    ),
    time = case_when(
      is.na(days_to_death) ~ days_to_last_follow_up,
      TRUE ~ days_to_death
    ),
    time = time/29.66667,
    stage = ajcc_pathologic_stage
  ) %>% 
  dplyr::select(
    row_id, gender, age_at_index, status, time, stage
  )

kidney_data_list <- lapply(kidney_list, function(x){
    sig_scores$tgca %>% lapply(
      function(y){
        y %>% left_join(., x)
      })
})

saveRDS(
  kidney_data_list,
  "data/kidney_data_list.rds"
)

##
## multi-variate cox regression
## using gender and age as covariates

cox_fit_list <- lapply(kidney_data_list, function(x){
  x %>% lapply(., function(y){
    coxph(
      Surv(time, status) ~ gender + age_at_index + score, 
      data =  y
      )
    })
  })

# a positive sign of coefficient means that the hazard (risk of death) is higher, 
# and thus the prognosis worse, for subjects with higher values of that variable.

# results
cox_res_tb <- lapply(cox_fit_list, function(x){
  x %>% lapply(., function(y){
    tidy(y)
  }) %>% bind_rows(
    .id = "signature"
  ) %>% 
    group_by(term) %>% 
    mutate(
      p.adj = p.adjust(p.value, method = "BH")
    ) %>% 
    ungroup()
  }) %>% 
  bind_rows(
    .id = "cohort"
  )

# get global p-value
pval_tb <- lapply(cox_fit_list, function(x){
  x %>% lapply(., function(y){
    y %>% summary() %>% .$waldtest %>% 
      t() %>% as_tibble() %>% 
      dplyr::rename(
        "p.value" = pvalue,
        "statistic" = test
      ) %>% 
      mutate(
        term = "cox-wald"
      )
  }) %>% bind_rows(
    .id = "signature"
  ) %>% 
    mutate(
      p.adj = p.adjust(p.value, method = "BH")
    ) 
  }) %>% 
  bind_rows(
    .id = "cohort"
  )

# combine results
cox_res <- bind_rows(
  cox_res_tb, pval_tb
  ) %>% 
  arrange(
    cohort, signature, term
    )

# test proportional hazards assumption
# if p-val signif. then the variable
# deviates from assumption

assum <- lapply(cox_fit_list, function(x){
  x %>% lapply(., cox.zph)
})

## 
## Kaplan-Meier estimate and log rank test 

# get p-values
km_res <- lapply(kidney_data_list, function(x){
  x %>% lapply(., function(y){
    survdiff(
      Surv(time, status) ~ score_cat, 
      data =  y
    ) %>% glance()
  }) %>% bind_rows(
    .id = "signature"
  ) %>% 
    mutate(
      p.adj = p.adjust(p.value, method = "BH")
    ) 
  }) %>% 
  bind_rows(
    .id = "cohort"
  ) %>% 
  mutate(
    term = "KM-logrank"
  )

# combine both results
surv_res <- bind_rows(
  cox_res %>% 
    filter(term=="score"), 
  km_res
  ) %>% 
  mutate(
    term = gsub("score", "Cox-Wald", term),
  ) %>% 
  arrange(
    cohort, signature, term
  )

# save obj
saveRDS(
  surv_res, "data/surv_res.rds"
)

# fit using groups
surv_fit <- lapply(kidney_data_list, function(x){
  x %>% lapply(., function(y){
    survfit(
      Surv(time, status) ~ score_cat, 
      data =  y
    )
  })
})

# extract fitted data
surv_fit_tb <- lapply(seq(surv_fit), function(i){
  z <- lapply(seq(surv_fit[[1]]), function(j){
    surv_summary(surv_fit[[i]][[j]], kidney_data_list[[i]][[j]])
  })
  names(z) <- names(surv_fit[[i]])
  z %>% bind_rows(.id = "signature")
  }) %>% setNames(names(surv_fit)) %>% 
  bind_rows(.id = "cohort") %>% 
  mutate(
    strata = gsub(".*=", "", strata)
  )

# save obj
saveRDS(
  surv_fit_tb, "data/surv_fit_tb.rds"
)


