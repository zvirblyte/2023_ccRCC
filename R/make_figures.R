
# libraries
library(SummarizedExperiment)
library(tidyverse)
library(patchwork)
library(rstatix)
library(scales)
library(ggpubr)
library(ggdist)

## load data
##

# exprs data
data_list <- readRDS(
  "data/data_list.rds"
)

# sig scores and meta data
kidney_data_list <- readRDS(
  "data/kidney_data_list.rds"
)

# survival results
surv_res <- readRDS(
  "data/surv_res.rds"
)

# fitted survival data
surv_fit_tb <- readRDS(
  "data/surv_fit_tb.rds"
)

# ora the results
ora_list_hallmark <- readRDS(
  "data/ora_list_hallmark.rds"
)

# load sig scores obj
sig_scores <- readRDS(
  "data/sig_scores.rds"
  )

###
### make plots
###

##
## survival plots
##

# define plot vars

signatures <- c(
  "Tumor_immune", "Endothelial_immune", "Stromal_immune",
  "Tumor_AVR_like_vasculature.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "Tumor_vasculature_3.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "Mesangial_vSMCs.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "Myofibroblasts.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "vSMCs.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "Tumor_vasculature_1.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "Tumor_vasculature_2.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "Tumor_vasculature_4.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "TAM_1", "TAM_2", "TAM_3", "TAM_4"
)

cohorts <- surv_fit_tb$cohort %>% 
  unique()

# make plot list
p_list <- vector(
  "list", length = length(signatures)
  )
names(p_list) <- signatures

p_list <- lapply(p_list, function(x){
  x <- vector("list", length = length(cohorts))
  names(x) <- cohorts
  return(x)
})

layout <- "
  aaaa
  aaaa
  aaaa
  aaaa
  bbbb
  "

# define layout for KM plot
for(signature in signatures){
  for(cohort in cohorts){
    
    # KM plot
    p1 <- surv_fit_tb %>% 
      filter(
        cohort %in% {{cohort}} & signature %in% {{signature}}
        ) %>%
      ggplot(
        aes(x=time, y=surv, fill = strata, color=strata)
        ) +
      geom_step() + 
      geom_point(shape="+") +
      geom_text(
        data = surv_res %>% 
          filter(cohort %in% {{cohort}} & signature %in% {{signature}} & term == "KM-logrank"),
        aes(label = paste("p[log-rank] ==", scientific(p.value, digits = 3)), x=1, y=0.1),
        inherit.aes = FALSE, parse = TRUE, size = 4,
        hjust=0, vjust = 1
      ) +
      geom_text(
        data = surv_res %>%
          filter(cohort %in% {{cohort}} & signature %in% {{signature}} & term== "Cox-Wald"),
        aes(label = paste("p[Cox] ==", scientific(p.value, digits = 3)), x=1, y=0.1),
        inherit.aes = FALSE, parse = TRUE, size = 4,
        hjust=0, vjust = -.2
      ) +
      theme(
        legend.title = element_blank()
      ) +
      ylim(c(0,1)) +
      labs(
        title = "TCGA - KIRC cohort",
        subtitle = case_when(
          grepl("_i|_e", signature) ~ paste(
            tolower(gsub("_([i,e])", " - \\1", signature)), "cell interaction signature score"
            ),
          grepl("EPIT", signature) ~ paste(
            gsub("_", " ", gsub("\\..*", " ", signature)),
            "EMT signature score"
            ),
          grepl("HYPO", signature) ~ paste(
            gsub("_", " ", gsub("\\..*", " ", signature)),
            "hypoxia signature score"
          ),
          TRUE ~ paste(
            gsub("_", " ", signature), "cell signature score"
            )
        ),
        x="Months",
        y="Survival probability"
      )
    
    add_events <- function(x, group){
        bind_rows(
          x, tibble(
            time = rep(c(50, 100, 150),2),
            strata = rep(c("high", "low"), each=3)
          )
        ) %>% 
        arrange(strata, time) %>% 
        fill(n.risk, .direction = "up") %>% 
        filter(
          time %in% c(0, 50, 100, 150)
        )
    }
    
    # table of events
    p2 <- surv_fit_tb %>% 
      filter(cohort %in% {{cohort}} & signature %in% {{signature}}) %>%
      add_events() %>% 
      mutate(
        n.risk = case_when(
          {{cohort}} %in% "tgca" & !{{signature}} %in% "TAM_4" & time == 150 & strata == "high" ~ 0,
          {{cohort}} %in% "tgca" & {{signature}} %in% "TAM_4" & time == 150 & strata == "low" ~ 0,
          TRUE ~ n.risk
        )
      ) %>% 
      ggplot(aes(time, strata)) +
      geom_text(
        aes(color = strata,label = n.risk),
        show.legend = FALSE
      ) +
      scale_x_continuous(
        sec.axis = dup_axis(),
        limits = c(0, max(surv_fit_tb[surv_fit_tb$cohort %in% cohort,]$time))
      ) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.y.left = element_text(hjust = 0)
      ) +
      labs(
        x="Number at risk (number of events)",
        y=" "
      )
    # add plot to a list
    p_list[[signature]][[cohort]] <- p1 + p2 + 
      plot_layout(design = layout, tag_level = 'new') &
      plot_annotation(
        theme = theme(plot.background = element_rect(fill ="transparent"))
        )
  }
}

###
### output figures
###

# output TCGA 
for(signature in signatures){
  file <- ifelse(
    grepl("EPIT", signature),
    paste0("tcga_", gsub("\\..*", "_EMT", tolower(signature)), "_signature"),
    paste0("tcga_", gsub(".hallmark", "", tolower(signature)), "_signature")
  )
  ggsave(
    p_list[[signature]]$tgca,
    width = 4, height = 4, dpi = 600, bg = "transparent",
    filename = paste0(file, ".pdf")
  )
}

###
### boxplots
###

##
## signature scores

# sig scores in tcga stages
#

tgca_stage <- kidney_data_list %>% 
  lapply(function(x){
    x[signatures] %>% 
      bind_rows(.id = "signature")
  }) %>% 
  bind_rows(.id = "cohort") %>% 
  filter(
    !is.na(stage)
  ) %>% 
  mutate(
    signature = case_when(
      grepl("_i|_e", signature) ~ paste(
        gsub("_([i,e])", " - \\U\\1", signature, perl = TRUE), "interaction"
      ),
      grepl("EPIT", signature) ~ paste(
        gsub("_", " ", gsub("\\..*", " ", signature)),
        "EMT"
      ),
      grepl("HYPO", signature) ~ paste(
        gsub("_", " ", gsub("\\..*", " ", signature)),
        "hypoxia"
      ),
      TRUE ~ gsub("_", " ", signature)
    ),
    stage = gsub("Stage ", "", stage)
  )

stat_tgca <- tgca_stage %>% 
  group_by(signature) %>%
  wilcox_test(
    score ~ stage, paired = F, p.adjust.method = "BH"
    ) %>% 
  ungroup() %>% 
  add_xy_position(
    step.increase = 0.2,
    scales = "free_y"
  ) %>%
  mutate(
    p.signif = case_when(
      0.01 < p & p <= 0.05 ~ "*",
      0.001 < p & p <= 0.01 ~ "**",
      0 < p & p <= 0.001 ~ "***",
      TRUE ~ "ns"
    )) %>% 
  filter(!p.adj.signif %in% "ns") %>% 
  ungroup()


p1 <- tgca_stage %>% 
  filter(grepl("Tumor - Immune", signature)) %>% 
  ggplot(
    aes(stage, score, fill= stage) 
  ) +
  stat_slab(aes(thickness = stat(pdf*n)), scale = 0.7) +
  stat_dotsinterval(aes(color=stage), side = "bottom", scale = 0.7, slab_size = NA) +
  stat_pvalue_manual(
    data=stat_tgca %>% 
      filter(grepl("Tumor - Immune", signature)), 
    label = "p.adj.signif", inherit.aes = F
  ) +
  theme(
    strip.text.x = element_text(size = 12)) +
  labs(
    title = "Tumor - Immune Cell Interaction Signature", 
    y = "signature score"
  )

ggsave(
  p1,
  filename = "Fig2e.pdf",
  width = 4,
  height = 4.5,
  dpi = 600,
  bg = "white"
)

p2 <- tgca_stage %>% 
  filter(grepl("Stromal - Immune", signature)) %>% 
  ggplot(
    aes(stage, score, fill= stage) 
  ) +
  stat_slab(aes(thickness = stat(pdf*n)), scale = 0.7) +
  stat_dotsinterval(aes(color=stage), side = "bottom", scale = 0.7, slab_size = NA) +
  stat_pvalue_manual(
    data=stat_tgca %>% 
      add_xy_position(
        step.increase = 0.35,
        scales = "free_y"
        ) %>% 
      filter(grepl("Stromal - Immune", signature)), 
    label = "p.adj.signif", inherit.aes = F,
    tip.length = 0.01
  ) +
  theme(
    strip.text.x = element_text(size = 12)
  ) +
  labs(
    title = "Stromal - Immune Cell Interaction Signature", 
    y = "signature score"
  )

ggsave(
  p2,
  filename = "Fig5e.pdf",
  width = 4,
  height = 4.5,
  dpi = 600,
  bg = "white"
)

###
### ORA figures
###

# results to table for plotting (top 10 in each element)
ora_list_hallmark_tb <- ora_list_hallmark %>%
  lapply(., as_tibble) %>%
  bind_rows(.id = "cell_type") %>%
  filter(p.adjust<0.05) %>% 
  filter(!grepl("TAM", cell_type)) %>% 
  group_by(cell_type) %>% 
  top_n(-10, wt=p.adjust) %>% 
  mutate(
    cell_type = gsub("_", " ", cell_type),
    cell_type = gsub("Tumor AVR like vasculature", "AVR like TV", cell_type),
    cell_type = as.factor(gsub("Tumor vasculature", "TV", cell_type)),
    Description = gsub("_", " ", gsub("HALLMARK_", "", Description))
  ) %>% 
  ungroup()

p3 <- ora_list_hallmark_tb %>% 
  dplyr::rename(gene_count = "Count") %>%
  ggplot(
    aes(
      x = cell_type,
      y = fct_reorder(Description, as.numeric(cell_type))
    )
  ) +
  geom_point(
    aes(size = gene_count, fill = p.adjust),
    shape = 21, color = "black"
    ) +
  scale_fill_viridis_c(
    "FDR",
    direction = 1,
    begin = 0.15,
    option = 3
  ) +
  guides(
    size=guide_legend("Gene Count"),
    fill = guide_colorbar(
      frame.colour = "black",
      frame.linewidth = 1,
      barwidth = unit(1, "line")
    )
  ) +
  theme(
    axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90),
    panel.border = element_rect(size = 1, fill = NA),
    plot.title = element_text(hjust = 1),
    plot.subtitle = element_text(hjust = 1),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"),
    axis.line = element_blank(),
    axis.title = element_blank()
  )

ggsave(
  p3, 
  filename = "Fig4a.pdf",
  width = 7, height = 6, dpi = 600,
  limitsize = FALSE
)





