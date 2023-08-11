
# libraries
library(SummarizedExperiment)
library(TCGAbiolinks)
library(tidyverse)
library(biomaRt)

##
## data
##

##
## TGCA

# link to repo
# https://portal.gdc.cancer.gov/repository

# get GDC project IDs
project_ids <- getGDCprojects()

# get kidney projects
kidney_ids <- project_ids %>% 
  filter(
    grepl("Kidney", name) & tumor == "KIRC"
    )

# query kidney data
query <- GDCquery(
  project = kidney_ids$id,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)

# download kidney data
GDCdownload(
  query = query, 
  method = "api", 
  files.per.chunk = 10
)

# prepare data
data <- GDCprepare(
  query = query,
  save = TRUE,
  save.filename = "data/kidney_data_hg38.rda"
  )

# load
load("data/kidney_data_hg38.rda")

# gene names as idx
rownames(data) <- rowData(data)$gene_name

# tumor samples
tumor_samples <- colData(data) %>% 
  as_tibble(rownames = "row_id") %>%
  # get tumor samples
  filter(
    grepl("Tumor", definition)
  ) %>% 
  pull(row_id)

# retain qced samples
data_qced <- data[,tumor_samples]

# dataset to list
data_list <- list(
  tgca = data_qced
)

# save obj
saveRDS(
  data_list,"data/data_list.rds"
  )
