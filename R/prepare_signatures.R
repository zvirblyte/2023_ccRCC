
# libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(biomaRt)
library(msigdbr)

##
## load data

# exprs data
data_list <- readRDS(
  "data/data_list.rds"
  )

##
## gene lists
##

## cell interaction
##

##
## this information is also provided in the
## Supplementary Table S7

interaction_list <- openxlsx::read.xlsx(
  "data/cell_interaction_summary.xlsx",
  startRow = 1,
  sheet = 1,
  check.names = FALSE
  ) %>% 
  gather(interaction, gene) %>% 
  filter(!is.na(gene)) %>% 
  mutate(
    interaction = gsub("[.]-[.]", "_", interaction),
    direction = 1,
    gene = case_when(
      gene %in% "CCL3L1" ~ "CCL3L3",
      TRUE ~ gene
    )
  ) %>% 
  split(.$interaction)

# check if all genes overlap tcga
lapply(seq(data_list), function(i){
  lapply(interaction_list, function(x){
    setdiff(x$gene, rownames(data_list[[i]]))
  }) %>% unlist() %>% unique() 
})

## top 100 DE genes per cell type
##
gene_list <- openxlsx::read.xlsx(
  "data/lists_enriched_genes_top_100_sp_cl.xlsx",
  startRow = 2, check.names = TRUE
  ) %>% 
  dplyr::select(-starts_with("X")) %>% 
  rename_all( ~gsub("[.]", "_", .x)) %>%
  rename_all(~gsub("_$|__", "",.x)) %>% 
  gather(cell_type, gene, -Annotated_cell_type) %>% 
  filter(Annotated_cell_type<100) %>% # top 100 genes
  # manually rename gene names, which do not match
  mutate(
    gene = case_when(
      gene %in% "CTGF" ~ "CCN2",
      gene %in% "SMIM25" ~ "PELATON",
      gene %in% "CCL3L1" ~ "CCL3L3",
      gene %in% "CYR61" ~ "CCN1",
      gene %in% "WARS" ~ "WARS1",
      gene %in% "CXorf36" ~ "DIPK2B",
      gene %in% "AC133644.2" ~ "LINC01943",
      gene %in% "HIST1H4C" ~ "H4C3",
      gene %in% "H2AFV" ~ "H2AZ2",
      gene %in% "H2AFZ" ~ "H2AZ1",
      gene %in% "C2orf40" ~ "ECRG4",
      gene %in% "AC020916.1" ~ "MIR24-2",
      gene %in% "AC020571.1" ~ "HECW2-AS1",
      gene %in% "AC114760.2" ~ "H4C3",
      gene %in% "KIAA1551" ~ "RESF1",
      gene %in% "H2AFX" ~ "H2AX",
      gene %in% "H1F0" ~ "H1-0",
      gene %in% "HIST1H1C" ~ "H1-2",
      gene %in% "H3F3B" ~ "H3-3B",
      gene %in% "SEPT4" ~ "SEPTIN4",
      gene %in% "SEPT7" ~ "SEPTIN7",
      gene %in% "H2AFJ" ~ "H2AJ",
      gene %in% "AES" ~ "TLE5",
      gene %in% "RARRES3" ~ "PLAAT4",
      gene %in% "H3F3A" ~ "H3-3A",
      gene %in% "HIST1H2AE" ~ "H2AC8",
      gene %in% "HIST1H2BG" ~ "H2BC8",
      gene %in% "FAM129A" ~ "NIBAN1",
      gene %in% "FAM49B" ~ "CYRIB",
      gene %in% "ADSSL1" ~ "ADSS1",
      gene %in% "PLA2G16" ~ "PLAAT3",
      gene %in% "MINOS1" ~ "MICOS10",
      TRUE ~ gene
    )
  ) %>% 
  filter(
    grepl("Tumor|vSMCs|Myofibroblasts|TAM", cell_type)
  ) %>% 
  split(.$cell_type)

# check if all genes overlap
lapply(seq(data_list), function(i){
  lapply(gene_list, function(x){
    setdiff(x$gene, rownames(data_list[[i]]))
  }) %>% unlist() %>% unique() 
})

##
## ORA signatures

# non-zero genes in kidney scRNA-seq data
non_zero <- read_csv(
  "data/all_nonzero_genes.csv.gz"
  ) %>% 
  pull(`0`)

non_zero_ens <- read_tsv(
  "data/features_ref.tsv.gz",
  col_names = FALSE,
  ) %>% 
  dplyr::select(-X3) %>% 
  dplyr::rename(
    "ensembl_gene_id" = X1,
    "hgnc_symbol" = X2 
  ) %>% 
  filter(
    hgnc_symbol %in% non_zero
  )

# get entrez IDs for universe
ensembl <- useEnsembl(
  biomart = "genes"
  )
ensembl <- useDataset(
  dataset = "hsapiens_gene_ensembl",
  mart = ensembl
  )
universe <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
  filters="ensembl_gene_id",
  values = non_zero_ens$ensembl_gene_id,
  mart=ensembl
) 

# translate to entrez IDs
entrez_list <- lapply(gene_list, function(x){
  bitr(
    x$gene,
    fromType="SYMBOL",
    toType="ENTREZID",
    OrgDb = org.Hs.eg.db,
    drop = FALSE
    ) %>%
    na.omit() %>% 
    dplyr::rename(gene = "SYMBOL")
})

# MSigDB hallmark signatures
ms_hallmark <- msigdbr(
  species = "Homo sapiens", category = "H"
  ) %>% 
  dplyr::select(gs_name, entrez_gene)

# perform enrichment analysis
ora_list_hallmark <- lapply(entrez_list, function(x){
  enricher(
    gene =unique(x$ENTREZID),
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    TERM2GENE=ms_hallmark,
    universe = as.character(universe$entrezgene_id)
  )
})

# save the results
saveRDS(
  ora_list_hallmark, "data/ora_list_hallmark.rds"
  )

# extract genes from the top 10 
ora_hallmark_tsig <-ora_list_hallmark %>% 
  lapply(., function(x){
    x %>% as_tibble() %>% 
      filter(p.adjust<0.05) %>% 
      top_n(-10, wt=p.adjust) %>% 
      split(.$Description) %>% 
      lapply(., function(y){
        y %>% .$geneID %>% 
          strsplit(., "/") %>% 
          unlist() %>% 
          bitr(
            fromType="ENTREZID",
            toType="SYMBOL",
            OrgDb = org.Hs.eg.db,
            drop = FALSE
          ) %>% dplyr::rename(
            gene = "SYMBOL"
          )
      })
  }) %>% unlist(recursive = FALSE)

###
### scaled and centered (z-score) median expression
###

# cell interactions
sig_scores_int <- lapply(interaction_list, function(x){
  assay(data_list[[1]], "fpkm_uq_unstrand")[
    intersect(x$gene,rownames(data_list[[1]])),
    ] %>% t() %>% 
    as_tibble(rownames = "row_id") %>% 
    mutate_if(is.numeric, scale) %>% 
    column_to_rownames("row_id") %>% 
    as.matrix() %>% rowMeans(na.rm=TRUE) %>% 
    setNames(., colnames(data_list[[1]])) %>% 
    as_tibble(rownames = "row_id") %>% 
    dplyr::rename(score = "value") %>% 
    mutate(
      score_cat = case_when(
        median(score) >= score ~ "low",
        TRUE ~ "high"
      )
    )
})

# top DE genes
sig_scores_tde <- lapply(gene_list, function(x){
  assay(data_list[[1]], "fpkm_uq_unstrand")[
    intersect(x$gene,rownames(data_list[[1]])),
    ] %>% t() %>% 
    as_tibble(rownames = "row_id") %>% 
    mutate_if(is.numeric, scale) %>% 
    column_to_rownames("row_id") %>% 
    as.matrix() %>% rowMeans(na.rm=TRUE) %>% 
    setNames(., colnames(data_list[[1]])) %>% 
    as_tibble(rownames = "row_id") %>% 
    dplyr::rename(score = "value") %>% 
    mutate(
      score_cat = case_when(
        median(score) >= score ~ "low",
        TRUE ~ "high"
        )
      )
  })

# top ORA genes
sig_scores_hallmark <- lapply(ora_hallmark_tsig, function(x){
  assay(data_list[[1]], "fpkm_uq_unstrand")[
    intersect(x$gene,rownames(data_list[[1]])),
    ] %>% t() %>% 
    as_tibble(rownames = "row_id") %>% 
    mutate_if(is.numeric, scale) %>% 
    column_to_rownames("row_id") %>% 
    as.matrix() %>% rowMeans(na.rm=TRUE) %>% 
    setNames(., colnames(data_list[[1]])) %>% 
    as_tibble(rownames = "row_id") %>% 
    dplyr::rename(score = "value") %>% 
    mutate(
      score_cat = case_when(
        median(score) >= score ~ "low",
        TRUE ~ "high"
        )
      )
  })

# combine sig scores
sig_scores <- list(c(
  sig_scores_tde,
  sig_scores_int,
  sig_scores_hallmark
  )
)
names(sig_scores) <- names(data_list)

# save obj
saveRDS(
  sig_scores, "data/sig_scores.rds"
)
