This repository contains Jupyter notebooks used for scRNA-seq data analysis starting with cell x gene matrices in paper "Single-cell transcriptional profiling of clear cell renal cell carcinoma reveals an invasive tumor vasculature phenotype." [1]

Most of the notebooks presented in this directory were taken directly and adapted from https://github.com/AllonMKlein/Pfirschke_et_al_2021 [2]. Sincere gratitude to Rapolas Zilionis for guidance on using this resource. 


Description of files:

| No. | File name | Brief description |
| -- | ------------- | ------------- |
| 1 | [1_Import_data_make_adata.ipynb](https://github.com/zvirblyte/2023_ccRCC/blob/main/notebooks/1_Import_data_make_adata.ipynb) | Importing data from storage by recursive search and concatenating into one adata file  |
| 2 | [2_qc_total_counts_mito.ipynb](https://github.com/zvirblyte/2023_ccRCC/blob/main/notebooks/2_qc_total_counts_mito.ipynb)  | Quality control: filtering of the libraries by total counts and mitochondrial gene expression fraction |
| 3 | [3_adding_clinical_and_technical_information.ipynb](https://github.com/zvirblyte/2023_ccRCC/blob/main/notebooks/3_adding_clinical_and_technical_information.ipynb)  | Addition of available clinical information from Table S1  |
| 4 | [4_umap_and_spring_with_batch_correction.ipynb](https://github.com/zvirblyte/2023_ccRCC/blob/main/notebooks/4_umap_and_spring_with_batch_correction.ipynb)  | Initial FDL and UMAP representations with batch correction using Harmony, loading to interactive SPRING environment  |
| 5 | [5_doublet_detection_in_Scrublet_using_precomputed_PCA.ipynb](https://github.com/zvirblyte/2023_ccRCC/blob/main/notebooks/5_doublet_detection_in_Scrublet_using_precomputed_PCA.ipynb) | Doublet detection using Scrublet |
| 6 | [6_ultrafine_clustering_and_RBC_removal.ipynb](https://github.com/zvirblyte/2023_ccRCC/blob/main/notebooks/6_ultrafine_clustering_and_RBC_removal.ipynb)  | Detecting RBC contaminationa and louvain clustering of the graph with ultra-high resolution parameter  |
| 7 | [7_doublet_detection_1round.ipynb](https://github.com/zvirblyte/2023_ccRCC/blob/main/notebooks/7_doublet_detection_1round.ipynb) | Calculating average doublet score and fraction of preficted doublets, removing clusters that are likely to be doublet based on these parameters and interactive exploration in SPRING viewer  |
| 8 | [8_no_dblt1_umap_batch_correction.ipynb](https://github.com/zvirblyte/2023_ccRCC/blob/main/notebooks/8_no_dblt1_umap_batch_correction.ipynb)  | Second iteration for FDL and UMAP after first round of doublet removal  |
| 9 | [9_ultrafine_clustering_2round.ipynb](https://github.com/zvirblyte/2023_ccRCC/blob/main/notebooks/9_ultrafine_clustering_2round.ipynb)  | Ultra-fine clustering of the new graph in attempt to find more doublets hidden in previous iteration  |
| 10 | [10_doublet_detection_2round.ipynb](https://github.com/zvirblyte/2023_ccRCC/blob/main/notebooks/10_doublet_detection_2round.ipynb)  | Same as in 7 |
| 11 | [11_cleanup_final_umap_batch_correction.ipynb](https://github.com/zvirblyte/2023_ccRCC/blob/main/notebooks/11_cleanup_final_umap_batch_correction.ipynb)  | Constructing the final UMAP embedding  |
| 12 | [12_cluster_final_graph.ipynb](https://github.com/zvirblyte/2023_ccRCC/blob/main/notebooks/12_cluster_final_graph.ipynb)  | Louvain, Pheno-graph and Spectral clustering of the graph, addition of this information to SPRING viewer  |
| 13 | [13_DGE_global_spectral_clusters.ipynb](https://github.com/zvirblyte/2023_ccRCC/blob/main/notebooks/13_DGE_global_spectral_clusters.ipynb)  | Differential gene expression analysis to determine cluster marker genes. Cluster vs the rest of cells, Mann-Whitney U test with Benjamini-Hochberg correction  |
| 14 | [14_annotation.ipynb](https://github.com/zvirblyte/2023_ccRCC/blob/main/notebooks/14_annotation.ipynb)  | Annotating cell types for the selected cluster configuration. Justification provided in Supplementary file Table 1. From this point on, the rest of notebooks can be used in any order  |
| 15 | [15_prepare_data_for_cellphoneDB.ipynb](https://github.com/zvirblyte/2023_ccRCC/blob/main/notebooks/15_prepare_data_for_cellphoneDB.ipynb)  | Preparing the data and meta files for CellPhoneDB analysis  |
| 16 | [16_DGE_annotated_Figure1f.ipynb](https://github.com/zvirblyte/2023_ccRCC/blob/main/notebooks/16_DGE_annotated_Figure1f.ipynb)  | Differential gene expression analysis for refined cell types and plotting of the nice heatmap in Figure 1f  |
| 17 | [17_DGE_within_vasculature_and_stromal.ipynb](https://github.com/zvirblyte/2023_ccRCC/blob/main/notebooks/17_DGE_within_vasculature_and_stromal.ipynb)  | Differential gene expression analysis within vasculature and stromal cell subpopulations. Cluster vs the rest of cells in group (i.e. vasculature), Mann-Whitney U test with Benjamini-Hochberg correction  |
| 18 | [18_DGE_tumor_vs_healthy_vasculature.ipynb](https://github.com/zvirblyte/2023_ccRCC/blob/main/notebooks/18_DGE_tumor_vs_healthy_vasculature.ipynb)  | Differential gene expression analysis for tumor vs healthy endothelial cells, Mann-Whitney U test with Benjamini-Hochberg correction. Plotting of a volcano in Figure 2b  |
| 19 | [19_DGE_tumor_vs_healthy_mesangial_vsmc_volcano.ipynb](https://github.com/zvirblyte/2023_ccRCC/blob/main/notebooks/19_DGE_tumor_vs_healthy_mesangial_vsmc_volcano.ipynb) | Differential gene expression analysis for tumor vs healthy mesangial/vSMC cells, Mann-Whitney U test with Benjamini-Hochberg correction. Plotting of a volcano in Supplementary figure S5  |
| 20 | [20_DGE_within_tumor_cells.ipynb](https://github.com/zvirblyte/2023_ccRCC/blob/main/notebooks/20_DGE_within_tumor_cells.ipynb)  | Differential gene expression analysis for tumor cell subpopulations, Mann-Whitney U test with Benjamini-Hochberg correction. Plotting of a heatmap in Supplementary figure S1b  |
| 21 | [21_sample_heterogeneity_measurement_by_entropy.ipynb](https://github.com/zvirblyte/2023_ccRCC/blob/main/notebooks/21_sample_heterogeneity_measurement_by_entropy.ipynb)  | Quantitative analysis of sample heterogeneity as measured by Shannon entropy values for major cell groups (i.e. stromal, endothelial etc)  |
| 22 | [22_explore_cellphonedb_results_and_plot.ipynb](https://github.com/zvirblyte/2023_ccRCC/blob/main/notebooks/22_explore_cellphonedb_results_and_plot.ipynb)  | Exploration of the cell-cell communication analysis results and plotting of the dot-plots presented in figures 2c, 3d, 5c and Supplementary figure 2b  |
| 23 | [Main_figures.ipynb](https://github.com/zvirblyte/2023_ccRCC/blob/main/notebooks/Main_figures.ipynb)  | Plotting of remaining main figures not plotted in previous notebooks  |
| 24 | [Supplementary_figures.ipynb](https://github.com/zvirblyte/2023_ccRCC/blob/main/notebooks/Supplementary_figures.ipynb)  | Plotting of remaining supplementary figures not plotted in previous notebooks |



[1] To be updated after publication
[2] Zilionis R, Engblom C, Pfirschke C, Savova V, Zemmour D, Saatcioglu HD, Krishnan I, Maroni G, Meyerovitz CV, Kerwin CM, Choi S, Richards WG, De Rienzo A, Tenen DG, Bueno R, Levantini E, Pittet MJ, Klein AM. Single-Cell Transcriptomics of Human and Mouse Lung Cancers Reveals Conserved Myeloid Populations across Individuals and Species. Immunity. 2019 May 21;50(5):1317-1334.e10. doi: 10.1016/j.immuni.2019.03.009. Epub 2019 Apr 9. PMID: 30979687; PMCID: PMC6620049.



