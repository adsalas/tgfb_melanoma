#####  DERIVING THE TREATMENT SIGNATURES FROM THE BULK RNASEQ ANALYSES IN A CELL LINE-WISE MANNER  #########

# Load the libraries ----
library(ComplexHeatmap)
library(dplyr)
library(RColorBrewer)

# Load the pre-computed outputs
load("path/to/folder/where/file/is/stored/outputs.RData") 

# Get the gene names for the significantly expressed genes  ----
# Keeping just the genes with "p ajusted < 0.05" and "log2FoldChange > X"
# M170117
TGFb_vs_DMSO.M170 <- rownames(res_TGFb_vs_DMSO.M170117_05padj[res_TGFb_vs_DMSO.M170117_05padj$log2FoldChange > 1, ])
TGFb_vs_DMSO.M170.order <- res_TGFb_vs_DMSO.M170117_05padj[res_TGFb_vs_DMSO.M170117_05padj$log2FoldChange > 1, ]
TGFb_vs_DMSO.M170.order <- TGFb_vs_DMSO.M170.order[order(TGFb_vs_DMSO.M170.order$log2FoldChange, decreasing = T), ]
TGFb_vs_DMSO.M170.order$ranking <- 1:dim(TGFb_vs_DMSO.M170.order)[1]

MEKi_vs_DMSO.M170 <- rownames(res_MEKi_vs_DMSO.M170117_05padj[res_MEKi_vs_DMSO.M170117_05padj$log2FoldChange > 1, ])
MEKi_vs_DMSO.M170.order <- res_MEKi_vs_DMSO.M170117_05padj[res_MEKi_vs_DMSO.M170117_05padj$log2FoldChange > 1, ]
MEKi_vs_DMSO.M170.order <- MEKi_vs_DMSO.M170.order[order(MEKi_vs_DMSO.M170.order$log2FoldChange, decreasing = T), ]
MEKi_vs_DMSO.M170.order$ranking <- 1:dim(MEKi_vs_DMSO.M170.order)[1]

TGFb.MEKi_vs_DMSO.M170 <- rownames(res_TGFb.MEKi_vs_DMSO.M170117_05padj[res_TGFb.MEKi_vs_DMSO.M170117_05padj$log2FoldChange > 1, ])
TGFb.MEKi_vs_DMSO.M170.order <- res_TGFb.MEKi_vs_DMSO.M170117_05padj[res_TGFb.MEKi_vs_DMSO.M170117_05padj$log2FoldChange > 1, ]
TGFb.MEKi_vs_DMSO.M170.order <- TGFb.MEKi_vs_DMSO.M170.order[order(TGFb.MEKi_vs_DMSO.M170.order$log2FoldChange, decreasing = T), ]
TGFb.MEKi_vs_DMSO.M170.order$ranking <- 1:dim(TGFb.MEKi_vs_DMSO.M170.order)[1]

# M130830
TGFb_vs_DMSO.M130 <- rownames(res_TGFb_vs_DMSO.M130830_05padj[res_TGFb_vs_DMSO.M130830_05padj$log2FoldChange > 1, ])
TGFb_vs_DMSO.M130.order <- res_TGFb_vs_DMSO.M130830_05padj[res_TGFb_vs_DMSO.M130830_05padj$log2FoldChange > 1, ]
TGFb_vs_DMSO.M130.order <- TGFb_vs_DMSO.M130.order[order(TGFb_vs_DMSO.M130.order$log2FoldChange, decreasing = T), ]
TGFb_vs_DMSO.M130.order$ranking <- 1:dim(TGFb_vs_DMSO.M130.order)[1]

MEKi_vs_DMSO.M130 <- rownames(res_MEKi_vs_DMSO.M130830_05padj[res_MEKi_vs_DMSO.M130830_05padj$log2FoldChange > 1, ])
MEKi_vs_DMSO.M130.order <- res_MEKi_vs_DMSO.M130830_05padj[res_MEKi_vs_DMSO.M130830_05padj$log2FoldChange > 1, ]
MEKi_vs_DMSO.M130.order <- MEKi_vs_DMSO.M130.order[order(MEKi_vs_DMSO.M130.order$log2FoldChange, decreasing = T), ]
MEKi_vs_DMSO.M130.order$ranking <- 1:dim(MEKi_vs_DMSO.M130.order)[1]

TGFb.MEKi_vs_DMSO.M130 <- rownames(res_TGFb.MEKi_vs_DMSO.M130830_05padj[res_TGFb.MEKi_vs_DMSO.M130830_05padj$log2FoldChange > 1, ])
TGFb.MEKi_vs_DMSO.M130.order <- res_TGFb.MEKi_vs_DMSO.M130830_05padj[res_TGFb.MEKi_vs_DMSO.M130830_05padj$log2FoldChange > 1, ]
TGFb.MEKi_vs_DMSO.M130.order <- TGFb.MEKi_vs_DMSO.M130.order[order(TGFb.MEKi_vs_DMSO.M130.order$log2FoldChange, decreasing = T), ]
TGFb.MEKi_vs_DMSO.M130.order$ranking <- 1:dim(TGFb.MEKi_vs_DMSO.M130.order)[1]

# Check for the sets of genes that are unique DEGs for each treatment ----
# Create a list with the genes to intersect
#  M170117
list.genes.1.UP.M170117 <- list(
  TGFb_vs_DMSO.M170 = TGFb_vs_DMSO.M170,
  MEKi_vs_DMSO.M170 = MEKi_vs_DMSO.M170,
  TGFb.MEKi_vs_DMSO.M170 = TGFb.MEKi_vs_DMSO.M170
)
#  M130830
list.genes.1.UP.M130830 <- list(
  TGFb_vs_DMSO.M130 = TGFb_vs_DMSO.M130,
  MEKi_vs_DMSO.M130 = MEKi_vs_DMSO.M130,
  TGFb.MEKi_vs_DMSO.M130 = TGFb.MEKi_vs_DMSO.M130
)
# Create a binary matrix from the list
#  M170117
matrix.genes.1.UP.M170117 <- list_to_matrix(list.genes.1.UP.M170117)
dim(matrix.genes.1.UP.M170117)
m.M170117 <- make_comb_mat(matrix.genes.1.UP.M170117)
#  M130830
matrix.genes.1.UP.M130830 <- list_to_matrix(list.genes.1.UP.M130830)
dim(matrix.genes.1.UP.M130830)
m.M130830 <- make_comb_mat(matrix.genes.1.UP.M130830)

# Make the UpSet plot representation
#  M170117
UpSet(m.M170117, 
      comb_col = brewer.pal(8, "Set2")[comb_degree(m.M170117)]
)

#  M130830
UpSet(m.M130830, 
      comb_col = brewer.pal(8, "Set2")[comb_degree(m.M130830)]
)

# Check the structure of the combination matrix to determine the column position of the conditions of interest
m.M170117
m.M130830

# Subset the matrix containing only unique sets
m.M170117.unique <- m.M170117[comb_degree(m.M170117) == 1]
m.M170117.unique

m.M130830.unique <- m.M130830[comb_degree(m.M130830) == 1]
m.M130830.unique

# Extract the full gene sets for each comparison
#############  M170117  #############
#             TGFB VS DMSO
# Extract the genes
TGFb_vs_DMSO.M170117.signature <- extract_comb(m.M170117.unique, "100")
TGFb_vs_DMSO.M170117.signature
# Get the gene symbols
TGFb_vs_DMSO.M170117.signature.symbol <- ext.gene.names[ext.gene.names$ensembl_gene_id %in% TGFb_vs_DMSO.M170117.signature, ] %>% pull(external_gene_name)
TGFb_vs_DMSO.M170117.signature.symbol <- TGFb_vs_DMSO.M170117.signature.symbol[TGFb_vs_DMSO.M170117.signature.symbol != ""]
TGFb_vs_DMSO.M170117.signature.symbol

#             MEKi VS DMSO
# Extract the genes
MEKi_vs_DMSO.M170117.signature <- extract_comb(m.M170117.unique, "010")
MEKi_vs_DMSO.M170117.signature
# Get the gene symbols
MEKi_vs_DMSO.M170117.signature.symbol <- ext.gene.names[ext.gene.names$ensembl_gene_id %in% MEKi_vs_DMSO.M170117.signature, ] %>% pull(external_gene_name)
MEKi_vs_DMSO.M170117.signature.symbol <- MEKi_vs_DMSO.M170117.signature.symbol[MEKi_vs_DMSO.M170117.signature.symbol != ""]
MEKi_vs_DMSO.M170117.signature.symbol

#             TGFb.MEKi VS DMSO
# Extract the genes
TGFb.MEKi_vs_DMSO.M170117.signature <- extract_comb(m.M170117.unique, "001")
TGFb.MEKi_vs_DMSO.M170117.signature
# Get the gene symbols
TGFb.MEKi_vs_DMSO.M170117.signature.symbol <- ext.gene.names[ext.gene.names$ensembl_gene_id %in% TGFb.MEKi_vs_DMSO.M170117.signature, ] %>% pull(external_gene_name)
TGFb.MEKi_vs_DMSO.M170117.signature.symbol <- TGFb.MEKi_vs_DMSO.M170117.signature.symbol[TGFb.MEKi_vs_DMSO.M170117.signature.symbol != ""]
TGFb.MEKi_vs_DMSO.M170117.signature.symbol

#############  M130830  #############
#             TGFB VS DMSO
# Extract the genes
TGFb_vs_DMSO.M130830.signature <- extract_comb(m.M130830.unique, "100")
TGFb_vs_DMSO.M130830.signature
# Get the gene symbols
TGFb_vs_DMSO.M130830.signature.symbol <- ext.gene.names[ext.gene.names$ensembl_gene_id %in% TGFb_vs_DMSO.M130830.signature, ] %>% pull(external_gene_name)
TGFb_vs_DMSO.M130830.signature.symbol <- TGFb_vs_DMSO.M130830.signature.symbol[TGFb_vs_DMSO.M130830.signature.symbol != ""]
TGFb_vs_DMSO.M130830.signature.symbol

#             MEKi VS DMSO
# Extract the genes
MEKi_vs_DMSO.M130830.signature <- extract_comb(m.M130830.unique, "010")
MEKi_vs_DMSO.M130830.signature
# Get the gene symbols
MEKi_vs_DMSO.M130830.signature.symbol <- ext.gene.names[ext.gene.names$ensembl_gene_id %in% MEKi_vs_DMSO.M130830.signature, ] %>% pull(external_gene_name)
MEKi_vs_DMSO.M130830.signature.symbol <- MEKi_vs_DMSO.M130830.signature.symbol[MEKi_vs_DMSO.M130830.signature.symbol != ""]
MEKi_vs_DMSO.M130830.signature.symbol

#             TGFb.MEKi VS DMSO
# Extract the genes
TGFb.MEKi_vs_DMSO.M130830.signature <- extract_comb(m.M130830.unique, "001")
TGFb.MEKi_vs_DMSO.M130830.signature
# Get the gene symbols
TGFb.MEKi_vs_DMSO.M130830.signature.symbol <- ext.gene.names[ext.gene.names$ensembl_gene_id %in% TGFb.MEKi_vs_DMSO.M130830.signature, ] %>% pull(external_gene_name)
TGFb.MEKi_vs_DMSO.M130830.signature.symbol <- TGFb.MEKi_vs_DMSO.M130830.signature.symbol[TGFb.MEKi_vs_DMSO.M130830.signature.symbol != ""]
TGFb.MEKi_vs_DMSO.M130830.signature.symbol

# Export the signatures as an R object
treatment.signatures.per.CellLine <- list("TGFb_M170117.signature" = TGFb_vs_DMSO.M170117.signature.symbol,
                                         "MEKi_M170117.signature" = MEKi_vs_DMSO.M170117.signature.symbol,
                                         "TGFb.MEKi_M170117.signature" = TGFb.MEKi_vs_DMSO.M170117.signature.symbol,
                                         "TGFb_M130830.signature" = TGFb_vs_DMSO.M130830.signature.symbol,
                                         "MEKi_M130830.signature" = MEKi_vs_DMSO.M130830.signature.symbol,
                                         "TGFb.MEKi_M130830.signature" = TGFb.MEKi_vs_DMSO.M130830.signature.symbol)

# Export the object
saveRDS(treatment.signatures.per.CellLine, "path/to/folder/where/file/will/be/stored/treatment_signatures_per_CellLine.RDS")

# Retrieve the session information
# sessionInfo()

# R version 4.0.2 (2020-06-22)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: OS X  13.3
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] RColorBrewer_1.1-2   dplyr_1.0.5          ComplexHeatmap_2.4.3
# 
# loaded via a namespace (and not attached):
#   [1] locfit_1.5-9.4              Rcpp_1.0.6                  lattice_0.20-41            
# [4] circlize_0.4.13             png_0.1-7                   assertthat_0.2.1           
# [7] utf8_1.2.1                  R6_2.5.0                    GenomeInfoDb_1.24.2        
# [10] stats4_4.0.2                RSQLite_2.2.7               ggplot2_3.3.3              
# [13] pillar_1.6.0                GlobalOptions_0.1.2         zlibbioc_1.34.0            
# [16] rlang_0.4.10                annotate_1.66.0             blob_1.2.1                 
# [19] S4Vectors_0.26.1            GetoptLong_1.0.5            Matrix_1.3-2               
# [22] splines_4.0.2               BiocParallel_1.22.0         geneplotter_1.66.0         
# [25] munsell_0.5.0               RCurl_1.98-1.3              bit_4.0.4                  
# [28] DelayedArray_0.14.1         compiler_4.0.2              pkgconfig_2.0.3            
# [31] BiocGenerics_0.34.0         shape_1.4.6                 tidyselect_1.1.0           
# [34] SummarizedExperiment_1.18.2 tibble_3.1.1                GenomeInfoDbData_1.2.3     
# [37] IRanges_2.22.2              matrixStats_0.58.0          XML_3.99-0.6               
# [40] fansi_0.4.2                 crayon_1.4.1                bitops_1.0-7               
# [43] xtable_1.8-4                gtable_0.3.0                lifecycle_1.0.0            
# [46] DBI_1.1.1                   magrittr_2.0.1              scales_1.1.1               
# [49] cachem_1.0.4                XVector_0.28.0              genefilter_1.70.0          
# [52] ellipsis_0.3.1              generics_0.1.0              vctrs_0.3.7                
# [55] rjson_0.2.20                tools_4.0.2                 bit64_4.0.5                
# [58] Biobase_2.48.0              glue_1.4.2                  DESeq2_1.28.1              
# [61] purrr_0.3.4                 parallel_4.0.2              fastmap_1.1.0              
# [64] survival_3.2-11             clue_0.3-60                 AnnotationDbi_1.50.3       
# [67] colorspace_2.0-0            cluster_2.1.2               GenomicRanges_1.40.0       
# [70] memoise_2.0.0



