
# Load the libraries
library(GEOquery)
library(Seurat)
library(dplyr)

# Define the directories ----
pathOutput = "path/to/output/"
# Set the output directory as working directory (where outputs are to be saved)
wd = pathOutput
setwd(wd)

# Computing the scores in the Rambow dataset for the signatures derived from bulk RNAseq for individual cell lines  ----
# Import the Seurat object (available in Github under the "files" folder)
seurat.rambow.subset <- readRDS("path/to/folder/where/file/is/stored/seurat_rambow_subset.RDS")
# Import the pre-computed signatures
treatment.signatures.per.CellLine <- readRDS("path/to/folder/where/file/is/stored/treatment_signatures_per_CellLine.RDS") #"path/to/folder/where/file/is/stored/"

# Set the desired order for the distinct cell states
seurat.rambow.subset$cell.state <- factor(seurat.rambow.subset$cell.state, levels = c("pigmented","SMC","invasive","NCSC"))
unique(seurat.rambow.subset$cell.state)
seurat.rambow.subset <- SetIdent(seurat.rambow.subset, value = "cell.state")
unique(Idents(seurat.rambow.subset))

######  M170117  #########
# Compute the signatures
seurat.rambow.subset <- AddModuleScore(
  object = seurat.rambow.subset,
  features = list(treatment.signatures.per.CellLine$TGFb_M170117.signature),
  ctrl = 5,
  name = 'TGFb.M170117.sig.score',
  seed = 123
)

# Plotting the score
# somePDFPath = "TGFb.M170117_sig_UpregulatedGenes_1L2FC.pdf"
# pdf(file=somePDFPath, height = 3, width = 4)
VlnPlot(seurat.rambow.subset, 
        features = "TGFb.M170117.sig.score1", 
        sort = FALSE, 
        cols = c("lightblue", "orange", "gray", "salmon")
) + NoLegend()
# dev.off()

# Compute the signatures
seurat.rambow.subset <- AddModuleScore(
  object = seurat.rambow.subset,
  features = list(treatment.signatures.per.CellLine$MEKi_M170117.signature),
  ctrl = 5,
  name = 'MEKi.M170117.sig.score',
  seed = 123
)

# Plotting the score
# somePDFPath = "MEKi.M170117_sig_UpregulatedGenes_1L2FC.pdf"
# pdf(file=somePDFPath, height = 3, width = 4)
VlnPlot(seurat.rambow.subset, 
        features = "MEKi.M170117.sig.score1", 
        sort = FALSE,
        cols = c("salmon", "gray", "orange",  "lightblue")
) + NoLegend()
# dev.off()

# Compute the signatures
seurat.rambow.subset <- AddModuleScore(
  object = seurat.rambow.subset,
  features = list(treatment.signatures.per.CellLine$TGFb.MEKi_M170117.signature),
  ctrl = 5,
  name = 'TGFb.MEKi.M170117.sig.score',
  seed = 123
)

# Plotting the score
# somePDFPath = "TGFbMEKi.M170117_sig_UpregulatedGenes_1L2FC.pdf"
# pdf(file=somePDFPath, height = 3, width = 4)
VlnPlot(seurat.rambow.subset, 
        features = "TGFb.MEKi.M170117.sig.score1", 
        sort = FALSE,
        cols = c("lightblue", "orange", "salmon", "gray")
) + NoLegend()
# dev.off()

######  M130830  #########
# Compute the signatures
seurat.rambow.subset <- AddModuleScore(
  object = seurat.rambow.subset,
  features = list(treatment.signatures.per.CellLine$TGFb_M130830.signature),
  ctrl = 5,
  name = 'TGFb.M130830.sig.score',
  seed = 123
)

# Plotting the score
# somePDFPath = "TGFb.M130830_sig_UpregulatedGenes_1L2FC.pdf"
# pdf(file=somePDFPath, height = 3, width = 4)
VlnPlot(seurat.rambow.subset, 
        features = "TGFb.M130830.sig.score1", 
        sort = FALSE,
        cols = c("lightblue", "orange", "gray", "salmon")
) + NoLegend()
# dev.off()

# Compute the signatures
seurat.rambow.subset <- AddModuleScore(
  object = seurat.rambow.subset,
  features = list(treatment.signatures.per.CellLine$MEKi_M130830.signature),
  ctrl = 5,
  name = 'MEKi.M130830.sig.score',
  seed = 123
)

# Plotting the score
# somePDFPath = "MEKi.M130830_sig_UpregulatedGenes_1L2FC.pdf"
# pdf(file=somePDFPath, height = 3, width = 4)
VlnPlot(seurat.rambow.subset, 
        features = "MEKi.M130830.sig.score1", 
        sort = FALSE,
        cols = c("salmon", "gray", "orange",  "lightblue")
) + NoLegend()
# dev.off()

# Compute the signatures
seurat.rambow.subset <- AddModuleScore(
  object = seurat.rambow.subset,
  features = list(treatment.signatures.per.CellLine$TGFb.MEKi_M130830.signature),
  ctrl = 5,
  name = 'TGFb.MEKi.M130830.sig.score',
  seed = 123
)

# Plotting the score
# somePDFPath = "TGFbMEKi.M130830_sig_UpregulatedGenes_1L2FC.pdf"
# pdf(file=somePDFPath, height = 3, width = 4)
VlnPlot(seurat.rambow.subset, 
        features = "TGFb.MEKi.M130830.sig.score1", 
        sort = FALSE,
        cols = c("lightblue", "orange", "salmon", "gray")
) + NoLegend()
# dev.off()

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
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] dplyr_1.0.5         SeuratObject_4.0.0  Seurat_4.0.1        GEOquery_2.56.0    
# [5] Biobase_2.48.0      BiocGenerics_0.34.0
# 
# loaded via a namespace (and not attached):
#   [1] nlme_3.1-152          spatstat.sparse_2.0-0 matrixStats_0.58.0    RcppAnnoy_0.0.18     
# [5] RColorBrewer_1.1-2    httr_1.4.2            sctransform_0.3.2     tools_4.0.2          
# [9] utf8_1.2.1            R6_2.5.0              irlba_2.3.3           rpart_4.1-15         
# [13] KernSmooth_2.23-18    uwot_0.1.10           mgcv_1.8-35           DBI_1.1.1            
# [17] lazyeval_0.2.2        colorspace_2.0-0      withr_2.4.2           tidyselect_1.1.0     
# [21] gridExtra_2.3         compiler_4.0.2        xml2_1.3.2            plotly_4.9.3         
# [25] labeling_0.4.2        scales_1.1.1          spatstat.data_2.1-0   lmtest_0.9-38        
# [29] readr_1.4.0           ggridges_0.5.3        pbapply_1.4-3         goftest_1.2-2        
# [33] stringr_1.4.0         digest_0.6.27         spatstat.utils_2.1-0  pkgconfig_2.0.3      
# [37] htmltools_0.5.1.1     parallelly_1.24.0     fastmap_1.1.0         limma_3.44.3         
# [41] htmlwidgets_1.5.3     rlang_0.4.10          shiny_1.6.0           farver_2.1.0         
# [45] generics_0.1.0        zoo_1.8-9             jsonlite_1.7.2        ica_1.0-2            
# [49] magrittr_2.0.1        patchwork_1.1.1       Matrix_1.3-2          Rcpp_1.0.6           
# [53] munsell_0.5.0         fansi_0.4.2           abind_1.4-5           reticulate_1.19      
# [57] lifecycle_1.0.0       stringi_1.5.3         MASS_7.3-53.1         Rtsne_0.15           
# [61] plyr_1.8.6            grid_4.0.2            listenv_0.8.0         promises_1.2.0.1     
# [65] ggrepel_0.9.1         crayon_1.4.1          deldir_0.2-10         miniUI_0.1.1.1       
# [69] lattice_0.20-41       cowplot_1.1.1         splines_4.0.2         tensor_1.5           
# [73] hms_1.0.0             pillar_1.6.0          igraph_1.2.6          spatstat.geom_2.1-0  
# [77] future.apply_1.7.0    reshape2_1.4.4        codetools_0.2-18      leiden_0.3.7         
# [81] glue_1.4.2            data.table_1.14.0     vctrs_0.3.7           png_0.1-7            
# [85] httpuv_1.6.0          polyclip_1.10-0       gtable_0.3.0          RANN_2.6.1           
# [89] purrr_0.3.4           spatstat.core_2.1-2   tidyr_1.1.3           scattermore_0.7      
# [93] future_1.21.0         assertthat_0.2.1      ggplot2_3.3.3         mime_0.10            
# [97] xtable_1.8-4          later_1.2.0           survival_3.2-11       viridisLite_0.4.0    
# [101] tibble_3.1.1          cluster_2.1.2         globals_0.14.0        fitdistrplus_1.1-3   
# [105] ellipsis_0.3.1        ROCR_1.0-11  


