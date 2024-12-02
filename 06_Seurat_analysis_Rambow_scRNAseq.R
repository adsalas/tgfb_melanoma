# Reanalysis of the Rmnbow scRNAseq dataset: GSE116237

# NOTE:
#      1. The main  goal of this script is to create a Seurat object using as input the Rambow scRNAseq data, as detailed below.
#      2. Alternatively, a pre-computed Seurat object is available in Github under the "files" folder. The object can be used as input for running script 07. 


# Load the libraries ----
library(GEOquery)
library(Seurat)
library(dplyr)
library(biomaRt)

# Define the directories ----
pathData = "path/to/data/"
pathOutput = "path/to/folder/where/file/will/be/stored/"
# Set the output directory as working directory (where outputs are to be saved)
wd = pathOutput
setwd(wd)

# Retrieve the metadata using the series matrix file downloaded from GEO accession number GSE116237
rambow.metadata <- getGEO(filename= paste0(pathData, "GSE116237-GPL18573_series_matrix.txt"))
# Explore the structure of the object
str(rambow.metadata)
# Extract the required information
rambow.metadata.df <- rambow.metadata@phenoData@data
dim(rambow.metadata.df)
colnames(rambow.metadata.df)
rambow.metadata.df <- rambow.metadata.df[, c("title" ,"source_name_ch1" ,"characteristics_ch1","characteristics_ch1.1","characteristics_ch1.2","age:ch1","genotype:ch1","tissue:ch1")]
dim(rambow.metadata.df)
rambow.metadata.df

# Import the expression matrix downloaded from GEO accession number GSE116237 
rambow.data <- read.table(paste0(pathData, "GSE116237_scRNAseq_expressionMatrix.txt"), sep = ",", header = T, row.names = 1)
dim(rambow.data) 

rambow.data[1:10, 1:5]
rambow.data[1:10, 672:674]
colnames(rambow.data)
# Check if there is a match between the colnames of the count matrix and the "title" column of the metadata
table(colnames(rambow.data) %in% rambow.metadata.df$title)
# Filter the metadata to match the cells in the count matrix
rownames(rambow.metadata.df) <- rambow.metadata.df$title
rambow.metadata.df <- rambow.metadata.df[colnames(rambow.data), ]
dim(rambow.metadata.df)
# Check that the names are matching
all(colnames(rambow.data) == rownames(rambow.metadata.df))

# NOTES:
#      1. The data is already filtered, it contains 674 cells corresponding to the same number mentioned in the paper figure
#      2. The authors of the original paper analyzed the data with a customized and complex pipeline. Here we are interested only in performing scoring using our 
#         treatment derived signatures and therefore we create a basic Seurat object keeping all cells present in the input expression matrix, and use the original 
#         cluster annotations provided by the authors for our analysis
    

# Analysis with Seurat ----
# Initialize the Seurat object with the raw (non-normalized data).
seurat.rambow <- CreateSeuratObject(counts = rambow.data, project = "seurat.rambow", min.cells = 3, min.features = 0)
seurat.rambow

seurat.rambow[["percent.mt"]] <- PercentageFeatureSet(seurat.rambow, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(seurat.rambow, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(seurat.rambow, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat.rambow, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Run data normalization
seurat.rambow <- NormalizeData(seurat.rambow, normalization.method = "LogNormalize", scale.factor = 10000)

# Identifying HVGs
seurat.rambow <- FindVariableFeatures(seurat.rambow, 
                                      selection.method = "vst", 
                                      nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat.rambow), 10)
# Plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat.rambow)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Perform scaling
all.genes <- rownames(seurat.rambow)
seurat.rambow <- ScaleData(seurat.rambow, features = all.genes) 

# Run linear dimensionality reduction
seurat.rambow <- RunPCA(seurat.rambow, features = VariableFeatures(object = seurat.rambow))

# Examine and visualize PCA results a few different ways
print(seurat.rambow[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(seurat.rambow, dims = 1:2, reduction = "pca")

DimPlot(seurat.rambow, reduction = "pca")

DimHeatmap(seurat.rambow, dims = 1:20, cells = 500, balanced = TRUE)

# Deterrmine the dimensionality of the data

ElbowPlot(seurat.rambow)

# Perform clustering analysis
seurat.rambow <- FindNeighbors(seurat.rambow, dims = 1:4) 
seurat.rambow <- FindClusters(seurat.rambow, resolution = 0.1)

# Running non linear dimensionality reduction 
seurat.rambow <- RunUMAP(seurat.rambow, dims = 1:4)

DimPlot(seurat.rambow, reduction = "umap")

# Add the cell annotation from the paper
seurat.rambow@meta.data$annotation <- rambow.metadata.df$source_name_ch1
seurat.rambow@meta.data$treatment <- rambow.metadata.df$`age:ch1`

DimPlot(seurat.rambow, reduction = "umap", label = T)
DimPlot(seurat.rambow, reduction = "umap", group.by = "annotation")
DimPlot(seurat.rambow, reduction = "umap", group.by = "treatment")

# Convert the gene IDs from mouse to human ----
# Define the database used for the gene mappings
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "dec2021.archive.ensembl.org")

# Create a dataframe that has information of the gene_id and gene_name
my.genes <- rownames(seurat.rambow)
length(my.genes)

# Retrieve the data frame with the converted gene names
rambow.human <- getBM(attributes = c('ensembl_gene_id',
                                     'external_gene_name'),
                      filters = 'ensembl_gene_id', 
                      values = my.genes,
                      mart = human)

dim(rambow.human)
# Remove duplicated genes
rambow.human.FTD <- rambow.human[!duplicated(rambow.human$external_gene_name), ]
dim(rambow.human.FTD)
length(unique(rambow.human.FTD$external_gene_name))
# Remove genes with missing external names
rambow.human.FTD <- rambow.human.FTD[rambow.human.FTD$external_gene_name != "", ]
dim(rambow.human.FTD)
# Subset the Seurat object by keeping only the retrieved genes
seurat.rambow.subset <- seurat.rambow[rambow.human.FTD$ensembl_gene_id, ]
seurat.rambow.subset
all(rownames(seurat.rambow.subset) == rambow.human.FTD$ensembl_gene_id)
# Transfer the human gene names to the distinct expression matrices
rownames(seurat.rambow.subset@assays$RNA@counts) <-  rambow.human.FTD$external_gene_name
rownames(seurat.rambow.subset@assays$RNA@data) <-  rambow.human.FTD$external_gene_name
rownames(seurat.rambow.subset@assays$RNA@scale.data) <-  rambow.human.FTD$external_gene_name

# Adding the cell states annotation as provided by the authors (available in Github under the "files" folder)
# Load the annotation file
cell.states.ann <- read.table(paste0(pathData, "674_cell state_call.txt"), header = T)
cell.states.ann
dim(cell.states.ann)
tail(cell.states.ann)
# Remove the extra characters at the end of the file
tail(cell.states.ann$phenotype)
cell.states.ann[674, ]$phenotype <- "pigmented"
tail(cell.states.ann$phenotype)
tail(cell.states.ann)
rownames(cell.states.ann) <- cell.states.ann$cell
# Check if the cells are in the appropriate order and re-arrange them otherwise
all(colnames(seurat.rambow.subset) == rownames(cell.states.ann))
cell.states.ann <- cell.states.ann[colnames(seurat.rambow.subset), ]
all(colnames(seurat.rambow.subset) == rownames(cell.states.ann))
# Transfer the cell state annotation
seurat.rambow.subset$cell.state <- cell.states.ann$phenotype
seurat.rambow.subset$cell.state <- as.factor(seurat.rambow.subset$cell.state)
unique(seurat.rambow.subset$cell.state)
# Set the Idents
seurat.rambow.subset <- SetIdent(seurat.rambow.subset, value = seurat.rambow.subset[["cell.state"]])
unique(Idents(seurat.rambow.subset))
# Subset the Seurat object to include only the 4 main cell states
seurat.rambow.subset <- subset(seurat.rambow.subset, idents = c("SMC","NCSC","invasive","pigmented"))

# Export the Seurat object
saveRDS(seurat.rambow.subset, "seurat_rambow_subset.RDS")

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
#   [1] biomaRt_2.44.4      dplyr_1.0.5         SeuratObject_4.0.0  Seurat_4.0.1       
# [5] GEOquery_2.56.0     Biobase_2.48.0      BiocGenerics_0.34.0
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15            colorspace_2.0-0      deldir_0.2-10         ellipsis_0.3.1       
# [5] ggridges_0.5.3        rstudioapi_0.13       spatstat.data_2.1-0   farver_2.1.0         
# [9] leiden_0.3.7          listenv_0.8.0         bit64_4.0.5           ggrepel_0.9.1        
# [13] RSpectra_0.16-0       AnnotationDbi_1.50.3  fansi_0.4.2           xml2_1.3.2           
# [17] codetools_0.2-18      splines_4.0.2         cachem_1.0.4          polyclip_1.10-0      
# [21] jsonlite_1.7.2        ica_1.0-2             dbplyr_2.1.1          cluster_2.1.2        
# [25] png_0.1-7             uwot_0.1.10           shiny_1.6.0           sctransform_0.3.2    
# [29] spatstat.sparse_2.0-0 readr_1.4.0           compiler_4.0.2        httr_1.4.2           
# [33] assertthat_0.2.1      Matrix_1.3-2          fastmap_1.1.0         lazyeval_0.2.2       
# [37] cli_3.6.2             limma_3.44.3          later_1.2.0           prettyunits_1.1.1    
# [41] htmltools_0.5.1.1     tools_4.0.2           igraph_1.2.6          gtable_0.3.0         
# [45] glue_1.4.2            RANN_2.6.1            reshape2_1.4.4        rappdirs_0.3.3       
# [49] Rcpp_1.0.6            scattermore_0.7       vctrs_0.3.7           nlme_3.1-152         
# [53] lmtest_0.9-38         stringr_1.4.0         globals_0.14.0        mime_0.10            
# [57] miniUI_0.1.1.1        lifecycle_1.0.0       irlba_2.3.3           XML_3.99-0.6         
# [61] goftest_1.2-2         future_1.21.0         MASS_7.3-53.1         zoo_1.8-9            
# [65] scales_1.1.1          spatstat.core_2.1-2   hms_1.0.0             promises_1.2.0.1     
# [69] spatstat.utils_2.1-0  RColorBrewer_1.1-2    curl_4.3              memoise_2.0.0        
# [73] reticulate_1.19       pbapply_1.4-3         gridExtra_2.3         ggplot2_3.3.3        
# [77] rpart_4.1-15          stringi_1.5.3         RSQLite_2.2.7         S4Vectors_0.26.1     
# [81] rlang_0.4.10          pkgconfig_2.0.3       matrixStats_0.58.0    lattice_0.20-41      
# [85] ROCR_1.0-11           purrr_0.3.4           tensor_1.5            labeling_0.4.2       
# [89] patchwork_1.1.1       htmlwidgets_1.5.3     bit_4.0.4             cowplot_1.1.1        
# [93] tidyselect_1.1.0      parallelly_1.24.0     RcppAnnoy_0.0.18      plyr_1.8.6           
# [97] magrittr_2.0.1        R6_2.5.0              IRanges_2.22.2        generics_0.1.0       
# [101] DBI_1.1.1             withr_2.4.2           pillar_1.6.0          mgcv_1.8-35          
# [105] fitdistrplus_1.1-3    survival_3.2-11       abind_1.4-5           tibble_3.1.1         
# [109] future.apply_1.7.0    crayon_1.4.1          KernSmooth_2.23-18    utf8_1.2.1           
# [113] BiocFileCache_1.12.1  spatstat.geom_2.1-0   plotly_4.9.3          progress_1.2.2       
# [117] grid_4.0.2            data.table_1.14.0     blob_1.2.1            digest_0.6.27        
# [121] xtable_1.8-4          tidyr_1.1.3           httpuv_1.6.0          openssl_1.4.3        
# [125] stats4_4.0.2          munsell_0.5.0         viridisLite_0.4.0     askpass_1.1 


