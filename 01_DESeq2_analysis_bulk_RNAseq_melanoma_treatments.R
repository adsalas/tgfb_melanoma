#####  Analysis bulk RNAseq melanoma cell lines treated with MEKi, TGFb and MEKi + TGFb  #######

# Load the libraries ----
library(DESeq2)
library(dplyr)
library(biomaRt)

# Import and format the data ----
# Load the count data table deposited in GEO
counttablefull <- read.table("path/to/folder/where/file/is/stored/raw.counts.txt")

# Load the metadata (available in Github under the "files" folder)
metadata <- readRDS("path/to/folder/where/file/is/stored/metadata.RDS") 
metadata
# Add the full name of the sample in the metadata
metadata$sample_full <- paste0(metadata$treatment, "_", metadata$replicate, ".", metadata$cell_line)

# Convert the Ensembl names to gene IDs with bioMart ----
# NOTES:
#     1. The database used for querying the external gene names was updated time after 
#        the analyses were performed, therefore, we provided in the GitHub repository 
#        an R object that contains the information retrieved at the time.
#     2. If one would re-run the query again the results will slightly change (there are 
#        more genes for which external names are available), however, that will not affect
#        the main results but some of the downstream analyses/representations in which genes
#        with no external names reported were filtered out. The code used for the query is 
#        included below but has been commented out
#     3. The queried external gene names cam be load using the following line of code after
#        providing the path to the R object "ext_gene_names.RDS"

# Import the external gene names queried from Ensembl (available in Github under the "files" folder)
ext.gene.names <- readRDS("path/to/folder/where/file/is/stored/ext_gene_names.RDS") 

# # Define the database for the query
# mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# # Define the set of genes to perform the query for
# mygenesEnsem <- counttablefull %>% rownames()
# mygenesEnsem
# length(mygenesEnsem)
# 
# # Make the query of the required information
# ext.gene.names <- getBM(attributes = c('ensembl_gene_id',
#                                        'external_gene_name'),
#                         filters = 'ensembl_gene_id', 
#                         values = mygenesEnsem,
#                         mart = mart)
# 
# dim(ext.gene.names)
# length(unique(ext.gene.names$external_gene_name))

# Check that the genes are matched in order between the query and the counts table
all(rownames(ext.gene.names) == rownames(counttablefull))
head(counttablefull)
head(ext.gene.names)

# Running the DESeq2 analysis on the subset of samples that have replicates ----
# Assign the metadata as coldata
coldata <- metadata
# Set the coldata rows in the same order as the colnames of the count table
coldata <- coldata[colnames(counttablefull), ]
# Check the cell lines that have the 3 replicates
table(coldata$cell_line, coldata$replicate)
# Filter the coldata table
cell.line <- coldata %>% filter(replicate == "REP3") %>% pull(cell_line) %>% unique() 
samples.keep <- coldata %>% filter(cell_line %in% cell.line) %>% rownames()
length(samples.keep)
coldata.2 <- coldata[samples.keep, ]
coldata.2$treatment <- factor(coldata.2$treatment, levels=c("DMSO", "TGFb", "MEKi", "TGFb + MEKi"))
# Filter the count table
counttable <- counttablefull[, rownames(coldata.2)]
counttable[1:10, 1:5]

# Create the DESeq object
dds2 <- DESeqDataSetFromMatrix(countData = counttable,
                               colData = coldata.2,
                               design = ~ cell_line + treatment + cell_line:treatment)

dds2
mcols(dds2)

# Pre-filtering of low count genes
keep <- rowSums(counts(dds2)) >= 10
dds2 <- dds2[keep,]

dds2

# Differential expression analysis
dds2 <- DESeq(dds2)
res <- results(dds2)
res

# Perform data normalization for the PCA representation
rld.dds2 <- rlog(dds2, blind=FALSE)

# Check the separation of the samples in a PCA
plotPCA(rld.dds2, intgroup="treatment")
plotPCA(rld.dds2, intgroup="cell_line")

# Assign the model matrix to an object to store the data
mod.mat <- model.matrix(design(dds2), colData(dds2))

# Define coefficient vectors for each condition per cell line
M170117_DMSO <- colMeans(mod.mat[dds2$cell_line == "M170117" & dds2$treatment == "DMSO", ])
M170117_TGFb <- colMeans(mod.mat[dds2$cell_line == "M170117" & dds2$treatment == "TGFb", ])
M170117_MEKi <- colMeans(mod.mat[dds2$cell_line == "M170117" & dds2$treatment == "MEKi", ])
M170117_TGFb.MEKi <- colMeans(mod.mat[dds2$cell_line == "M170117" & dds2$treatment == "TGFb + MEKi", ])

M130830_DMSO <- colMeans(mod.mat[dds2$cell_line == "M130830" & dds2$treatment == "DMSO", ])
M130830_TGFb <- colMeans(mod.mat[dds2$cell_line == "M130830" & dds2$treatment == "TGFb", ])
M130830_MEKi <- colMeans(mod.mat[dds2$cell_line == "M130830" & dds2$treatment == "MEKi", ])
M130830_TGFb.MEKi <- colMeans(mod.mat[dds2$cell_line == "M130830" & dds2$treatment == "TGFb + MEKi", ])

# Define any contrasts of interest and obtain results for each pairwise contrast
res_TGFb_vs_DMSO.M170117 <- results(dds2, contrast = M170117_TGFb - M170117_DMSO) # Effect of TGFb
res_MEKi_vs_DMSO.M170117 <- results(dds2, contrast = M170117_MEKi - M170117_DMSO) # Effect of MEKi
res_TGFb.MEKi_vs_DMSO.M170117 <- results(dds2, contrast = M170117_TGFb.MEKi - M170117_DMSO) # Effect of combining TGFb and MEKi

res_TGFb_vs_DMSO.M130830 <- results(dds2, contrast = M130830_TGFb - M130830_DMSO) # Effect of TGFb
res_MEKi_vs_DMSO.M130830 <- results(dds2, contrast = M130830_MEKi - M130830_DMSO) # Effect of MEKi
res_TGFb.MEKi_vs_DMSO.M130830 <- results(dds2, contrast = M130830_TGFb.MEKi - M130830_DMSO) # Effect of combining TGFb and MEKi

# Filter the results based on a padj threshold
res_TGFb_vs_DMSO.M170117_05padj <- filter(as.data.frame(res_TGFb_vs_DMSO.M170117), padj < 0.05)
res_MEKi_vs_DMSO.M170117_05padj <- filter(as.data.frame(res_MEKi_vs_DMSO.M170117), padj < 0.05)
res_TGFb.MEKi_vs_DMSO.M170117_05padj <- filter(as.data.frame(res_TGFb.MEKi_vs_DMSO.M170117), padj < 0.05) 

res_TGFb_vs_DMSO.M130830_05padj <- filter(as.data.frame(res_TGFb_vs_DMSO.M130830), padj < 0.05)
res_MEKi_vs_DMSO.M130830_05padj <- filter(as.data.frame(res_MEKi_vs_DMSO.M130830), padj < 0.05)
res_TGFb.MEKi_vs_DMSO.M130830_05padj <- filter(as.data.frame(res_TGFb.MEKi_vs_DMSO.M130830), padj < 0.05) 


# Running the DESeq2 analysis on all samples (including the ones with 1N) for the representative heatmaps ----
# Assign the metadata as coldata
coldata <- metadata
# Set the coldata rows in the same order as the colnames of the count table
coldata <- coldata[colnames(counttablefull), ]
# Create the DESeq object
dds <- DESeqDataSetFromMatrix(countData = counttablefull,
                              colData = coldata,
                              design = ~ cell_line + treatment + cell_line:treatment)

dds
mcols(dds)

# Pre-filtering of low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

# Differential expression analysis
dds <- DESeq(dds)
res <- results(dds)
res

# Perform data normalization for representation
rld.dds <- rlog(dds, blind=FALSE)
# NOTE:
#       This normalized result is used for the representative heatmaps made for 
#       the expression of single genes displaying all samples

# Save the outputs/results required for further analyses
save(dds, dds2, rld.dds, coldata, coldata.2, res_TGFb_vs_DMSO.M170117_05padj, 
     res_MEKi_vs_DMSO.M170117_05padj,res_TGFb.MEKi_vs_DMSO.M170117_05padj, 
     res_TGFb_vs_DMSO.M130830_05padj, res_MEKi_vs_DMSO.M130830_05padj,
     res_TGFb.MEKi_vs_DMSO.M130830_05padj, ext.gene.names, 
     file = "path/to/folder/where/file/will/be/stored/outputs.RData")

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
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] biomaRt_2.44.4              dplyr_1.0.5                 DESeq2_1.28.1              
# [4] SummarizedExperiment_1.18.2 DelayedArray_0.14.1         matrixStats_0.58.0         
# [7] Biobase_2.48.0              GenomicRanges_1.40.0        GenomeInfoDb_1.24.2        
# [10] IRanges_2.22.2              S4Vectors_0.26.1            BiocGenerics_0.34.0        
# 
# loaded via a namespace (and not attached):
#   [1] fgsea_1.14.0           colorspace_2.0-0       ellipsis_0.3.1         ggridges_0.5.3        
# [5] circlize_0.4.13        qvalue_2.20.0          XVector_0.28.0         GlobalOptions_0.1.2   
# [9] farver_2.1.0           urltools_1.7.3         graphlayouts_0.7.1     ggrepel_0.9.1         
# [13] bit64_4.0.5            AnnotationDbi_1.50.3   fansi_0.4.2            scatterpie_0.1.7      
# [17] xml2_1.3.2             splines_4.0.2          cachem_1.0.4           GOSemSim_2.14.2       
# [21] geneplotter_1.66.0     polyclip_1.10-0        jsonlite_1.7.2         annotate_1.66.0       
# [25] dbplyr_2.1.1           GO.db_3.11.4           ggforce_0.3.3          BiocManager_1.30.12   
# [29] compiler_4.0.2         httr_1.4.2             rvcheck_0.2.1          assertthat_0.2.1      
# [33] Matrix_1.3-2           fastmap_1.1.0          tweenr_1.0.2           prettyunits_1.1.1     
# [37] tools_4.0.2            igraph_1.2.6           gtable_0.3.0           glue_1.4.2            
# [41] GenomeInfoDbData_1.2.3 reshape2_1.4.4         DO.db_2.9              rappdirs_0.3.3        
# [45] fastmatch_1.1-3        Rcpp_1.0.6             enrichplot_1.8.1       cellranger_1.1.0      
# [49] vctrs_0.3.7            ggraph_2.0.5           stringr_1.4.0          lifecycle_1.0.0       
# [53] clusterProfiler_3.16.1 XML_3.99-0.6           DOSE_3.14.0            europepmc_0.4.1       
# [57] MASS_7.3-53.1          zlibbioc_1.34.0        scales_1.1.1           tidygraph_1.2.0       
# [61] hms_1.0.0              RColorBrewer_1.1-2     curl_4.3               memoise_2.0.0         
# [65] gridExtra_2.3          ggplot2_3.3.3          downloader_0.4         ggfun_0.0.9           
# [69] yulab.utils_0.0.4      triebeard_0.3.0        stringi_1.5.3          RSQLite_2.2.7         
# [73] genefilter_1.70.0      BiocParallel_1.22.0    shape_1.4.6            rlang_0.4.10          
# [77] pkgconfig_2.0.3        bitops_1.0-7           lattice_0.20-41        purrr_0.3.4           
# [81] labeling_0.4.2         cowplot_1.1.1          bit_4.0.4              tidyselect_1.1.0      
# [85] plyr_1.8.6             magrittr_2.0.1         R6_2.5.0               generics_0.1.0        
# [89] DBI_1.1.1              pillar_1.6.0           survival_3.2-11        RCurl_1.98-1.3        
# [93] tibble_3.1.1           crayon_1.4.1           utf8_1.2.1             BiocFileCache_1.12.1  
# [97] viridis_0.6.0          progress_1.2.2         locfit_1.5-9.4         grid_4.0.2            
# [101] readxl_1.3.1           data.table_1.14.0      blob_1.2.1             digest_0.6.27         
# [105] xtable_1.8-4           tidyr_1.1.3            gridGraphics_0.5-1     openssl_1.4.3         
# [109] munsell_0.5.0          viridisLite_0.4.0      ggplotify_0.1.0        askpass_1.1




