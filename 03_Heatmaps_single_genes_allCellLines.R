
# Representative heatmaps with single genes per cell line and treatment ----
# Load the libraries
library(DESeq2)
library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)

# Load the pre-computed outputs
load("path/to/folder/where/file/is/stored/outputs.RData")

# Define the set of genes of interest
genes.oi <- c("ACTN4", 
              "ANK2",
              "APLP1",
              "APP",  
              "ATF3",  
              "ATP2A3",
              "BAP1",
              "BCL2L11",
              "CD24", 
              "DAPK3",
              "DDIT4",   
              "DYRK2", 
              "EPC2", 
              "GADD45B",
              "HTATIP2",
              "ITGA1",
              "JADE1", 
              "KDELR1",
              "MEF2A",
              "MEF2D",
              "MMP2", 
              "NOTCH2",
              "NPR2",
              "PMP22",   
              "PPP1R13L",
              "PPP2R1B",
              "PRKD1", 
              "RAPGEF2",
              "RTKN", 
              "RXRA", 
              "SCN2A", 
              "SGPL1",
              "SHC4", 
              "SMURF2",
              "TCIRG1",
              "TGFB1",
              "TIAM1", 
              "TRAF1", 
              "TRIM2",
              "TRIM24",
              "TRPV1",
              "TWIST1",  
              "UBE4B", 
              "ZNF274"
)

# Double check that ALL the genes are present in the table
length(genes.oi)
table(genes.oi %in% ext.gene.names$external_gene_name)
# Get the Ensembl IDs
genes.oi.enseml.ordered <- ext.gene.names
genes.oi.enseml.ordered <- genes.oi.enseml.ordered[genes.oi.enseml.ordered$external_gene_name %in% genes.oi, ] 
rownames(genes.oi.enseml.ordered) <- genes.oi.enseml.ordered$external_gene_name
genes.oi.enseml <- genes.oi.enseml.ordered[genes.oi, ] %>% pull(ensembl_gene_id)
genes.oi.enseml

# Set the coldata
# Set the order of the treatments/samples in the coldata df
coldata

coldata$treatment.order <- 4
coldata[coldata$treatment == "DMSO", ]$treatment.order <- 1
coldata[coldata$treatment == "TGFb", ]$treatment.order <- 2
coldata[coldata$treatment == "MEKi", ]$treatment.order <- 3
coldata

coldata <- arrange(coldata, treatment.order)
coldata

# Retrieve the rlog normalized values for the representative heatmaps
df_norm <- assay(rld.dds)
df_norm <- as.data.frame(df_norm)
colnames(df_norm)
coldata$Sample
all(colnames(df_norm) == coldata$Sample)
# Re-order the colnames of the df_norm
df_norm <- df_norm[ , coldata$Sample]
colnames(df_norm)
coldata$Sample
all(colnames(df_norm) == coldata$Sample)
colnames(df_norm) <- coldata$sample_full

# Create the dataframe where to store the average expression
table.mean.ALL <- data.frame(matrix(0, nrow = length(rownames(df_norm)), ncol = length(unique(coldata$cell_line))*length(unique(coldata$treatment))))
rownames(table.mean.ALL) <- rownames(df_norm)

# Compute the average and populate the table
cell_line_name <- unique(coldata$cell_line)
treatment_type <- unique(coldata$treatment)

treatment_type_full <- rep(treatment_type, length(cell_line_name))
cell_line_name_full <- rep(cell_line_name, length(treatment_type))
colnames(table.mean.ALL) <- paste0(treatment_type_full, "_", cell_line_name_full, "_", "mean")
# Define a counter for the loop computations
z <- 1
# Populate the table with the results
for (i in 1:length(cell_line_name)) {
  for (j in 1:length(treatment_type)) {
    x <- filter(coldata, cell_line == cell_line_name[i] & treatment == treatment_type[j]) %>% pull(sample_full)
    if (length(x) > 1) {
      y <- apply(df_norm[ , x], 1, mean)
    }
    else
      y <- df_norm[ , x]
    
    table.mean.ALL[ , z] <- y
    colnames(table.mean.ALL)[z] <- paste0(treatment_type[j], "_", cell_line_name[i], "_", "mean")
    z <- z + 1
  }
}

table(rowSums(table.mean.ALL) == 0)

# Inspect the dataframe
head(table.mean.ALL)

# Filtering by the geneset of interest
df <- table.mean.ALL
# Check if all genes are present in the data
table(rownames(df) %in% genes.oi.enseml)
genes.oi.enseml[!genes.oi.enseml %in% rownames(df)]
df_zsc <- df[genes.oi.enseml, ]
dim(df_zsc)
# Create a dataframe with the structure required for storing the data for the heatmap representation
df_zsc <- t(as.data.frame(t(df_zsc)))

# Define a variable containing the cell lines information
my.cell.line <- as.character(unique(coldata$cell_line))
# Create the dataframe where to store the data
data.temp.full <- data.frame(matrix(0, nrow = length(my.cell.line), ncol = length(unique(coldata$treatment))))
rownames(data.temp.full) <- my.cell.line
colnames(data.temp.full) <- unique(coldata$treatment)

# Arrange the data and export the plots 
for (i in 1:length(genes.oi.enseml)) {
  
  for (j in my.cell.line) {
    
    data.temp <- t(df_zsc[genes.oi.enseml[i], grep(j, colnames(df_zsc))])
    rownames(data.temp) <- unique(sapply(strsplit(colnames(data.temp),"\\_"), function(x) x[2]))
    colnames(data.temp) <- unique(sapply(strsplit(colnames(data.temp),"\\_"), function(x) x[1]))
    data.temp <- t(as.data.frame(scale(t(data.temp))))
    
    data.temp.full[j, ] <- data.temp
    
  }
  
  # Fix the order of the columns
  data.temp.full <- data.temp.full[ , c("DMSO", "TGFb", "MEKi", "TGFb + MEKi")]
  # Fix the order of the rows
  data.temp.full <- data.temp.full[c("M130903","M010817",  "M170117","M161201", "MM150543", "M130830"), ]
  
  # Define the annotations
  treatment.ann <- colnames(data.temp.full)
  # Color space for the annotation
  my.col <- brewer.pal(length(unique(treatment.ann)), "Set3")
  
  cell.ann <- rownames(data.temp.full)
  cell.ann <- factor(cell.ann, levels = c("M130903","M010817",  "M170117","M161201", "MM150543", "M130830"))
  # Color space for the annotation
  my.col2 <- brewer.pal(6, "Set2")
  
  # Heatmap annotations ----
  # Column annotations
  names(my.col) <- treatment.ann
  
  column_ha <- HeatmapAnnotation(treatment = treatment.ann,
                                 col = list(treatment = my.col))
  
  # Row annotations
  names(my.col2) <- cell.ann
  
  row_ha = rowAnnotation(cell_line = cell.ann, 
                         col = list(cell_line = my.col2))
  
  # Plot the heatmaps
  set.seed(1234)
  ht1 <- Heatmap(data.temp.full, 
                 show_row_names = T,
                 show_column_names = F, 
                 cluster_rows = F,
                 cluster_columns = F, 
                 column_title = genes.oi[i],
                 name = "zscore", 
                 na_col = "grey",
                 top_annotation = column_ha,
                 right_annotation = row_ha
  )
  
  # Export the plots
  somePDFPath = paste0("path/to/folder/where/file/will/be/stored/", "Heatmap_", genes.oi[i], "_CellLines_by_treatment.pdf")
  pdf(file=somePDFPath, height = 2.5, width = 6)
  print(ht1)
  dev.off()
}

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
#   [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets  methods  
# [10] base     
# 
# other attached packages:
#   [1] RColorBrewer_1.1-2          ComplexHeatmap_2.4.3        dplyr_1.0.5                
# [4] DESeq2_1.28.1               SummarizedExperiment_1.18.2 DelayedArray_0.14.1        
# [7] matrixStats_0.58.0          Biobase_2.48.0              GenomicRanges_1.40.0       
# [10] GenomeInfoDb_1.24.2         IRanges_2.22.2              S4Vectors_0.26.1           
# [13] BiocGenerics_0.34.0        
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.6             locfit_1.5-9.4         circlize_0.4.13        lattice_0.20-41       
# [5] png_0.1-7              assertthat_0.2.1       utf8_1.2.1             R6_2.5.0              
# [9] RSQLite_2.2.7          ggplot2_3.3.3          pillar_1.6.0           GlobalOptions_0.1.2   
# [13] zlibbioc_1.34.0        rlang_0.4.10           annotate_1.66.0        blob_1.2.1            
# [17] GetoptLong_1.0.5       Matrix_1.3-2           splines_4.0.2          BiocParallel_1.22.0   
# [21] geneplotter_1.66.0     RCurl_1.98-1.3         bit_4.0.4              munsell_0.5.0         
# [25] compiler_4.0.2         pkgconfig_2.0.3        shape_1.4.6            tidyselect_1.1.0      
# [29] tibble_3.1.1           GenomeInfoDbData_1.2.3 XML_3.99-0.6           fansi_0.4.2           
# [33] crayon_1.4.1           bitops_1.0-7           xtable_1.8-4           gtable_0.3.0          
# [37] lifecycle_1.0.0        DBI_1.1.1              magrittr_2.0.1         scales_1.1.1          
# [41] cachem_1.0.4           XVector_0.28.0         genefilter_1.70.0      ellipsis_0.3.1        
# [45] vctrs_0.3.7            generics_0.1.0         rjson_0.2.20           tools_4.0.2           
# [49] bit64_4.0.5            glue_1.4.2             purrr_0.3.4            fastmap_1.1.0         
# [53] survival_3.2-11        clue_0.3-60            AnnotationDbi_1.50.3   colorspace_2.0-0      
# [57] cluster_2.1.2          memoise_2.0.0  




