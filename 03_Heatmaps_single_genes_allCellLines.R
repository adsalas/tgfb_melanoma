
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
