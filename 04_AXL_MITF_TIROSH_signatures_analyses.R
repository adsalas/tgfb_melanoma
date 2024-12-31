######   GENE SET ENRICHMENT ANALYSES WITH THE 3 SIGNATURES: AXL, MITF, TIROSH RESISTANCE AT THE SAME TIME ###########

# Load the libraries ----
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)
library(DESeq2)

# Define the directories ----
# pathScripts = "path/to/scripts/"
pathData = "path/to/data/"
pathOutput = "path/to/output/"
# Set the output directory as working directory (where outputs are to be saved)
wd = pathOutput
setwd(wd)

# Load the data ----
# Load the pre-computed outputs
load("path/to/folder/where/file/is/stored/outputs.RData") 

# Load and format the gene signatures (from Tirosh et al., 2016. Science. DOI: 10.1126/science.aad0501)
# AXL signature
AXL.signature <- readxl::read_xlsx(paste0(pathData, "aad0501_table_s8.xlsx"), skip = 4, col_names = FALSE)
AXL.signature
AXL.signature$Signature <- "AXL"
colnames(AXL.signature) <- c("Gene", "Signature")
AXL.signature <- AXL.signature[ , c(2,1)]
# MITF signature
MITF.signature <- readxl::read_xlsx(paste0(pathData, "aad0501_table_s7.xlsx"), skip = 4, col_names = FALSE)
MITF.signature
MITF.signature$Signature <- "MITF"
colnames(MITF.signature) <- c("Gene", "Signature")
MITF.signature <- MITF.signature[ , c(2,1)]
# Tirosh resistance signature
tirosh.resistance <- readxl::read_xlsx(paste0(pathData, "aad0501_table_s6.xlsx"), skip = 4, col_names = TRUE)
tirosh.resistance
tirosh.resistance <- tirosh.resistance[ , 1]
tirosh.resistance
tirosh.resistance$Signature <- "Tirosh Resistance"
colnames(tirosh.resistance) <- c("Gene", "Signature")
tirosh.resistance <- tirosh.resistance[ , c(2,1)]
tirosh.resistance

# Combine all the Signatures in a dataframe to perform the analysis simultaneously
AXL.MITF.TIROSH.Signatures <- rbind(tirosh.resistance, AXL.signature, MITF.signature)
AXL.MITF.TIROSH.Signatures
unique(AXL.MITF.TIROSH.Signatures$Signature)

# Performing Gene set enrichment analysis ----
# Define the genesets I want to used for the analyses
# Define the set of genes to test for
################### M170117
######## DEGs TGFb vs DMSO
#   Check which genes don't have an external name and remove them
TGFb.M170117 <- ext.gene.names[rownames(res_TGFb_vs_DMSO.M170117_05padj), ]
TGFb.M170117 <- filter(TGFb.M170117, external_gene_name != "")
dim(TGFb.M170117)
#   Create the named vector with the log2FC values for the selected set of genes
TGFb.M170117$log2FC <- res_TGFb_vs_DMSO.M170117_05padj[rownames(TGFb.M170117), "log2FoldChange"]
TGFb.M170117.genes <- TGFb.M170117$external_gene_name
TGFb.M170117 <- TGFb.M170117$log2FC
names(TGFb.M170117) <- TGFb.M170117.genes
#    Order the vector in decreasing order based on log2FC
TGFb.M170117 <- TGFb.M170117[order(TGFb.M170117, decreasing = T)]
TGFb.M170117
tail(TGFb.M170117)
######## DEGs MEKi vs DMSO
#   Check which genes don't have an external name and remove them
MEKi.M170117 <- ext.gene.names[rownames(res_MEKi_vs_DMSO.M170117_05padj), ]
MEKi.M170117 <- filter(MEKi.M170117, external_gene_name != "")
dim(MEKi.M170117)
#   Create the named vector with the log2FC values for the selected set of genes
MEKi.M170117$log2FC <- res_MEKi_vs_DMSO.M170117_05padj[rownames(MEKi.M170117), "log2FoldChange"]
MEKi.M170117.genes <- MEKi.M170117$external_gene_name
MEKi.M170117 <- MEKi.M170117$log2FC
names(MEKi.M170117) <- MEKi.M170117.genes
#    Order the vector in decreasing order based on log2FC
MEKi.M170117 <- MEKi.M170117[order(MEKi.M170117, decreasing = T)]
MEKi.M170117
tail(MEKi.M170117)
######## DEGs TGFb.MEKi vs DMSO
#   Check which genes don't have an external name and remove them
TGFb.MEKi.M170117 <- ext.gene.names[rownames(res_TGFb.MEKi_vs_DMSO.M170117_05padj), ]
TGFb.MEKi.M170117 <- filter(TGFb.MEKi.M170117, external_gene_name != "")
dim(TGFb.MEKi.M170117)
#   Create the named vector with the log2FC values for the selected set of genes
TGFb.MEKi.M170117$log2FC <- res_TGFb.MEKi_vs_DMSO.M170117_05padj[rownames(TGFb.MEKi.M170117), "log2FoldChange"]
TGFb.MEKi.M170117.genes <- TGFb.MEKi.M170117$external_gene_name
TGFb.MEKi.M170117 <- TGFb.MEKi.M170117$log2FC
names(TGFb.MEKi.M170117) <- TGFb.MEKi.M170117.genes
#    Order the vector in decreasing order based on log2FC
TGFb.MEKi.M170117 <- TGFb.MEKi.M170117[order(TGFb.MEKi.M170117, decreasing = T)]
TGFb.MEKi.M170117
tail(TGFb.MEKi.M170117)

################### M130830
######## DEGs TGFb vs DMSO
#   Check which genes don't have an external name and remove them
TGFb.M130830 <- ext.gene.names[rownames(res_TGFb_vs_DMSO.M130830_05padj), ]
TGFb.M130830 <- filter(TGFb.M130830, external_gene_name != "")
dim(TGFb.M130830)
#   Create the named vector with the log2FC values for the selected set of genes
TGFb.M130830$log2FC <- res_TGFb_vs_DMSO.M130830_05padj[rownames(TGFb.M130830), "log2FoldChange"]
TGFb.M130830.genes <- TGFb.M130830$external_gene_name
TGFb.M130830 <- TGFb.M130830$log2FC
names(TGFb.M130830) <- TGFb.M130830.genes
#    Order the vector in decreasing order based on log2FC
TGFb.M130830 <- TGFb.M130830[order(TGFb.M130830, decreasing = T)]
TGFb.M130830
tail(TGFb.M130830)
######## DEGs MEKi vs DMSO
#   Check which genes don't have an external name and remove them
MEKi.M130830 <- ext.gene.names[rownames(res_MEKi_vs_DMSO.M130830_05padj), ]
MEKi.M130830 <- filter(MEKi.M130830, external_gene_name != "")
dim(MEKi.M130830)
#   Create the named vector with the log2FC values for the selected set of genes
MEKi.M130830$log2FC <- res_MEKi_vs_DMSO.M130830_05padj[rownames(MEKi.M130830), "log2FoldChange"]
MEKi.M130830.genes <- MEKi.M130830$external_gene_name
MEKi.M130830 <- MEKi.M130830$log2FC
names(MEKi.M130830) <- MEKi.M130830.genes
#    Order the vector in decreasing order based on log2FC
MEKi.M130830 <- MEKi.M130830[order(MEKi.M130830, decreasing = T)]
MEKi.M130830
tail(MEKi.M130830)
######## DEGs TGFb.MEKi vs DMSO
#   Check which genes don't have an external name and remove them
TGFb.MEKi.M130830 <- ext.gene.names[rownames(res_TGFb.MEKi_vs_DMSO.M130830_05padj), ]
TGFb.MEKi.M130830 <- filter(TGFb.MEKi.M130830, external_gene_name != "")
dim(TGFb.MEKi.M130830)
#   Create the named vector with the log2FC values for the selected set of genes
TGFb.MEKi.M130830$log2FC <- res_TGFb.MEKi_vs_DMSO.M130830_05padj[rownames(TGFb.MEKi.M130830), "log2FoldChange"]
TGFb.MEKi.M130830.genes <- TGFb.MEKi.M130830$external_gene_name
TGFb.MEKi.M130830 <- TGFb.MEKi.M130830$log2FC
names(TGFb.MEKi.M130830) <- TGFb.MEKi.M130830.genes
#    Order the vector in decreasing order based on log2FC
TGFb.MEKi.M130830 <- TGFb.MEKi.M130830[order(TGFb.MEKi.M130830, decreasing = T)]
TGFb.MEKi.M130830
tail(TGFb.MEKi.M130830)

# Running the GSE analysis with the defined signatures ----
##########################      M170117     ####################################
# TGFb
gsea.TGFb.M170117.TIROSH.AXL.MITF <- GSEA(TGFb.M170117,
                                      pvalueCutoff = 1,
                                      TERM2GENE = AXL.MITF.TIROSH.Signatures,
                                      seed = TRUE)

# MEKi
gsea.MEKi.M170117.TIROSH.AXL.MITF <- GSEA(MEKi.M170117,
                                               pvalueCutoff = 1,
                                               TERM2GENE = AXL.MITF.TIROSH.Signatures,
                                               seed = TRUE)

# TGFb.MEKi
gsea.TGFb.MEKi.M170117.TIROSH.AXL.MITF <- GSEA(TGFb.MEKi.M170117,
                                          pvalueCutoff = 1,
                                          TERM2GENE = AXL.MITF.TIROSH.Signatures,
                                          seed = TRUE)

##########################      M130830     ####################################
# TGFb
gsea.TGFb.M130830.TIROSH.AXL.MITF <- GSEA(TGFb.M130830,
                                          pvalueCutoff = 1,
                                          TERM2GENE = AXL.MITF.TIROSH.Signatures,
                                          seed = TRUE)

# MEKi
gsea.MEKi.M130830.TIROSH.AXL.MITF <- GSEA(MEKi.M130830,
                                          pvalueCutoff = 1,
                                          TERM2GENE = AXL.MITF.TIROSH.Signatures,
                                          seed = TRUE)

# TGFb.MEKi
gsea.TGFb.MEKi.M130830.TIROSH.AXL.MITF <- GSEA(TGFb.MEKi.M130830,
                                               pvalueCutoff = 1,
                                               TERM2GENE = AXL.MITF.TIROSH.Signatures,
                                               seed = TRUE)


# Making the plots ----
# Import the modified version of the gseaplot2 function, including a dashed line at the zero value as reference
source("~/Documents/scripts/gSeaplot2_modif.R")
##########################      M170117     ####################################
# TGFb
# somePath = "GSEA_AXL_MITF_TiroshResistance_signature_M170117_TGFb.pdf"
# pdf(file=somePath, height = 3, width = 5)
gseaplot2.test(gsea.TGFb.M170117.TIROSH.AXL.MITF, 
               geneSetID = 1:3, 
               subplots = 1:2, 
               color = c("orange", "navy", "lightblue"),
               pvalue_table = F, color_zero_line = "black", size_zero_line = 1)
# dev.off()

# MEKi
# somePath = "GSEA_AXL_MITF_TiroshResistance_signature_M170117_MEKi.pdf"
# pdf(file=somePath, height = 3, width = 5)
gseaplot2.test(gsea.MEKi.M170117.TIROSH.AXL.MITF, 
               geneSetID = 1:3, 
               subplots = 1:2, 
               color = c("orange", "navy", "lightblue"),
               pvalue_table = F, color_zero_line = "black", size_zero_line = 1)
# dev.off()

# TGFb + MEKi
# somePath = "GSEA_AXL_MITF_TiroshResistance_signature_M170117_TGFb.MEKi.pdf"
# pdf(file=somePath, height = 3, width = 5)
gseaplot2.test(gsea.TGFb.MEKi.M170117.TIROSH.AXL.MITF, 
               geneSetID = 1:3, 
               subplots = 1:2, 
               color = c("orange", "navy", "lightblue"),
               pvalue_table = F, color_zero_line = "black", size_zero_line = 1)
# dev.off()

##########################      M130830     ####################################
# TGFb
# somePath = "GSEA_AXL_MITF_TiroshResistance_signature_M130830_TGFb.pdf"
# pdf(file=somePath, height = 3, width = 5)
gseaplot2.test(gsea.TGFb.M130830.TIROSH.AXL.MITF, 
               geneSetID = 1:3, 
               subplots = 1:2, 
               color = c("orange", "navy", "lightblue"),
               pvalue_table = F, color_zero_line = "black", size_zero_line = 1)
# dev.off()

# MEKi
# somePath = "GSEA_AXL_MITF_TiroshResistance_signature_M130830_MEKi.pdf"
# pdf(file=somePath, height = 3, width = 5)
gseaplot2.test(gsea.MEKi.M130830.TIROSH.AXL.MITF, 
               geneSetID = 1:3, 
               subplots = 1:2, 
               color = c("orange", "navy", "lightblue"),
               pvalue_table = F, color_zero_line = "black", size_zero_line = 1)
# dev.off()

# TGFb + MEKi
# somePath = "GSEA_AXL_MITF_TiroshResistance_signature_M130830_TGFb.MEKi.pdf"
# pdf(file=somePath, height = 3, width = 5)
gseaplot2.test(gsea.TGFb.MEKi.M130830.TIROSH.AXL.MITF, 
               geneSetID = 1:3, 
               subplots = 1:2, 
               color = c("orange", "navy", "lightblue"),
               pvalue_table = F, color_zero_line = "black", size_zero_line = 1)
# dev.off()

# Export the results as a combined table
# Make a data frame with all the results indicating the condition and the cell line
###  M170117  #####
df1 <- gsea.TGFb.M170117.TIROSH.AXL.MITF@result
df1$Cell_Line <- "M170117"
df1$Condition <- "TGFb"
df2 <- gsea.MEKi.M170117.TIROSH.AXL.MITF@result
df2$Cell_Line <- "M170117"
df2$Condition <- "MEKi"
df3 <- gsea.TGFb.MEKi.M170117.TIROSH.AXL.MITF@result
df3$Cell_Line <- "M170117"
df3$Condition <- "TGFb_MEKi"
###  M130830  #####
df4 <- gsea.TGFb.M130830.TIROSH.AXL.MITF@result
df4$Cell_Line <- "M130830"
df4$Condition <- "TGFb"
df5 <- gsea.MEKi.M130830.TIROSH.AXL.MITF@result
df5$Cell_Line <- "M130830"
df5$Condition <- "MEKi"
df6 <- gsea.TGFb.MEKi.M130830.TIROSH.AXL.MITF@result
df6$Cell_Line <- "M130830"
df6$Condition <- "TGFb_MEKi"
# Combine the data frames
combined.df <- rbind(df1, 
                     df2,
                     df3,
                     df4,
                     df5,
                     df6)

# Filter only the required data
colnames(combined.df)
combined.df <- combined.df[ , c("Cell_Line", "Condition", "Description","enrichmentScore","NES","p.adjust")]
rownames(combined.df) <- NULL
combined.df
# Export the results
write.csv(combined.df, "path/to/folder/where/file/is/stored/Results_GSEA_AXL_MITF_TiroshResistance_signatures.csv")


# Heatmap representations for the signature gene sets in the M170117 and M130830 bulk RNAseq ----
# Retrieve the Ensembl IDs for the genes in the distinct signatures
AXL.Ensembl <- ext.gene.names[ext.gene.names$external_gene_name %in% AXL.signature$Gene, ] %>% rownames()
length(AXL.Ensembl)

MITF.Ensembl <- ext.gene.names[ext.gene.names$external_gene_name %in% MITF.signature$Gene, ] %>% rownames()
length(MITF.Ensembl)

tirosh.resistance.Ensembl <- ext.gene.names[ext.gene.names$external_gene_name %in% tirosh.resistance$Gene, ] %>% rownames()
length(tirosh.resistance.Ensembl)
# Create a vector that contains the unique genes present in the signatures
genes.signatures <- c(AXL.Ensembl, MITF.Ensembl, tirosh.resistance.Ensembl) %>% unique()

############## M170117   ##################
# Retrieve the set of DEGs for each treatment vs DMSO
genes.DEG <- c(rownames(res_TGFb_vs_DMSO.M170117_05padj), 
               rownames(res_MEKi_vs_DMSO.M170117_05padj),
               rownames(res_TGFb.MEKi_vs_DMSO.M170117_05padj))
genes.DEG <- unique(genes.DEG)
length(genes.DEG)


# Subset the coldata by cell line of interest
coldata.M170117 <- filter(coldata.2, cell_line == "M170117")
# Set the order of the treatments/samples in the coldata df
coldata.M170117$treatment.order <- 4
coldata.M170117[coldata.M170117$treatment == "DMSO", ]$treatment.order <- 1
coldata.M170117[coldata.M170117$treatment == "TGFb", ]$treatment.order <- 2
coldata.M170117[coldata.M170117$treatment == "MEKi", ]$treatment.order <- 3
coldata.M170117

coldata.M170117 <- arrange(coldata.M170117, treatment.order)
coldata.M170117

# Retrieve the normalized expression values
df_norm <- counts(dds2, normalized= T)[, ]
all(colnames(df_norm) == rownames(coldata.2))
colnames(df_norm) <- coldata.2$sample_full
dim(df_norm)
# Subset the Cell line of interest
df_norm_M170 <- df_norm[ , grep("M170117", colnames(df_norm), value = T)]
dim(df_norm_M170)
# Filter the expression matrix using the gene set of interest
df_norm_M170 <- df_norm_M170[rownames(df_norm_M170) %in% genes.DEG, ]
dim(df_norm_M170)
# Keep only the genes that are found in the signatures
df_norm_M170 <- df_norm_M170[rownames(df_norm_M170) %in% genes.signatures, ]
dim(df_norm_M170)
# Set the order of the treatments/samples in the expression matrix
df_norm_M170 <- df_norm_M170[ , coldata.M170117$sample_full]
# Scale the data
df_zsc_M170 <- t(as.data.frame((scale(t(df_norm_M170)))))
colnames(df_zsc_M170)
# Format the column names
colnames(df_zsc_M170) <- sub(".M170117", "", colnames(df_zsc_M170))
colnames(df_zsc_M170)
dim(df_zsc_M170)
# Remove genes with NA entries
df_zsc_M170 <- df_zsc_M170[!is.na(rowSums(df_zsc_M170)), ]
dim(df_zsc_M170)

# Define the annotations
treatment.ann <- sapply(strsplit(colnames(df_zsc_M170),"\\_"), function(x) x[1])
treatment.ann <- factor(treatment.ann, levels = c("DMSO", "TGFb", "MEKi", "TGFb + MEKi"))
# Defien the color for the annotations
my.col <- brewer.pal(length(unique(treatment.ann)), "Set3")

# Column annotation
column_ha <- HeatmapAnnotation(treatment = treatment.ann, 
                               col = list(treatment = c("DMSO" = my.col[1], 
                                                        "TGFb" = my.col[2], 
                                                        "MEKi" = my.col[3], 
                                                        "TGFb + MEKi" = my.col[4])))
# Create the individual heatmaps
####  AXL signature
set.seed(1234)
heatm.AXL.M170117 <- Heatmap(df_zsc_M170[rownames(df_zsc_M170) %in% AXL.Ensembl, ], 
                             show_row_names = F, 
                             show_column_names = F,
                             column_title = "M170117", 
                             name = "zscore", 
                             row_title = "AXL program",
                             column_order = colnames(df_zsc_M170),
                             top_annotation = column_ha) 


heatm.AXL.M170117

####  MITF signature
set.seed(1234)
heatm.MITF.M170117 <- Heatmap(df_zsc_M170[rownames(df_zsc_M170) %in% MITF.Ensembl, ], 
                              show_row_names = F, 
                              show_column_names = F,
                              column_title = "M170117", 
                              name = "zscore", 
                              row_title = "MITF program",
                              column_order = colnames(df_zsc_M170)
) 


heatm.MITF.M170117

####  Tirosh resistance signature
set.seed(1234)
heatm.tirosh.RES.M170117 <- Heatmap(df_zsc_M170[rownames(df_zsc_M170) %in% tirosh.resistance.Ensembl, ], 
                                    show_row_names = F, 
                                    show_column_names = F,
                                    column_title = "M170117", 
                                    name = "zscore", 
                                    row_title = "TIROSH Res. program",
                                    column_order = colnames(df_zsc_M170)
) 


heatm.tirosh.RES.M170117

# Combine and export the heatmaps
ht_list <- heatm.AXL.M170117 %v% heatm.MITF.M170117 %v% heatm.tirosh.RES.M170117
# somePDFPath = "Heatmap_DEGs_EachTreatmentvsDMSO_inSignatures_M170117.pdf"
# pdf(file=somePDFPath, height = 7, width = 4)
draw(ht_list)
# dev.off()


############## M130830   ##################
# Retrieve the set of DEGs for each treatment vs DMSO
genes.DEG <- c(rownames(res_TGFb_vs_DMSO.M130830_05padj), 
               rownames(res_MEKi_vs_DMSO.M130830_05padj),
               rownames(res_TGFb.MEKi_vs_DMSO.M130830_05padj))
genes.DEG <- unique(genes.DEG)
length(genes.DEG)

# Subset the coldata by cell line of interest
coldata.M130830 <- filter(coldata.2, cell_line == "M130830")
# Set the order of the treatments/samples in the coldata df
coldata.M130830$treatment.order <- 4
coldata.M130830[coldata.M130830$treatment == "DMSO", ]$treatment.order <- 1
coldata.M130830[coldata.M130830$treatment == "TGFb", ]$treatment.order <- 2
coldata.M130830[coldata.M130830$treatment == "MEKi", ]$treatment.order <- 3
coldata.M130830

coldata.M130830 <- arrange(coldata.M130830, treatment.order)
coldata.M130830

# Retrieve the normalized expression values
df_norm <- counts(dds2, normalized= T)[, ]
all(colnames(df_norm) == rownames(coldata.2))
colnames(df_norm) <- coldata.2$sample_full
dim(df_norm)
# Subset the Cell line of interest
df_norm_M130 <- df_norm[ , grep("M130830", colnames(df_norm), value = T)]
dim(df_norm_M130)
# Filter the expression matrix using the gene set of interest
df_norm_M130 <- df_norm_M130[rownames(df_norm_M130) %in% genes.DEG, ]
dim(df_norm_M130)
# Keep only the genes that are found in the signatures
df_norm_M130 <- df_norm_M130[rownames(df_norm_M130) %in% genes.signatures, ]
dim(df_norm_M130)
# Set the order of the treatments/samples in the expression matrix
df_norm_M130 <- df_norm_M130[ , coldata.M130830$sample_full]
# Scale the data
df_zsc_M130 <- t(as.data.frame((scale(t(df_norm_M130)))))
colnames(df_zsc_M130)
# Format the column names
colnames(df_zsc_M130) <- sub(".M130830", "", colnames(df_zsc_M130))
colnames(df_zsc_M130)
dim(df_zsc_M130)
# Remove genes with NA entries
df_zsc_M130 <- df_zsc_M130[!is.na(rowSums(df_zsc_M130)), ]
dim(df_zsc_M130)

# Define the annotations
treatment.ann <- sapply(strsplit(colnames(df_zsc_M130),"\\_"), function(x) x[1])
treatment.ann <- factor(treatment.ann, levels = c("DMSO", "TGFb", "MEKi", "TGFb + MEKi"))
# Defien the color for the annotations
my.col <- brewer.pal(length(unique(treatment.ann)), "Set3")

# Column annotation
column_ha <- HeatmapAnnotation(treatment = treatment.ann, 
                               col = list(treatment = c("DMSO" = my.col[1], 
                                                        "TGFb" = my.col[2], 
                                                        "MEKi" = my.col[3], 
                                                        "TGFb + MEKi" = my.col[4])))
# Create the individual heatmaps
####  AXL signature
set.seed(1234)
heatm.AXL.M130830 <- Heatmap(df_zsc_M130[rownames(df_zsc_M130) %in% AXL.Ensembl, ], 
                             show_row_names = F, 
                             show_column_names = F,
                             column_title = "M130830", 
                             name = "zscore", 
                             row_title = "AXL program",
                             column_order = colnames(df_zsc_M130),
                             top_annotation = column_ha) 


heatm.AXL.M130830

####  MITF signature
set.seed(1234)
heatm.MITF.M130830 <- Heatmap(df_zsc_M130[rownames(df_zsc_M130) %in% MITF.Ensembl, ], 
                              show_row_names = F, 
                              show_column_names = F,
                              column_title = "M130830", 
                              name = "zscore", 
                              row_title = "MITF program",
                              column_order = colnames(df_zsc_M130)
) 


heatm.MITF.M130830

####  Tirosh resistance signature
set.seed(1234)
heatm.tirosh.RES.M130830 <- Heatmap(df_zsc_M130[rownames(df_zsc_M130) %in% tirosh.resistance.Ensembl, ], 
                                    show_row_names = F, 
                                    show_column_names = F,
                                    column_title = "M130830", 
                                    name = "zscore", 
                                    row_title = "TIROSH Res. program",
                                    column_order = colnames(df_zsc_M130)
) 


heatm.tirosh.RES.M130830

# Combine and export the heatmaps
ht_list <- heatm.AXL.M130830 %v% heatm.MITF.M130830 %v% heatm.tirosh.RES.M130830
# somePDFPath = "Heatmap_DEGs_EachTreatmentvsDMSO_inSignatures_M130830.pdf"
# pdf(file=somePDFPath, height = 7, width = 4)
draw(ht_list)
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
#   [1] parallel  stats4    grid      stats     graphics  grDevices utils     datasets  methods  
# [10] base     
# 
# other attached packages:
#   [1] DESeq2_1.28.1               SummarizedExperiment_1.18.2 DelayedArray_0.14.1        
# [4] matrixStats_0.58.0          Biobase_2.48.0              GenomicRanges_1.40.0       
# [7] GenomeInfoDb_1.24.2         IRanges_2.22.2              S4Vectors_0.26.1           
# [10] BiocGenerics_0.34.0         ComplexHeatmap_2.4.3        RColorBrewer_1.1-2         
# [13] ggplot2_3.3.3               dplyr_1.0.5                 enrichplot_1.8.1           
# [16] clusterProfiler_3.16.1     
# 
# loaded via a namespace (and not attached):
#   [1] fgsea_1.14.0           colorspace_2.0-0       rjson_0.2.20           ellipsis_0.3.1        
# [5] ggridges_0.5.3         circlize_0.4.13        qvalue_2.20.0          XVector_0.28.0        
# [9] GlobalOptions_0.1.2    rstudioapi_0.13        clue_0.3-60            farver_2.1.0          
# [13] urltools_1.7.3         graphlayouts_0.7.1     ggrepel_0.9.1          bit64_4.0.5           
# [17] AnnotationDbi_1.50.3   fansi_0.4.2            scatterpie_0.1.7       xml2_1.3.2            
# [21] splines_4.0.2          cachem_1.0.4           GOSemSim_2.14.2        geneplotter_1.66.0    
# [25] polyclip_1.10-0        jsonlite_1.7.2         annotate_1.66.0        cluster_2.1.2         
# [29] GO.db_3.11.4           png_0.1-7              ggforce_0.3.3          BiocManager_1.30.12   
# [33] compiler_4.0.2         httr_1.4.2             rvcheck_0.2.1          assertthat_0.2.1      
# [37] Matrix_1.3-2           fastmap_1.1.0          cli_3.6.2              tweenr_1.0.2          
# [41] prettyunits_1.1.1      tools_4.0.2            igraph_1.2.6           gtable_0.3.0          
# [45] glue_1.4.2             GenomeInfoDbData_1.2.3 reshape2_1.4.4         DO.db_2.9             
# [49] fastmatch_1.1-3        Rcpp_1.0.6             cellranger_1.1.0       vctrs_0.3.7           
# [53] ggraph_2.0.5           stringr_1.4.0          lifecycle_1.0.0        XML_3.99-0.6          
# [57] DOSE_3.14.0            europepmc_0.4.1        MASS_7.3-53.1          zlibbioc_1.34.0       
# [61] scales_1.1.1           tidygraph_1.2.0        hms_1.0.0              memoise_2.0.0         
# [65] gridExtra_2.3          downloader_0.4         ggfun_0.0.9            yulab.utils_0.0.4     
# [69] triebeard_0.3.0        stringi_1.5.3          RSQLite_2.2.7          genefilter_1.70.0     
# [73] BiocParallel_1.22.0    shape_1.4.6            rlang_0.4.10           pkgconfig_2.0.3       
# [77] bitops_1.0-7           lattice_0.20-41        purrr_0.3.4            labeling_0.4.2        
# [81] cowplot_1.1.1          bit_4.0.4              tidyselect_1.1.0       plyr_1.8.6            
# [85] magrittr_2.0.1         R6_2.5.0               generics_0.1.0         DBI_1.1.1             
# [89] pillar_1.6.0           withr_2.4.2            survival_3.2-11        RCurl_1.98-1.3        
# [93] tibble_3.1.1           crayon_1.4.1           utf8_1.2.1             viridis_0.6.0         
# [97] GetoptLong_1.0.5       progress_1.2.2         locfit_1.5-9.4         readxl_1.3.1          
# [101] data.table_1.14.0      blob_1.2.1             digest_0.6.27          xtable_1.8-4          
# [105] tidyr_1.1.3            gridGraphics_0.5-1     munsell_0.5.0          viridisLite_0.4.0     
# [109] ggplotify_0.1.0  



