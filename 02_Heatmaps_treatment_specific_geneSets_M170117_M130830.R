
# Heatmaps DEGs per cell line for M130830 and M1170117 ----
# Load the libraries
library(DESeq2)
library(dplyr)
library(genefilter)
library(ComplexHeatmap)
library(RColorBrewer)
library(biomaRt)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

# Define the directories ----
# pathScripts = "path/to/scripts/"
pathOutput = "path/to/output/"
# Set the output directory as working directory (where outputs are to be saved)
wd = pathOutput
setwd(wd)

# Load the pre-computed outputs
load("path/to/folder/where/file/is/stored/outputs.RData") 

# M1170117
# Retrieve the set of DEGs between each condition vs DMSO
genes.DMSO.All.M170117 <- c(rownames(res_TGFb_vs_DMSO.M170117_05padj),
                            rownames(res_MEKi_vs_DMSO.M170117_05padj),
                            rownames(res_TGFb.MEKi_vs_DMSO.M170117_05padj))
length(genes.DMSO.All.M170117)
genes.DMSO.All.M170117 <- unique(genes.DMSO.All.M170117)
length(genes.DMSO.All.M170117)

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
df_norm_M170 <- df_norm_M170[genes.DMSO.All.M170117, ]
dim(df_norm_M170)
# Set the order of the treatments/samples in the expression matrix
df_norm_M170 <- df_norm_M170[ , coldata.M170117$sample_full]
# Scale the data
df_zsc_M170 <- t(as.data.frame((scale(t(df_norm_M170)))))
# Format the column names
colnames(df_zsc_M170)
colnames(df_zsc_M170) <- sub(".M170117", "", colnames(df_zsc_M170))
colnames(df_zsc_M170)

# Define the annotations
treatment.ann <- sapply(strsplit(colnames(df_zsc_M170),"\\_"), function(x) x[1])
treatment.ann <- factor(treatment.ann, levels = c("DMSO", "TGFb", "MEKi", "TGFb + MEKi"))
# Define the color for the annotations
my.col <- brewer.pal(length(unique(treatment.ann)), "Set3")

column_ha <- HeatmapAnnotation(treatment = treatment.ann, 
                               col = list(treatment = c("DMSO" = my.col[1], 
                                                        "TGFb" = my.col[2], 
                                                        "MEKi" = my.col[3], 
                                                        "TGFb + MEKi" = my.col[4])))
# Make the heatmap
set.seed(1234)
heatm.M170 <- Heatmap(df_zsc_M170, 
                  show_row_names = F, 
                  show_column_names = F,
                  column_title = "M170117", 
                  name = "zscore", 
                  row_split = 7, 
                  column_order = colnames(df_zsc_M170),
                  top_annotation = column_ha) 


# somePDFPath = "Heatmap_M170117.pdf"
# pdf(file=somePDFPath, height = 6, width = 4) 
heatm.M170
# dev.off()

# M130830
# Retrieve the set of DEGs between each condition vs DMSO
genes.DMSO.All.M130830 <- c(rownames(res_TGFb_vs_DMSO.M130830_05padj), 
                            rownames(res_MEKi_vs_DMSO.M130830_05padj), 
                            rownames(res_TGFb.MEKi_vs_DMSO.M130830_05padj))
length(genes.DMSO.All.M130830)
genes.DMSO.All.M130830 <- unique(genes.DMSO.All.M130830)
length(genes.DMSO.All.M130830)

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
df_norm_M130 <- df_norm_M130[genes.DMSO.All.M130830, ]
dim(df_norm_M130)
# Set the order of the treatments/samples in the expression matrix
df_norm_M130 <- df_norm_M130[ , coldata.M130830$sample_full]
# Scale the data
df_zsc_M130 <- t(as.data.frame((scale(t(df_norm_M130)))))
# Format the column names
colnames(df_zsc_M130)
colnames(df_zsc_M130) <- sub(".M130830", "", colnames(df_zsc_M130))
colnames(df_zsc_M130)

# Define the annotations
treatment.ann <- sapply(strsplit(colnames(df_zsc_M130),"\\_"), function(x) x[1])
treatment.ann <- factor(treatment.ann, levels = c("DMSO", "TGFb", "MEKi", "TGFb + MEKi"))
# Define the color for the annotations
my.col <- brewer.pal(length(unique(treatment.ann)), "Set3")

column_ha <- HeatmapAnnotation(treatment = treatment.ann, 
                               col = list(treatment = c("DMSO" = my.col[1], 
                                                        "TGFb" = my.col[2],
                                                        "MEKi" = my.col[3], 
                                                        "TGFb + MEKi" = my.col[4])))

# Make the heatmap
set.seed(1234)
heatm.M130 <- Heatmap(df_zsc_M130, 
                  show_row_names = F, 
                  show_column_names = F,
                  column_title = "M130830", 
                  name = "zscore", 
                  row_split = 5, 
                  column_order = colnames(df_zsc_M130),
                  top_annotation = column_ha) 


# somePDFPath = "Heatmap_M130830.pdf"
# pdf(file=somePDFPath, height = 6, width = 4) 
heatm.M130
# dev.off()

# Retrieve the genes that belong to heatmap slices representing DEGs Up in specific treatments ----
# Extract the genes for the specific slices per cell line

# M170117
# The zscore scale data was computed above
df_zsc_M170
# Retrieve the order of the genes in the hetamap as defined by the clustering
k.M170117 <- row_order(heatm.M170)
k.gene.names.M170117 <- lapply(k.M170117, function(x) rownames(df_zsc_M170[x, ]))
length(k.M170117)
# M170117 slice 6
k.gene.names.M170117.slice6 <- unlist(k.gene.names.M170117[6])
# Retrieve the gene IDs from the slice
ext.gene.names[k.gene.names.M170117.slice6, ]
write.table(ext.gene.names[k.gene.names.M170117.slice6, ], "Genes_M170slice6.txt", quote = FALSE, row.names = FALSE)
# M170117 slice 7
k.gene.names.M170117.slice7 <- unlist(k.gene.names.M170117[7])
# Retrieve the gene IDs from the slice
ext.gene.names[k.gene.names.M170117.slice7, ]
write.table(ext.gene.names[k.gene.names.M170117.slice7, ], "Genes_M170slice7.txt", quote = FALSE, row.names = FALSE)


# M130830
# The zscore scale data was computed above
df_zsc_M130
# Retrieve the order of the genes in the hetamap as defined by the clustering
k.M130830 <- row_order(heatm.M130)
k.gene.names.M130830 <- lapply(k.M130830, function(x) rownames(df_zsc_M130[x, ]))
length(k.M130830)
# M130830 slice 1
k.gene.names.M130830.slice1 <- unlist(k.gene.names.M130830[1])
# Retrieve the gene IDs from the slice
ext.gene.names[k.gene.names.M130830.slice1, ]
write.table(ext.gene.names[k.gene.names.M130830.slice1, ], "Genes_M130slice1.txt", quote = FALSE, row.names = FALSE)
# M130830 slice 4
k.gene.names.M130830.slice4 <- unlist(k.gene.names.M130830[4])
# Retrieve the gene IDs from the slice
ext.gene.names[k.gene.names.M130830.slice4, ]
write.table(ext.gene.names[k.gene.names.M130830.slice4, ], "Genes_M130slice4.txt", quote = FALSE, row.names = FALSE)

# Converting apoptosis GO term genes from mouse to Human homologs ----
# Full list of mouse apoptosis related genes (available in Github under the "files" folder)
mouse.GO.apopt <- readxl::read_excel("path/to/folder/where/file/is/stored/GO_term_summary_20230123_062954.xlsx", col_names = T) 
mouse.GO.apopt.genes <- unique(mouse.GO.apopt$Symbol)
length(mouse.GO.apopt.genes)
# Translate the genes from the mouse to the human
# Define the databases for the query
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "dec2021.archive.ensembl.org")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "dec2021.archive.ensembl.org")

# Define the genes to perform the query for
my.genes <- mouse.GO.apopt.genes
length(my.genes)

# Make the query of the required information
mouse.to.human.apopt.genes <- getLDS(attributes = c('ensembl_gene_id',
                                                    'external_gene_name'), 
                                     filters = 'external_gene_name', 
                                     values = my.genes,
                                     mart = mouse,
                                     attributesL = c('ensembl_gene_id',
                                                     'external_gene_name'), 
                                     martL = human, 
                                     uniqueRows=T)

dim(mouse.to.human.apopt.genes)
# Remove duplicated genes
mouse.to.human.apopt.genes.FTD <- mouse.to.human.apopt.genes[!duplicated(mouse.to.human.apopt.genes$Gene.name), ]
dim(mouse.to.human.apopt.genes.FTD)
mouse.to.human.apopt.genes.FTD <- mouse.to.human.apopt.genes.FTD[!duplicated(mouse.to.human.apopt.genes.FTD$Gene.name.1), ]
dim(mouse.to.human.apopt.genes.FTD)
# Check if there are still duplicated genes
length(unique(mouse.to.human.apopt.genes.FTD$Gene.name))
length(unique(mouse.to.human.apopt.genes.FTD$Gene.name.1))
# Export the file
xlsx::write.xlsx(mouse.to.human.apopt.genes.FTD, "Apoptosis_related_genes_mouse_to_human.xlsx")

# Import the list of mouse negative regulators of apoptosis (available in Github under the "files" folder)
mouse.GO.apopt.neg <- readxl::read_excel("path/to/folder/where/file/is/stored/GO_term_0043066_neg_apoptosis.xlsx", col_names = T) 
dim(mouse.GO.apopt.neg)
mouse.GO.apopt.neg.genes <- unique(mouse.GO.apopt.neg$Symbol)
length(mouse.GO.apopt.neg.genes)
# Retrieve the human gene names using the full set of apoptosis genes from above 
mouse.GO.apopt.neg.genes.to.human <- mouse.to.human.apopt.genes.FTD[mouse.to.human.apopt.genes.FTD$Gene.name %in% mouse.GO.apopt.neg.genes, ]
dim(mouse.GO.apopt.neg.genes.to.human)
# Check if there are any duplicated genes
length(unique(mouse.GO.apopt.neg.genes.to.human$Gene.name))
length(unique(mouse.GO.apopt.neg.genes.to.human$Gene.name.1))
# Export the file
xlsx::write.xlsx(mouse.GO.apopt.neg.genes.to.human, "Apoptosis_related_genes_mouse_to_human_negative_regulators.xlsx")

