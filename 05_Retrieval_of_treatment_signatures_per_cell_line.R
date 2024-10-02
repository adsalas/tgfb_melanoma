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
