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
