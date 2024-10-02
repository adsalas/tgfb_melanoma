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
