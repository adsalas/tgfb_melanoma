
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

