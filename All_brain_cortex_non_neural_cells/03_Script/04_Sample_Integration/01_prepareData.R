# ===========================================
# This script reads the sc10x data from
# QC and filtering analysis
# ===========================================



#######################
# Read data
#######################

## @knitr loadData

# Load seurat object
sc10x.rna.seurat.ju = readRDS(PATH_SEURAT_OBJECT_JUVENILE)
cat("<br>Number of cells in dataset from juvenile sample: ", length( Cells( sc10x.rna.seurat.ju)))

sc10x.rna.seurat.ad = readRDS(PATH_SEURAT_OBJECT_ADULT)
cat("<br>Number of cells in dataset from adult sample:", length( Cells( sc10x.rna.seurat.ad)))

sc10x.rna.seurat.ag = readRDS(PATH_SEURAT_OBJECT_AGED)
cat("<br>Number of cells in dataset from aged sample:", length( Cells( sc10x.rna.seurat.ag)))
