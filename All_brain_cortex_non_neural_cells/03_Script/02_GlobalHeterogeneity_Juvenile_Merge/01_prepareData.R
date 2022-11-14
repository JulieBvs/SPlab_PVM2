# ===========================================
# This script reads the sc10x data from
# QC and filtering analysis
# ===========================================



#######################
# Read data
#######################

## @knitr loadData

# # Load seurat object
# sc10x.rna.seurat.ad.rep1 = readRDS(PATH_SEURAT_OBJECT_ADULT_REP1)
# cat("<br>Number of cells in loaded cells from replicate 1:", length( Cells( sc10x.rna.seurat.ad.rep1)))
# sc10x.rna.seurat.ad.rep1[["sample"]] = "adult_rep1"

sc10x.rna.seurat.ju.rep1 = readRDS(PATH_SEURAT_OBJECT_JUVENILE_REP1)
cat("<br>Number of cells in loaded cells from replicate 1:", length( Cells( sc10x.rna.seurat.ju.rep1)))
sc10x.rna.seurat.ju.rep1[["sample"]] = "juvenile_rep1"

sc10x.rna.seurat.ju.rep2 = readRDS(PATH_SEURAT_OBJECT_JUVENILE_REP2)
cat("<br>Number of cells in loaded cells from replicate 2:", length( Cells( sc10x.rna.seurat.ju.rep2)))
sc10x.rna.seurat.ju.rep2[["sample"]] = "juvenile_rep2"

sc10x.rna.seurat.ju.rep3 = readRDS(PATH_SEURAT_OBJECT_JUVENILE_REP3)
cat("<br>Number of cells in loaded cells from replicate 3:", length( Cells( sc10x.rna.seurat.ju.rep3)))
sc10x.rna.seurat.ju.rep3[["sample"]] = "juvenile_rep3"

sc10x.rna.seurat.ju.rep4 = readRDS(PATH_SEURAT_OBJECT_JUVENILE_REP4)
cat("<br>Number of cells in loaded cells from replicate 4:", length( Cells( sc10x.rna.seurat.ju.rep4)))
sc10x.rna.seurat.ju.rep4[["sample"]] = "juvenile_rep4"


# Merge seurat object from replicate of a same time point
sc10x.rna.seurat.ju.all <- merge(sc10x.rna.seurat.ju.rep1, 
                                 y = c(sc10x.rna.seurat.ju.rep2, sc10x.rna.seurat.ju.rep3, sc10x.rna.seurat.ju.rep4), 
                                 add.cell.ids = c("Juvenile_rep1", "Juvenile_rep2", "Juvenile_rep3", "Juvenile_rep4"), 
                                 project = "Weirdophage")



# NORMALIZE DATA
# --------------

## @knitr normalizeData

# Normalize RNA data with log normalization
sc10x.rna.seurat.ju.all = NormalizeData( object = sc10x.rna.seurat.ju.all,
                                         normalization.method = DATA_NORM_METHOD,
                                         scale.factor = DATA_NORM_SCALEFACTOR,
                                         verbose = .VERBOSE)

sc10x.rna.seurat.ju.all = ScaleData( object    = sc10x.rna.seurat.ju.all,
                                     do.center = DATA_CENTER,
                                     do.scale  = DATA_SCALE,
                                     vars.to.regress = DATA_VARS_REGRESS,
                                     verbose = .VERBOSE)




