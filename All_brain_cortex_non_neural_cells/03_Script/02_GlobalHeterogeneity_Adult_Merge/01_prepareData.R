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

sc10x.rna.seurat.ad.rep2 = readRDS(PATH_SEURAT_OBJECT_ADULT_REP2)
cat("<br>Number of cells in loaded cells from replicate 2:", length( Cells( sc10x.rna.seurat.ad.rep2)))
sc10x.rna.seurat.ad.rep2[["sample"]] = "adult_rep2"

sc10x.rna.seurat.ad.rep3 = readRDS(PATH_SEURAT_OBJECT_ADULT_REP3)
cat("<br>Number of cells in loaded cells from replicate 3:", length( Cells( sc10x.rna.seurat.ad.rep3)))
sc10x.rna.seurat.ad.rep3[["sample"]] = "adult_rep3"


# Merge seurat object from replicate of a same time point
sc10x.rna.seurat.ad.all <- merge(sc10x.rna.seurat.ad.rep2, 
                                 y = c(sc10x.rna.seurat.ad.rep3), 
                                 add.cell.ids = c("Adult_rep2", "Adult_rep3"), 
                                 project = "Weirdophage")



# NORMALIZE DATA
# --------------

## @knitr normalizeData

# Normalize RNA data with log normalization
sc10x.rna.seurat.ad.all = NormalizeData( object = sc10x.rna.seurat.ad.all,
                                         normalization.method = DATA_NORM_METHOD,
                                         scale.factor = DATA_NORM_SCALEFACTOR,
                                         verbose = .VERBOSE)

sc10x.rna.seurat.ad.all = ScaleData( object    = sc10x.rna.seurat.ad.all,
                                     do.center = DATA_CENTER,
                                     do.scale  = DATA_SCALE,
                                     vars.to.regress = DATA_VARS_REGRESS,
                                     verbose = .VERBOSE)




