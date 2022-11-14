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

sc10x.rna.seurat.ag.rep1 = readRDS(PATH_SEURAT_OBJECT_AGED_REP1)
cat("<br>Number of cells in loaded cells from replicate 1:", length( Cells( sc10x.rna.seurat.ag.rep1)))
sc10x.rna.seurat.ag.rep1[["sample"]] = "aged_rep1"

sc10x.rna.seurat.ag.rep2 = readRDS(PATH_SEURAT_OBJECT_AGED_REP2)
cat("<br>Number of cells in loaded cells from replicate 2:", length( Cells( sc10x.rna.seurat.ag.rep2)))
sc10x.rna.seurat.ag.rep2[["sample"]] = "aged_rep2"


# Merge seurat object from replicate of a same time point
sc10x.rna.seurat.ag.all <- merge(sc10x.rna.seurat.ag.rep1, 
                                 y = c(sc10x.rna.seurat.ag.rep2), 
                                 add.cell.ids = c("Aged_rep1", "Aged_rep2"), 
                                 project = "Weirdophage")



# NORMALIZE DATA
# --------------

## @knitr normalizeData

# Normalize RNA data with log normalization
sc10x.rna.seurat.ag.all = NormalizeData( object = sc10x.rna.seurat.ag.all,
                                         normalization.method = DATA_NORM_METHOD,
                                         scale.factor = DATA_NORM_SCALEFACTOR,
                                         verbose = .VERBOSE)

sc10x.rna.seurat.ag.all = ScaleData( object    = sc10x.rna.seurat.ag.all,
                                     do.center = DATA_CENTER,
                                     do.scale  = DATA_SCALE,
                                     vars.to.regress = DATA_VARS_REGRESS,
                                     verbose = .VERBOSE)




