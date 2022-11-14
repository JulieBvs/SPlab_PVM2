###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used byt the current
# analysis
#

ANALYSIS_STEP_NAME = "04_Sample_Integration"
PATH_ANALYSIS_OUTPUT = file.path( PATH_PROJECT_OUTPUT, ANALYSIS_STEP_NAME)

LITERAL_TITLE = "Integration of the 3 sample (merged datasets)"


# Path
PATH_SEURAT_OBJECT_JUVENILE = file.path(PATH_PROJECT_OUTPUT,
                                    "02_GlobalHeterogeneity_Juvenile_Merge",
                                    paste0(outputFilesPrefix, "QCfiltered_sc10x.rna.seurat.ju.all.RDS"))

PATH_SEURAT_OBJECT_ADULT = file.path(PATH_PROJECT_OUTPUT,
                                    "02_GlobalHeterogeneity_Adult_Merge",
                                    paste0(outputFilesPrefix, "QCfiltered_sc10x.rna.seurat.ad.all.RDS"))

PATH_SEURAT_OBJECT_AGED = file.path(PATH_PROJECT_OUTPUT,
                                    "02_GlobalHeterogeneity_Aged_Merge",
                                    paste0(outputFilesPrefix, "QCfiltered_sc10x.rna.seurat.ag.all.RDS"))

# Select highly variable genes
HVG_METHOD = "vst"
VARIABLE_FEATURES_MAXNB = 2000
VARIABLE_FEATURES_SHOWTOP = 10

#
N_PCA = 30
CL_RESOLUTION = 0.4

#Marker genes
FINDMARKERS_PVAL_THR = 0.05
FINDCONSERVEDMARKERS_SHOWTOP = 10
FINDMARKERS_SHOWTOP = 50




