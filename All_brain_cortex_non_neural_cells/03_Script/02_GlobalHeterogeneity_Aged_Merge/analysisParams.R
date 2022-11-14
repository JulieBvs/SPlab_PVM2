###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used byt the current
# analysis
#

ANALYSIS_STEP_NAME = "02_GlobalHeterogeneity_Aged_Merge"
PATH_ANALYSIS_OUTPUT = file.path( PATH_PROJECT_OUTPUT, ANALYSIS_STEP_NAME)

LITERAL_TITLE = "Integration of the 2 aged replicate datasets"


# Path
PATH_SEURAT_OBJECT_ADULT_REP1 = file.path(PATH_PROJECT_OUTPUT,
                                    "01_QC_Adult_rep1",
                                    paste0(outputFilesPrefix, "sc10x.rna.seurat.RDS"))
PATH_SEURAT_OBJECT_ADULT_REP2 = file.path(PATH_PROJECT_OUTPUT,
                                    "01_QC_Adult_rep2",
                                    paste0(outputFilesPrefix, "sc10x.rna.seurat.RDS"))
PATH_SEURAT_OBJECT_ADULT_REP3 = file.path(PATH_PROJECT_OUTPUT,
                                    "01_QC_Adult_rep3",
                                    paste0(outputFilesPrefix, "sc10x.rna.seurat.RDS"))


PATH_SEURAT_OBJECT_AGED_REP1 = file.path(PATH_PROJECT_OUTPUT,
                                          "01_QC_Aged_rep1",
                                          paste0(outputFilesPrefix, "sc10x.rna.seurat.RDS"))
PATH_SEURAT_OBJECT_AGED_REP2 = file.path(PATH_PROJECT_OUTPUT,
                                          "01_QC_Aged_rep2",
                                          paste0(outputFilesPrefix, "sc10x.rna.seurat.RDS"))


PATH_SEURAT_OBJECT_JUVENILE_REP1 = file.path(PATH_PROJECT_OUTPUT,
                                         "01_QC_Juvenile_rep1",
                                         paste0(outputFilesPrefix, "sc10x.rna.seurat.RDS"))
PATH_SEURAT_OBJECT_JUVENILE_REP2 = file.path(PATH_PROJECT_OUTPUT,
                                         "01_QC_Juvenile_rep2",
                                         paste0(outputFilesPrefix, "sc10x.rna.seurat.RDS"))
PATH_SEURAT_OBJECT_JUVENILE_REP3 = file.path(PATH_PROJECT_OUTPUT,
                                         "01_QC_Juvenile_rep3",
                                         paste0(outputFilesPrefix, "sc10x.rna.seurat.RDS"))
PATH_SEURAT_OBJECT_JUVENILE_REP4 = file.path(PATH_PROJECT_OUTPUT,
                                         "01_QC_Juvenile_rep4",
                                         paste0(outputFilesPrefix, "sc10x.rna.seurat.RDS"))


# Normalization parameters (see Seurat::NormalizeData())
DATA_NORM_METHOD      = "LogNormalize"
DATA_NORM_SCALEFACTOR = 10000

# Scaling parameters (see Seurat::ScaleData())
DATA_CENTER       = TRUE
DATA_SCALE        = FALSE
DATA_VARS_REGRESS = NULL  # c("nCount_RNA") for UMIs (NULL to ignore)


# Select highly variable genes
HVG_METHOD = "vst"
VARIABLE_FEATURES_MAXNB = 1500
VARIABLE_FEATURES_SHOWTOP = 10

# PCA parameters
PCA_NPC              = 50;  # Default number of dimensions to use for PCA (see Seurat::RunPCA())
PCA_PLOTS_NBDIMS     = 3;   # Number of dimensions to show in PCA-related plots
PCA_PLOTS_NBFEATURES = 15;  # Number of'top' features to show when plotting PCA loading


FINDCLUSTERS_RESOLUTION = 0.5
N_PCA = 10 # Number of PCA used for cell clustering
# Dimensionality reduction parameters (TSNE/UMAP)
DIMREDUC_USE_PCA_NBDIMS = 20;  # Number of dimensions to use from PCA results


# Parameters for identification of marker annotations for clusters (see Seurat::FindAllMarkers())
FINDMARKERS_METHOD    = "wilcox"  # Method used to identify markers
FINDMARKERS_ONLYPOS   = TRUE;     # Only consider overexpressed annotations for markers ? (if FALSE downregulated genes can also be markers)
FINDMARKERS_MINPCT    = 0.1;      # Only test genes that are detected in a minimum fraction of cells in either of the two populations. Speed up the function by not testing genes that are very infrequently expressed. Default is '0.1'.
FINDMARKERS_LOGFC_THR = 0.25;     # Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is '0.25'. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
FINDMARKERS_PVAL_THR  = 0.001;    # PValue threshold for identification of significative markers
FINDMARKERS_SHOWTOP   = 20;       # Number of marker genes to show in report and tables (NULL for all)





