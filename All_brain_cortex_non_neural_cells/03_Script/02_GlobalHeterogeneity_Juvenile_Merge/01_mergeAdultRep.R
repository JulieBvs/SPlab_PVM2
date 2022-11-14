






#######################
# Read data
#######################

## @knitr loadData

# Load seurat object
sc10x.rna.seurat.ad.rep1 = readRDS(PATH_SEURAT_OBJECT_ADULT_REP1)
cat("<br>Number of cells in loaded cells from replicate 1:", length( Cells( sc10x.rna.seurat.ad.rep1)))
sc10x.rna.seurat.ad.rep1[["sample"]] = "adult_rep1"

sc10x.rna.seurat.ad.rep2 = readRDS(PATH_SEURAT_OBJECT_ADULT_REP2)
cat("<br>Number of cells in loaded cells from replicate 2:", length( Cells( sc10x.rna.seurat.ad.rep2)))
sc10x.rna.seurat.ad.rep2[["sample"]] = "adult_rep2"

sc10x.rna.seurat.ad.rep3 = readRDS(PATH_SEURAT_OBJECT_ADULT_REP3)
cat("<br>Number of cells in loaded cells from replicate 3:", length( Cells( sc10x.rna.seurat.ad.rep3)))
sc10x.rna.seurat.ad.rep3[["sample"]] = "adult_rep3"


# Merge seurat object from replicate of a same time point
sc10x.rna.seurat.ad.all <- merge(sc10x.rna.seurat.ad.rep1, 
                                 y = c(sc10x.rna.seurat.ad.rep2, sc10x.rna.seurat.ad.rep3), 
                                 add.cell.ids = c("Adult_rep1", "Adult_rep2", "Adult_rep3"), 
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






#######################
# Identification of the highly variable genes
#######################

## @knitr findVariableGenes_seuratMethod

# Find most variable features
sc10x.rna.seurat.ad.all = FindVariableFeatures( object = sc10x.rna.seurat.ad.all, 
                                                selection.method = "vst",
                                                loess.span = 0.3,
                                                clip.max = "auto",
                                                mean.function = ExpMean, 
                                                dispersion.function = LogVMR,
                                                nfeatures = VARIABLE_FEATURES_MAXNB,
                                                verbose = .VERBOSE);

variablesGenesStats = paste(length( VariableFeatures( sc10x.rna.seurat.ad.all)), "/", nrow( sc10x.rna.seurat.ad.all));


## @knitr findVariableGenes_summaryPlot

# Prepare a Variance/Expression plot highlighting variable genes and add names of most variable genes
suppressMessages(suppressWarnings(LabelPoints( plot = VariableFeaturePlot( sc10x.rna.seurat.ad.all) + theme(legend.position = "none"), 
                                               points = head( VariableFeatures( sc10x.rna.seurat.ad.all), 10), 
                                               repel = TRUE)));


## @knitr findVariableGenes_summaryTable

# Extract variable genes info as data.frame
variableAnnotationsDT = head( HVFInfo( object = sc10x.rna.seurat.ad.all, assay = "RNA", selection.method = 'vst')[VariableFeatures( sc10x.rna.seurat.ad.all),], VARIABLE_FEATURES_SHOWTOP);
variableAnnotationsDT = cbind("Gene" = rownames(variableAnnotationsDT), variableAnnotationsDT);

# Create a table in report containing information about top variable genes
datatable( variableAnnotationsDT, # Set annotation names as column instead of rownames so datatable handles column search properly
           class = "compact",
           filter="top",
           rownames = FALSE,
           colnames = c("Gene", "Avg. Expression", "Variance", "Std. Variance"),
           extensions = c('Buttons', 'Select'),
           options = list(dom = "<'row'<'col-sm-8'B><'col-sm-4'f>> <'row'<'col-sm-12'l>> <'row'<'col-sm-12'rt>> <'row'<'col-sm-12'ip>>", # Set elements for CSS formatting ('<Blf><rt><ip>')
                          autoWidth = FALSE,
                          buttons = exportButtonsListDT,
                          columnDefs = list( 
                            list( # Center all columns except first one
                              targets = 1:(ncol( variableAnnotationsDT)-1),
                              className = 'dt-center'),
                            list( # Set renderer function for 'float' type columns
                              targets = 1:(ncol( variableAnnotationsDT)-1),
                              render = htmlwidgets::JS( "function ( data, type, row ) {return type === 'export' ? data : data.toFixed(4);}"))), 
                          #fixedColumns = TRUE, # Does not work well with filter on this column
                          #fixedHeader = TRUE, # Does not work well with 'scrollX'
                          lengthMenu = list(c( 10, 50, 100, -1),
                                            c( 10, 50, 100, "All")),
                          orderClasses = FALSE, # Disable flag for CSS to highlight columns used for ordering (for performance)
                          processing = TRUE, 
                          #search.regex= TRUE, # Does not work well with 'search.smart'
                          search.smart = TRUE,
                          select = TRUE, # Enable ability to select rows
                          scrollCollapse = TRUE,
                          scroller = TRUE,  # Only load visible data
                          scrollX = TRUE,
                          scrollY = "525px",
                          stateSave = TRUE));








#####################
# Global heterogeneity
#####################

sc10x.rna.seurat.ad.all <- RunPCA( object   = sc10x.rna.seurat.ad.all, 
                                   features = VariableFeatures( sc10x.rna.seurat.ad.all),
                                   npcs     = PCA_NPC,
                                   verbose  = .VERBOSE)

# # Identify clusters of cells by graph approach
# sc10x.rna.seurat.ad.all <- FindNeighbors(object  = sc10x.rna.seurat.ad.all, 
#                                          dims    = 1:10,
#                                          verbose = .VERBOSE);
# 
# sc10x.rna.seurat.ad.all <- FindClusters(object              = sc10x.rna.seurat.ad.all, 
#                                         # reduction.type     = "pca", 
#                                         resolution         = FINDCLUSTERS_RESOLUTION, 
#                                         random.seed        = 1511,
#                                         temp.file.location = "/tmp/",
#                                         verbose            = .VERBOSE);


sc10x.rna.seurat.ad.all = RunTSNE( sc10x.rna.seurat.ad.all, dims = 1:DIMREDUC_USE_PCA_NBDIMS, seed.use = 1511);
sc10x.rna.seurat.ad.all = RunUMAP( sc10x.rna.seurat.ad.all, dims = 1:DIMREDUC_USE_PCA_NBDIMS, seed.use = 1511);

DimPlot(sc10x.rna.seurat.ad.all, reduction = 'umap', group.by = 'sample')
DimPlot(sc10x.rna.seurat.ad.all, reduction = 'umap', group.by = 'sample', split.by = 'sample')


