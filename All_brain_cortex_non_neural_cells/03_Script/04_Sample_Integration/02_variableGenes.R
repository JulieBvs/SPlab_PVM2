# ===========================================
# This script finds the highly variable genes
# among the cell population for each replicates
# ===========================================



#######################
# Juvenile - seurat method
#######################

## @knitr findVariableGenes_seuratMethodData_rep1

sc10x.rna.seurat.ju = FindVariableFeatures( object = sc10x.rna.seurat.ju, 
                                              selection.method = HVG_METHOD,
                                              loess.span = 0.3,
                                              clip.max = "auto",
                                              mean.function = ExpMean, 
                                              dispersion.function = LogVMR,
                                              nfeatures = VARIABLE_FEATURES_MAXNB,
                                              verbose = .VERBOSE)

variablesGenesStats = paste(length( VariableFeatures( sc10x.rna.seurat.ju)), "/", nrow( sc10x.rna.seurat.ju))



## @knitr findVariableGenes_summaryPlot_rep1

# Prepare a Variance/Expression plot highlighting variable genes and add names of most variable genes
suppressMessages(suppressWarnings(LabelPoints( plot = VariableFeaturePlot( sc10x.rna.seurat.ju) + theme(legend.position = "none"), 
                                               points = head( VariableFeatures( sc10x.rna.seurat.ju), 10), 
                                               repel = TRUE)));


## @knitr findVariableGenes_summaryTable_rep1

# Extract variable genes info as data.frame
variableAnnotationsDT = head( HVFInfo( object = sc10x.rna.seurat.ju, assay = "RNA", selection.method = HVG_METHOD)[VariableFeatures( sc10x.rna.seurat.ju),], VARIABLE_FEATURES_SHOWTOP);
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





#######################
# Adult - seurat method
#######################

## @knitr findVariableGenes_seuratMethodData_rep2

sc10x.rna.seurat.ad = FindVariableFeatures( object = sc10x.rna.seurat.ad, 
                                              selection.method = HVG_METHOD,
                                              loess.span = 0.3,
                                              clip.max = "auto",
                                              mean.function = ExpMean, 
                                              dispersion.function = LogVMR,
                                              nfeatures = VARIABLE_FEATURES_MAXNB,
                                              verbose = .VERBOSE)

variablesGenesStats = paste(length( VariableFeatures( sc10x.rna.seurat.ad)), "/", nrow( sc10x.rna.seurat.ad))



## @knitr findVariableGenes_summaryPlot_rep2

# Prepare a Variance/Expression plot highlighting variable genes and add names of most variable genes
suppressMessages(suppressWarnings(LabelPoints( plot = VariableFeaturePlot( sc10x.rna.seurat.ad) + theme(legend.position = "none"), 
                                               points = head( VariableFeatures( sc10x.rna.seurat.ad), 10), 
                                               repel = TRUE)));


## @knitr findVariableGenes_summaryTable_rep2

# Extract variable genes info as data.frame
variableAnnotationsDT = head( HVFInfo( object = sc10x.rna.seurat.ad, assay = "RNA", selection.method = HVG_METHOD)[VariableFeatures( sc10x.rna.seurat.ad),], VARIABLE_FEATURES_SHOWTOP);
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





#######################
# Aged - seurat method
#######################

## @knitr findVariableGenes_seuratMethodData_rep3

sc10x.rna.seurat.ag = FindVariableFeatures( object = sc10x.rna.seurat.ag, 
                                              selection.method = HVG_METHOD,
                                              loess.span = 0.3,
                                              clip.max = "auto",
                                              mean.function = ExpMean, 
                                              dispersion.function = LogVMR,
                                              nfeatures = VARIABLE_FEATURES_MAXNB,
                                              verbose = .VERBOSE)

variablesGenesStats = paste(length( VariableFeatures( sc10x.rna.seurat.ag)), "/", nrow( sc10x.rna.seurat.ag))



## @knitr findVariableGenes_summaryPlot_rep3

# Prepare a Variance/Expression plot highlighting variable genes and add names of most variable genes
suppressMessages(suppressWarnings(LabelPoints( plot = VariableFeaturePlot( sc10x.rna.seurat.ag) + theme(legend.position = "none"), 
                                               points = head( VariableFeatures( sc10x.rna.seurat.ag), 10), 
                                               repel = TRUE)));


## @knitr findVariableGenes_summaryTable_rep3

# Extract variable genes info as data.frame
variableAnnotationsDT = head( HVFInfo( object = sc10x.rna.seurat.ag, assay = "RNA", selection.method = HVG_METHOD)[VariableFeatures( sc10x.rna.seurat.ag),], VARIABLE_FEATURES_SHOWTOP);
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






