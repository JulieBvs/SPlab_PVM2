# ===========================================
# This script combine the seurat object from the 
# replicates in one combined seurat object.
# ===========================================




## @knitr pca_tsne_umap

# Run PCA
sc10x.rna.seurat.ad.all <- RunPCA( object   = sc10x.rna.seurat.ad.all, 
                                   features = VariableFeatures( sc10x.rna.seurat.ad.all),
                                   npcs     = PCA_NPC,
                                   verbose  = .VERBOSE)

# Run T-sne and UMAP
sc10x.rna.seurat.ad.all = RunTSNE( sc10x.rna.seurat.ad.all, dims = 1:DIMREDUC_USE_PCA_NBDIMS, seed.use = 1511);
sc10x.rna.seurat.ad.all = RunUMAP( sc10x.rna.seurat.ad.all, dims = 1:DIMREDUC_USE_PCA_NBDIMS, seed.use = 1511);



## @knitr visualizationMergedData

# Plot the projection with relicate colors
print( suppressWarnings( DimPlot( sc10x.rna.seurat.ad.all, reduction = ifelse(exists("useReduction"), useReduction, "umap"), group.by = "sample")) + 
         theme( axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                plot.margin = margin(0, 0, 0, 0, "cm")))

print( suppressWarnings( DimPlot( sc10x.rna.seurat.ad.all, reduction = ifelse(exists("useReduction"), useReduction, "umap"), split.by = "sample", group.by = "sample")) + 
         theme( axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                legend.position = "none", 
                plot.margin = margin(0, 0, 0, 0, "cm")) +
         coord_fixed(1))
cat(" \n \n") # Required for '.tabset'



#######################
# Cell clustering
#######################

## @knitr cellClustering
sc10x.rna.seurat.ad.all <- FindNeighbors(sc10x.rna.seurat.ad.all, reduction = "pca", dims = 1:N_PCA)
sc10x.rna.seurat.ad.all <- FindClusters(sc10x.rna.seurat.ad.all, resolution = FINDCLUSTERS_RESOLUTION)


# Plot the projection with colored clusters
## @knitr cellClustering_visualization
print( suppressWarnings( DimPlot( sc10x.rna.seurat.ad.all, reduction = ifelse(exists("useReduction"), useReduction, "umap"), 
                                  label = FALSE, cols = col_rainbow[1:nlevels( Idents( sc10x.rna.seurat.ad.all))])) + 
         theme( axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                plot.margin = margin(0, 0, 0, 0, "cm")))

print( suppressWarnings( DimPlot( sc10x.rna.seurat.ad.all, reduction = ifelse(exists("useReduction"), useReduction, "umap"), split.by = "sample", cols = col_rainbow[1:nlevels( Idents( sc10x.rna.seurat.ad.all))])) + 
         theme( axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                legend.position = "none", 
                plot.margin = margin(0, 0, 0, 0, "cm")) +
         coord_fixed(1))
cat(" \n \n") # Required for '.tabset'






# ---------------------------------------------------------------
# ---------------------------------------------------------------
# MARKER GENES - IDENTIFICATION
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr heterogeneity_markerGenes

Idents( sc10x.rna.seurat.ad.all) <- "seurat_clusters"

# Identify marker genes
markers = FindAllMarkers( object          = sc10x.rna.seurat.ad.all,
                          test.use        = FINDMARKERS_METHOD,
                          only.pos        = FINDMARKERS_ONLYPOS,
                          min.pct         = FINDMARKERS_MINPCT,
                          logfc.threshold = FINDMARKERS_LOGFC_THR,
                          verbose         = .VERBOSE);

# #Print the datatable with all markers
# print( datatable( markers[ which( markers$p_val_adj < FINDMARKERS_PVAL_THR), ]), caption = paste( "All markers by clusters with adj.pval <", FINDMARKERS_PVAL_THR))

# Filter markers by cluster (TODO: check if downstream code works when no markers found)
topMarkers = by( markers, markers[["cluster"]], function(x)
{
  # Filter markers based on adjusted PValue
  x = x[ x[["p_val_adj"]] < FINDMARKERS_PVAL_THR, , drop = FALSE];
  # Sort by decreasing logFC
  x = x[ order(abs(x[["avg_log2FC"]]), decreasing = TRUE), , drop = FALSE ]
  # Return top ones
  return( if(is.null( FINDMARKERS_SHOWTOP)) x else head( x, n = FINDMARKERS_SHOWTOP));
});

# Filter markers by cluster (TODO: check if downstream code works when no markers found)
topMarkersHalf = by( markers, markers[["cluster"]], function(x)
{
  # Filter markers based on adjusted PValue
  x = x[ x[["p_val_adj"]] < FINDMARKERS_PVAL_THR, , drop = FALSE];
  # Sort by decreasing logFC
  x = x[ order(abs(x[["avg_log2FC"]]), decreasing = TRUE), , drop = FALSE ]
  # Return top ones
  return( if(is.null( FINDMARKERS_SHOWTOP)) x else head( x, n = floor( FINDMARKERS_SHOWTOP /2)));
});

# Merge marker genes in a single data.frame and render it as datatable
topMarkersDF = do.call( rbind, topMarkers);
topMarkersHalfDF = do.call( rbind, topMarkersHalf);
# Select and order columns to be shown in datatable
topMarkersDT = topMarkersDF[c("gene", "cluster", "avg_log2FC", "p_val_adj")]
topMarkersHalfDT = topMarkersHalfDF[c("gene", "cluster", "avg_log2FC", "p_val_adj")]



# ---------------------------------------------------------------
# ---------------------------------------------------------------
# MARKER GENES - TABLE
# ---------------------------------------------------------------
# ---------------------------------------------------------------
## @knitr heterogeneity_markerGenes_table

# Create datatable
htmltools::tagList(datatable( topMarkersDT, 
                              class = "compact",
                              filter="top",
                              rownames = FALSE,
                              colnames = c("Gene", "Cluster", "Avg. LogFC", "Adj. Pvalue"),
                              caption = paste(ifelse( is.null( FINDMARKERS_SHOWTOP), "All", paste("Top", FINDMARKERS_SHOWTOP)), "marker genes for each cluster"),
                              extensions = c('Buttons', 'Select'),
                              options = list(dom = "<'row'<'col-sm-8'B><'col-sm-4'f>> <'row'<'col-sm-12'l>> <'row'<'col-sm-12'rt>> <'row'<'col-sm-12'ip>>", # Set elements for CSS formatting ('<Blf><rt><ip>')
                                             autoWidth = FALSE,
                                             buttons = exportButtonsListDT,
                                             columnDefs = list( 
                                               list( # Center all columns except first one
                                                 targets = 1:(ncol( topMarkersDT)-1),
                                                 className = 'dt-center'),
                                               list( # Set renderer function for 'float' type columns (LogFC)
                                                 targets = ncol( topMarkersDT)-2,
                                                 render = htmlwidgets::JS("function ( data, type, row ) {return type === 'export' ? data : data.toFixed(4);}")),
                                               list( # Set renderer function for 'scientific' type columns (PValue)
                                                 targets = ncol( topMarkersDT)-1,
                                                 render = htmlwidgets::JS( "function ( data, type, row ) {return type === 'export' ? data : data.toExponential(4);}"))), 
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
                                             stateSave = TRUE)))



# ---------------------------------------------------------------
# ---------------------------------------------------------------
# MARKER GENES - HEATMAP
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr heterogeneity_markerGenes_heatmap

# DoHeatmap replaced by use of iHeatmapr after testing several options
#htmltools::tagList( ggplotly( DoHeatmap(object = sc10x.rna.seurat.ad.all, features = topMarkersDF[["gene"]])));

# Get the matrix of expression and associated clusters from Seurat object 
expMat = as.matrix( GetAssayData( sc10x.rna.seurat.ad.all));
clusterID = Idents( sc10x.rna.seurat.ad.all);

# Select marker genes and reorder cells to group clusters together
topMarkersGenes = topMarkersDF[["gene"]];
topMarkersHalfGenes = topMarkersHalfDF[["gene"]];
clusterOrdering = order( clusterID);

topMarkersGenes_expMat = expMat[topMarkersGenes, clusterOrdering];
topMarkersHalfGenes_expMat = expMat[topMarkersHalfGenes, clusterOrdering];
clusterID = clusterID[clusterOrdering];

# Create a matrix for text on mouse over events (customize 'text' that is supposed to show value for adding info)
hoverTextMarkers = topMarkersGenes_expMat;
hoverTextMarkers[] = paste( paste( "Cell ID:", sc10x.rna.seurat.ad.all[["numID", drop = FALSE]][colnames( topMarkersGenes_expMat)[col( topMarkersGenes_expMat)],]),
                            paste( "Value:", topMarkersGenes_expMat),
                            paste( "Cluster:", sort( clusterID)[col( topMarkersGenes_expMat)]),
                            paste( "Markers cl.:", topMarkersDF[["cluster"]][row( topMarkersGenes_expMat)]), 
                            sep = "<br>");


# Prepare rows and columns annotation bars (cluster groups)
rowsAnnot = data.frame(Cluster = topMarkersDF[["cluster"]]);
rowsAnnotHalf = data.frame(Cluster = topMarkersHalfDF[["cluster"]]);
colsAnnot = data.frame(Cluster = clusterID);

# Create heatmap
heatmapPlot = iheatmap(topMarkersGenes_expMat, 
                       row_labels = TRUE, 
                       col_labels = FALSE, 
                       text = hoverTextMarkers,
                       tooltip = setup_tooltip_options(row = TRUE,
                                                       col = TRUE,
                                                       value = TRUE,
                                                       prepend_row = "Gene: ",
                                                       prepend_col = "Cell: ",
                                                       prepend_value = ""),
                       row_annotation = rowsAnnot, 
                       col_annotation = colsAnnot,
                       row_annotation_colors = list(Cluster = hue_pal()(nlevels(Idents( sc10x.rna.seurat.ad.all)))),
                       col_annotation_colors = list(Cluster = hue_pal()(nlevels(Idents( sc10x.rna.seurat.ad.all)))));

# Hack to customize modebar buttons as for plotly native objects
# see: https://github.com/ropensci/iheatmapr/issues/38#issuecomment-445428166
heatmapPlot_widget = iheatmapr::to_widget(heatmapPlot);
# Edit the htmlwidget object itself
heatmapPlot_widget[["x"]][["config"]][["displaylogo"]] = FALSE;
heatmapPlot_widget[["x"]][["config"]][["toImageButtonOptions"]] = list( format='svg');  # This one does not seem to work... TODO: Check issue on GitHub
heatmapPlot_widget[["x"]][["config"]][["modeBarButtons"]]       = list( list( 'toImage'),
                                                                        list( 'zoom2d', 'zoomIn2d', 'zoomOut2d', 'pan2d', 'resetScale2d'));

# # Render the heatmap
# htmltools::tagList(heatmapPlot_widget)



# ---------------------------------------------------------------
# ---------------------------------------------------------------
# MARKER GENES - EXPRESSION on T-SNE and UMAP
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr heterogeneity_markerGenes_expression
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

# Plot expression values of marker genes (TODO: message if list empty)
resNULL=lapply(names(topMarkers), function(clusterName)
{
  cat("###### Cl.", clusterName, " \n");
  
  topMarkersGenes = topMarkers[[clusterName]][["gene"]]
  
  # Extract cells from current cluster to highlight them on a dimreduc plot
  clusterCells = which( Idents( sc10x.rna.seurat.ad.all) == clusterName);
  # Create a dimreduc plot with current cluster highlighted (useful for large number of clusters)
  print( DimPlot( sc10x.rna.seurat.ad.all, 
                  reduction = ifelse(exists("useReduction"), useReduction, "umap"), 
                  cols="#44444422", 
                  cells.highlight = clusterCells, 
                  cols.highlight = "#FF000088", 
                  sizes.highlight = 1.5, 
                  order = clusterCells,  # Plot highlighted cells last
                  group.by=NULL) + 
           theme( axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  legend.position = "none") +
           ggtitle( paste( "Cluster", clusterName)));
  
  # Prepare panels as list of ggplots and remove axes title
  markersPanels = lapply( topMarkersGenes, function(featuresName) 
  { 
    FeaturePlot( sc10x.rna.seurat.ad.all, features = featuresName, reduction = ifelse(exists("useReduction"), useReduction, "umap"), order = TRUE, cols = col_greyMagma) +
      theme( axis.title.x = element_blank(),
             axis.title.y = element_blank(),
             legend.position = "none");
    
  });
  
  resNULL = lapply( markersPanels, print);
  cat(" \n \n"); # Required for '.tabset'
  
  
});

# ---------------------------------------------------------------
# ---------------------------------------------------------------
# MARKER GENES - EXPRESSION ON VIOLIN PLOTS
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr heterogeneity_markerGenes_expression_violinplot

topMarkersGenes= sort( unique( do.call( c, lapply(names(topMarkers), function( clusterName){
  return( topMarkers[[clusterName]][["gene"]])
}))))

# Gather data to be visualized together (cell name + numID + metrics + Cluster)
cellsDataMarkers = cbind( "Cell" = colnames( sc10x.rna.seurat.ad.all), # cell names from rownames conflict with search field in datatable
                          sc10x.rna.seurat.ad.all[[c( "numID", "nCount_RNA", "nFeature_RNA")]],
                          "percent.mito" = as.numeric(format(sc10x.rna.seurat.ad.all[["percent.mito", drop = TRUE]], digits = 5)),
                          "Cluster" = Idents( sc10x.rna.seurat.ad.all));

# # Add gene expression information on cell data
original_names = names( cellsDataMarkers)
for( marker_gene in topMarkersGenes){
  cellsDataMarkers = cbind( cellsDataMarkers, sc10x.rna.seurat.ad.all@assays$RNA@data[ marker_gene, colnames( sc10x.rna.seurat.ad.all)])
}
names( cellsDataMarkers) = c( original_names, topMarkersGenes)

# Generate plotly violin/jitter panels for each marker gene
for( marker_gene in topMarkersGenes){
  
  # Compute the maximal expression value of the marker genen over the clusters
  max_expression = max( cellsDataMarkers[ , marker_gene])
  
  # Compute the number of cells with positive expression in the different clusters
  positive_cellsDataMarkers = cellsDataMarkers[ which( cellsDataMarkers[ , marker_gene] > 0), "Cluster"]
  positive_cellsDataMarkers = as.data.frame( table( positive_cellsDataMarkers))
  names( positive_cellsDataMarkers) = c( "Cluster", "NbCells")
  
  # Compute the number of cells with null expression in the different clusters
  null_cellsDataMarkers = cellsDataMarkers[ which( cellsDataMarkers[ , marker_gene] == 0), "Cluster"]
  null_cellsDataMarkers = as.data.frame( table( null_cellsDataMarkers))
  names( null_cellsDataMarkers) = c( "Cluster", "NbCells")
  
  # Plot the expression of the marker gene for each cluster separatly as violin + jitter
  # and adding the number of positive expressed cells and null expressed cells above each cluster plot
  print( 
    ggplot( cellsDataMarkers, aes_string ( x="Cluster", y=as.name( marker_gene))) +
      geom_violin( aes( col = Cluster, fill=Cluster)) +
      scale_color_manual(values=col_rainbow) +
      scale_fill_manual(values=col_rainbow) +
      geom_jitter( width=0.2, size = 0.2, shape = 1) + 
      geom_text( data = positive_cellsDataMarkers, aes( x=Cluster, y = 1.15*max_expression, label=NbCells), size = 3) +
      geom_text( data = null_cellsDataMarkers, aes( x=Cluster, y = 1.1*max_expression, label=NbCells), size = 3) +
      theme_minimal() +
      theme( legend.position = "None")
  )
}



# ---------------------------------------------------------------
# ---------------------------------------------------------------
# MONITORED GENES
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr heterogeneity_monitoredGenes

# Just remind the warning for genes names not in object
if(any( is.na( MARKERS_ALL))) warning( paste( "Following gene(s) to be monitored could not be found in experimental data and will be ignored:", paste( monitoredGenesNotFound, collapse=" - ")));




## @knitr heterogeneity_monitoredGenes_heatmap

# DoHeatmap replaced by use of iHeatmapr after testing several options
#htmltools::tagList( ggplotly( DoHeatmap(object = sc10x.rna.seurat.ad.all, features = unlist(MARKERS_PVM))));

# Get the matrix of expression and associated clusters from Seurat object 
expMat = as.matrix( GetAssayData( sc10x.rna.seurat.ad.all));
clusterID = Idents( sc10x.rna.seurat.ad.all);

# Select monitored genes and reorder cells to group clusters together
monitoredGenes = unlist( MARKERS_PVM);
clusterOrdering = order( clusterID);

expMat = expMat[monitoredGenes, clusterOrdering];
clusterID = clusterID[ clusterOrdering];

# Create a matrix for text on mouse over events (customize 'text' that is supposed to show value for adding info)
hoverTextMarkers = expMat;
hoverTextMarkers[] = paste( paste( "Cell ID:", sc10x.rna.seurat.ad.all[["numID", drop = FALSE]][colnames( expMat)[col( expMat)],]),
                            paste( "Value:", expMat),
                            paste( "Cluster:", sort( clusterID)[col( expMat)]),
                            paste( "Monitored Grp.:", rep(names(MARKERS_PVM), sapply(MARKERS_PVM, length))[row( expMat)]), 
                            sep = "<br>");

# Prepare rows and columns annotation bars (monitored group and cluster respectively)
rowsAnnot = data.frame( Monitored = rep( names( MARKERS_PVM), sapply( MARKERS_PVM, length)));
colsAnnot = data.frame( Cluster = clusterID);

# Create heatmap
heatmapPlot = iheatmap( expMat,
                        row_labels = TRUE,
                        col_labels = FALSE,
                        text = hoverTextMarkers,
                        tooltip = setup_tooltip_options(row = TRUE,
                                                        col = TRUE,
                                                        value = TRUE,
                                                        prepend_row = "Gene: ",
                                                        prepend_col = "Cell: ",
                                                        prepend_value = ""),
                        row_annotation = rowsAnnot,
                        col_annotation = colsAnnot,
                        col_annotation_colors = list( Cluster = hue_pal()( nlevels( Idents( sc10x.rna.seurat.ad.all)))));

# Hack to customize modebar buttons as for plotly native objects
# see: https://github.com/ropensci/iheatmapr/issues/38#issuecomment-445428166
heatmapPlot_widget = iheatmapr::to_widget(heatmapPlot);
# Edit the htmlwidget object itself
heatmapPlot_widget[["x"]][["config"]][["displaylogo"]] = FALSE;
heatmapPlot_widget[["x"]][["config"]][["toImageButtonOptions"]] = list( format='svg');  # This one does not seem to work... TODO: Check issue on GitHub
heatmapPlot_widget[["x"]][["config"]][["modeBarButtons"]]       = list( list( 'toImage'),
                                                                        list( 'zoom2d', 'zoomIn2d', 'zoomOut2d', 'pan2d', 'resetScale2d'));

# # Render the heatmap
# htmltools::tagList(heatmapPlot_widget)


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# MONITORED GENES - EXPRESSION - on T-SNE AND UMAP
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr heterogeneity_monitoredGenes_expression
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

MARKERS_ALL <- list( ALL = MARKERS_ALL)

# Plot expression values of monitored genes (TODO: message if list empty)
resNULL = lapply( names( MARKERS_ALL), function(monitoredGroup)
{
  cat("######", monitoredGroup, " \n");
  
  # List of plots (or error message if feature not found)
  expressionPlots = lapply( MARKERS_ALL[[monitoredGroup]], function(featuresName)
  {
    tryCatch( FeaturePlot( sc10x.rna.seurat.ad.all, features = featuresName, reduction = ifelse(exists("useReduction"), useReduction, "umap"), order = TRUE, cols = col_greyMagma)  +
                theme( axis.title.x = element_blank(),
                       axis.title.y = element_blank(),
                       legend.position = "none"),
              error=function(e){return(conditionMessage(e))}); # In case gene name could not be found in object (should have been removed earlier anyway)
  });
  
  resNULL = lapply( expressionPlots, print);
  cat(" \n \n"); # Required for '.tabset'
});


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# MONITORED GENES - EXPRESSION ON VIOLIN PLOTS
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr heterogeneity_monitoredGenes_expression_violinplot

# Gather data to be visualized together (cell name + numID + metrics + Cluster)
cellsDataMonitored = cbind( "Cell" = colnames( sc10x.rna.seurat.ad.all), # cell names from rownames conflict with search field in datatable
                            sc10x.rna.seurat.ad.all[[c( "numID", "nCount_RNA", "nFeature_RNA")]],
                            "percent.mito" = as.numeric(format(sc10x.rna.seurat.ad.all[["percent.mito", drop = TRUE]], digits = 5)),
                            "Cluster" = Idents( sc10x.rna.seurat.ad.all));

all_monitored_genes = sort( unique( unlist( MARKERS_PVM, use.names = FALSE)))

# # Add gene expression information on cell data
original_names = names( cellsDataMonitored)
for( monitored_gene in all_monitored_genes){
  cellsDataMonitored = cbind( cellsDataMonitored, sc10x.rna.seurat.ad.all@assays$RNA@data[ monitored_gene, colnames( sc10x.rna.seurat.ad.all)])
}
names( cellsDataMonitored) = c( original_names, all_monitored_genes)

# Generate plotly violin/jitter panels for each monitored gene
for( monitored_gene in all_monitored_genes){
  
  # Compute the maximal expression value of the monitored genen over the clusters
  max_expression = max( cellsDataMonitored[ , monitored_gene])
  
  # Compute the number of cells with positive expression in the different clusters
  positive_cellsDataMonitored = cellsDataMonitored[ which( cellsDataMonitored[ , monitored_gene] > 0), "Cluster"]
  positive_cellsDataMonitored = as.data.frame( table( positive_cellsDataMonitored))
  names( positive_cellsDataMonitored) = c( "Cluster", "NbCells")
  
  # Compute the number of cells with null expression in the different clusters
  null_cellsDataMonitored = cellsDataMonitored[ which( cellsDataMonitored[ , monitored_gene] == 0), "Cluster"]
  null_cellsDataMonitored = as.data.frame( table( null_cellsDataMonitored))
  names( null_cellsDataMonitored) = c( "Cluster", "NbCells")
  
  # Plot the expression of the monitored gene for each cluster separatly as violin + jitter
  # and adding the number of positive expressed cells and null expressed cells above each cluster plot
  print( 
    ggplot( cellsDataMonitored, aes_string ( x="Cluster", y=as.name( monitored_gene))) +
      geom_violin( aes( col = Cluster, fill=Cluster)) +
      scale_color_manual(values=col_rainbow) +
      scale_fill_manual(values=col_rainbow) +
      geom_jitter( width=0.2, size = 0.2, shape = 1) + 
      geom_text( data = positive_cellsDataMonitored, aes( x=Cluster, y = 1.15*max_expression, label=NbCells), size = 3) +
      geom_text( data = null_cellsDataMonitored, aes( x=Cluster, y = 1.1*max_expression, label=NbCells), size = 3) +
      theme_minimal() +
      theme( legend.position = "None")
  )
}


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# MONITORED GENES - EXPRESSION ON DOTPLOT
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr heterogeneity_monitoredGenes_expression_dotplot

monitored_genes_dist <- as.dist(1 - cor(t(as.matrix(GetAssayData(sc10x.rna.seurat.ad.all, slot = "data")[unlist(MARKERS_ALL),], method = "pearson"))))
monitored_genes_clust <- hclust(monitored_genes_dist, method = "average")
monitored_genes_ordered <- MARKERS_ALL$ALL[monitored_genes_clust$order]

DotPlot(sc10x.rna.seurat.ad.all, features = monitored_genes_ordered, cols = c("blue", "orange", "red"), dot.scale = 8) +
  RotatedAxis() +
  coord_flip()

DotPlot(sc10x.rna.seurat.ad.all, features = monitored_genes_ordered, cols = c("blue", "orange", "red"), dot.scale = 8, cluster.idents = TRUE) +
  RotatedAxis() +
  coord_flip()




# Density plot for monitored genes
## @knitr heterogeneity_monitoredGenes_densityplot

for (i in MARKERS_ALL$ALL) {
  print(plot_density(object = sc10x.rna.seurat.ad.all, features = i, reduction = "umap", size = 0.5))
}




#######################
# Save combined seurat object
#######################

## @knitr save_data
# saveRDS(sc10x.rna.seurat.ad.all, file = file.path(PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "sc10x.rna.seurat.ad.all.RDS")))
# print(paste0("Exporting combined seurat object to : ", file.path(PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "sc10x.rna.seurat.ad.all.RDS"))))
# 



