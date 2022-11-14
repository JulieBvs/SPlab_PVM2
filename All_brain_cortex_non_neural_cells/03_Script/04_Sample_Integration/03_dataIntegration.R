# ===========================================
# This script combine the seurat object from the 
# replicates in one combined seurat object.
# ===========================================



#######################
# Data integration using seurat method
#######################

## @knitr dataIntegration

# Add metadata from individual timepoint analysis
sc10x.rna.seurat.ju[["sample"]] <- "juvenile"
sc10x.rna.seurat.ad[["sample"]] <- "adult"
sc10x.rna.seurat.ag[["sample"]] <- "aged"

sc10x.rna.seurat.ju[["ju_cl"]] <- sc10x.rna.seurat.ju[["seurat_clusters"]]
sc10x.rna.seurat.ad[["ad_cl"]] <- sc10x.rna.seurat.ad[["seurat_clusters"]]
sc10x.rna.seurat.ag[["ag_cl"]] <- sc10x.rna.seurat.ag[["seurat_clusters"]]


# Create a list of seurat object containing the replicates
test <- list(juvenile = sc10x.rna.seurat.ju,
             adult = sc10x.rna.seurat.ad,
             aged = sc10x.rna.seurat.ag)

# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = test)

# Find anchors for data integration
anchors <- suppressMessages(FindIntegrationAnchors(object.list = test, anchor.features = features, verbose = TRUE))

# Combine replicate in one integrated seurat object
sc10x.rna.seurat.combined <- IntegrateData(anchorset = anchors, verbose = .VERBOSE)




#######################
# Integration analysis
#######################

## @knitr analysisIntegration

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(sc10x.rna.seurat.combined) <- "integrated"

# Run the standard workflow for visualization
sc10x.rna.seurat.combined <- ScaleData(sc10x.rna.seurat.combined, verbose = FALSE)
sc10x.rna.seurat.combined <- RunPCA(sc10x.rna.seurat.combined, npcs = N_PCA, verbose = FALSE)
sc10x.rna.seurat.combined <- RunUMAP(sc10x.rna.seurat.combined, reduction = "pca", dims = 1:N_PCA, verbose = .VERBOSE)
sc10x.rna.seurat.combined <- RunTSNE(sc10x.rna.seurat.combined, reduction = "pca", dims = 1:N_PCA, verbose = .VERBOSE)



## @knitr visualizationIntegration

# Plot the projection with relicate colors
print( suppressWarnings( DimPlot( sc10x.rna.seurat.combined, reduction = ifelse(exists("useReduction"), useReduction, "umap"), group.by = "sample")) + 
         theme( axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                plot.margin = margin(0, 0, 0, 0, "cm")))

print( suppressWarnings( DimPlot( sc10x.rna.seurat.combined, reduction = ifelse(exists("useReduction"), useReduction, "umap"), group.by = "sample", split.by = "sample")) + 
         theme( axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                legend.position = "none", 
                plot.margin = margin(0, 0, 0, 0, "cm")) +
         coord_fixed(1))
cat(" \n \n") # Required for '.tabset'





#######################
# Cell clustering from individual analysis (each time-point)
#######################

### Juvenile ###
## @knitr juvenile_clustering_umap
DimPlot(sc10x.rna.seurat.combined, reduction = 'umap', group.by = 'ju_cl', cols = col_rainbow, na.value = "lightgrey", order = TRUE)

## @knitr juvenile_clustering_umap_splitted
juv_cl <- unique(sc10x.rna.seurat.ju[[]][,"ju_cl"])
juv_cl <- juv_cl[order(juv_cl)]
for (cl in juv_cl) {
  ju_cl <- rownames(sc10x.rna.seurat.combined[[]][which(sc10x.rna.seurat.combined[["ju_cl"]] == cl),])
  print(DimPlot(sc10x.rna.seurat.combined, reduction = 'umap', cells.highlight = ju_cl) + 
          ggtitle(paste0("Cluster ", cl)) + 
          NoLegend())
}


### Adult ###
## @knitr adult_clustering_umap
DimPlot(sc10x.rna.seurat.combined, reduction = 'umap', group.by = 'ad_cl', cols = col_rainbow, na.value = "lightgrey", order = TRUE)

## @knitr adult_clustering_umap_splitted
adu_cl <- unique(sc10x.rna.seurat.ad[[]][,"ad_cl"])
adu_cl <- adu_cl[order(adu_cl)]
for (cl in adu_cl) {
  ad_cl <- rownames(sc10x.rna.seurat.combined[[]][which(sc10x.rna.seurat.combined[["ad_cl"]] == cl),])
  print(DimPlot(sc10x.rna.seurat.combined, reduction = 'umap', cells.highlight = ad_cl) + 
          ggtitle(paste0("Cluster ", cl)) + 
          NoLegend())
}


### Aged ###
## @knitr aged_clustering_umap
DimPlot(sc10x.rna.seurat.combined, reduction = 'umap', group.by = 'ag_cl', cols = col_rainbow, na.value = "lightgrey", order = TRUE)

## @knitr aged_clustering_umap_splitted
age_cl <- unique(sc10x.rna.seurat.ag[[]][,"ag_cl"])
age_cl <- age_cl[order(age_cl)]
for (cl in age_cl) {
  ag_cl <- rownames(sc10x.rna.seurat.combined[[]][which(sc10x.rna.seurat.combined[["ag_cl"]] == cl),])
  print(DimPlot(sc10x.rna.seurat.combined, reduction = 'umap', cells.highlight = ag_cl) + 
          ggtitle(paste0("Cluster ", cl)) + 
          NoLegend())
}

cat(" \n \n") # Required for '.tabset'






#######################
# Cell clustering
#######################

## @knitr cellClustering
sc10x.rna.seurat.combined <- FindNeighbors(sc10x.rna.seurat.combined, reduction = "pca", dims = 1:N_PCA)
sc10x.rna.seurat.combined <- FindClusters(sc10x.rna.seurat.combined, resolution = CL_RESOLUTION)


# Plot the projection with colored clusters
## @knitr cellClustering_visualization
print( suppressWarnings( DimPlot( sc10x.rna.seurat.combined, reduction = ifelse(exists("useReduction"), useReduction, "umap"), 
                                  label = FALSE, cols = c(col_rainbow, col_rainbow))) + 
         theme( axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                plot.margin = margin(0, 0, 0, 0, "cm")))

print( suppressWarnings( DimPlot( sc10x.rna.seurat.combined, reduction = ifelse(exists("useReduction"), useReduction, "umap"), split.by = "sample", cols = col_rainbow[1:nlevels( Idents( sc10x.rna.seurat.combined))])) + 
         theme( axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                legend.position = "none", 
                plot.margin = margin(0, 0, 0, 0, "cm")) +
         coord_fixed(1))
cat(" \n \n") # Required for '.tabset'


pvm <- rownames(sc10x.rna.seurat.combined[[]][which(sc10x.rna.seurat.combined[["integrated_snn_res.0.4"]] == 7),])
bam <- rownames(sc10x.rna.seurat.combined[[]][which(sc10x.rna.seurat.combined[["integrated_snn_res.0.4"]] == 13),])
microglia <- rownames(sc10x.rna.seurat.combined[[]][which(sc10x.rna.seurat.combined[["integrated_snn_res.0.4"]] == 4 |
                                                            sc10x.rna.seurat.combined[["integrated_snn_res.0.4"]] == 1 |
                                                            sc10x.rna.seurat.combined[["integrated_snn_res.0.4"]] == 15 |
                                                            sc10x.rna.seurat.combined[["integrated_snn_res.0.4"]] == 23),])
sc10x.rna.seurat.combined[["simplified_clusters"]] <- 0
sc10x.rna.seurat.combined[["simplified_clusters"]][bam,] <- 2
sc10x.rna.seurat.combined[["simplified_clusters"]][pvm,] <- 1
sc10x.rna.seurat.combined[["simplified_clusters"]][microglia,] <- 3



devoff_msg <- dev.off()
png(filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "UMAP_pvm_bam_microglia.png")), width = 900, height = 900)
print( suppressWarnings( DimPlot( sc10x.rna.seurat.combined, reduction = ifelse(exists("useReduction"), useReduction, "umap"), group.by = "simplified_clusters",
                                  label = FALSE, cols = c("lightgrey", col_rainbow[c(14,11,9)]))) + 
         theme( axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                plot.margin = margin(0, 0, 0, 0, "cm")))
devoff_msg <- dev.off()


#######################
# Identify conserved cell type markers
#######################

## @knitr heterogeneity_conservedmarkerGenes

DefaultAssay(sc10x.rna.seurat.combined) <- "RNA"
#BiocManager::install('multtest')
#install.packages('metap')

######## Get top 10 of conserved markers
markers <- data.frame()
for (i in 0:(length(unique(Idents(sc10x.rna.seurat.combined)))-1)) {
  markers_temp <- FindConservedMarkers(sc10x.rna.seurat.combined, ident.1 = i, grouping.var = "sample", verbose = FALSE)[1:FINDMARKERS_SHOWTOP,] #Run with seurat4.1.0
  markers_temp[,"cluster"] <- i
  markers_temp[,"gene"] <- rownames(markers_temp)
  markers <- rbind.fill(markers, markers_temp)
}

#markers_save <- markers

topMarkersDT = markers[c("gene", "cluster", "minimump_p_val", "juvenile_p_val_adj", "aged_p_val_adj", "adult_p_val_adj")]
topMarkersDT[,"minimump_p_val"] <- as.character(topMarkersDT[,"minimump_p_val"])
topMarkersDT[,"juvenile_p_val_adj"] <- as.character(topMarkersDT[,"juvenile_p_val_adj"])
topMarkersDT[,"aged_p_val_adj"] <- as.character(topMarkersDT[,"aged_p_val_adj"])
topMarkersDT[,"adult_p_val_adj"] <- as.character(topMarkersDT[,"adult_p_val_adj"])

topMarkersDT[is.na(topMarkersDT)] <- "no value"


######## Table of top 10 conserved markers
## @knitr heterogeneity_conservedmarkerGenes_table
datatable( topMarkersDT, 
           class = "compact",
           filter="top",
           rownames = FALSE,
           colnames = c("Gene", "Cluster", "Minimum Adj. Pvalue", "Rep1 - Adj. Pvalue", "Rep2 - Adj. Pvalue", "Rep3 - Adj. Pvalue"),
           caption = paste(ifelse( is.null( FINDMARKERS_SHOWTOP), "All", paste("Top", FINDMARKERS_SHOWTOP)), "marker genes for each cluster"),
           extensions = c('Buttons', 'Select'),
           options = list(dom = "<'row'<'col-sm-8'B><'col-sm-4'f>> <'row'<'col-sm-12'l>> <'row'<'col-sm-12'rt>> <'row'<'col-sm-12'ip>>", # Set elements for CSS formatting ('<Blf><rt><ip>')
                          autoWidth = FALSE,
                          buttons = exportButtonsListDT,
                          lengthMenu = list(c( 20, 50, 100, -1),
                                            c( 20, 50, 100, "All")),
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



######## Heatmap of conserved markers
## @knitr heterogeneity_conservedmarkerGenes_heatmap

sc10x.rna.seurat.combined <- ScaleData(sc10x.rna.seurat.combined, features = rownames(sc10x.rna.seurat.combined))
DoHeatmap(sc10x.rna.seurat.combined, features = topMarkersDT$gene, group.colors = col_rainbow[1:nlevels( Idents( sc10x.rna.seurat.combined))]) + NoLegend()




#######################
# Identify gene markers
#######################

## @knitr heterogeneity_markerGenes

DefaultAssay(sc10x.rna.seurat.combined) <- "RNA"
#BiocManager::install('multtest')
#install.packages('metap')

######## Get top 10 of conserved markers
#markers <- data.frame()

markers <- FindAllMarkers(object          = sc10x.rna.seurat.combined,
                          test.use        = "wilcox",
                          only.pos        = TRUE,
                          min.pct         = 0.1,
                          logfc.threshold = 0.25)
markers <- markers[,c("gene", colnames(markers)[1:6])]
write.table(markers, file = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "all_gene_markers.txt")), sep = "\t", quote = FALSE, dec = ".", row.names = FALSE)

topMarkers = by( markers, markers[["cluster"]], function(x)
{
  # Filter markers based on adjusted PValue
  x = x[ x[["p_val_adj"]] < 0.001, , drop = FALSE];
  # Sort by decreasing logFC
  x = x[ order(abs(x[["avg_logFC"]]), decreasing = TRUE), , drop = FALSE ]
  # Return top ones
  return( if(is.null( 30)) x else head( x, n = 30));
});

# Filter markers by cluster (TODO: check if downstream code works when no markers found)
topMarkersHalf = by( markers, markers[["cluster"]], function(x)
{
  # Filter markers based on adjusted PValue
  x = x[ x[["p_val_adj"]] < 0.001, , drop = FALSE];
  # Sort by decreasing logFC
  x = x[ order(abs(x[["avg_logFC"]]), decreasing = TRUE), , drop = FALSE ]
  # Return top ones
  return( if(is.null( 30)) x else head( x, n = floor( 30 /2)));
});

# Merge marker genes in a single data.frame and render it as datatable
topMarkersDF = do.call( rbind, topMarkers);
topMarkersHalfDF = do.call( rbind, topMarkersHalf);
# Select and order columns to be shown in datatable
topMarkersDT = topMarkersDF[c("gene", "cluster", "avg_logFC", "p_val_adj")]
topMarkersHalfDT = topMarkersHalfDF[c("gene", "cluster", "avg_logFC", "p_val_adj")]



######## Table of top 10 conserved markers
## @knitr heterogeneity_markerGenes_table
datatable( topMarkersDT, 
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
                          stateSave = TRUE));



######## Heatmap of conserved markers
## @knitr heterogeneity_markerGenes_heatmap

sc10x.rna.seurat.combined <- ScaleData(sc10x.rna.seurat.combined, features = rownames(sc10x.rna.seurat.combined))
DoHeatmap(sc10x.rna.seurat.combined, features = topMarkersDT$gene, group.colors = col_rainbow[1:nlevels( Idents( sc10x.rna.seurat.combined))]) + NoLegend()







#######################
# DotPlot by samples
#######################

## @knitr integartion_dotplotbysamples


#######################
# Marker genes
#######################

## @knitr heterogeneity_markerGenes





#######################
# Monitored genes
#######################

## @knitr heterogeneity_monitoredGenes_dotplot

# DotPlot for monitored genes
DotPlot(sc10x.rna.seurat.combined, features = MARKERS_ALL, cols = c("blue", "orange", "red"), dot.scale = 8) +
  RotatedAxis()
DotPlot(sc10x.rna.seurat.combined, features = MARKERS_ALL, cols = c("blue", "orange", "red"), dot.scale = 8) +
  RotatedAxis() +
  coord_flip()

#Same by sample
DotPlot(sc10x.rna.seurat.combined, features = MARKERS_ALL, cols = c("blue", "orange", "red"), dot.scale = 8, split.by = "sample") +
  RotatedAxis()
DotPlot(sc10x.rna.seurat.combined, features = MARKERS_ALL, cols = c("blue", "orange", "red"), dot.scale = 8, split.by = "sample") +
  RotatedAxis() + 
  coord_flip()


devoff_msg <- dev.off()
png(filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "dotplot1_pvm_bam_microglia.png")), width = 1200, height = 300)
DotPlot(sc10x.rna.seurat.combined, features = MARKERS_ALL, cols = c("blue", "orange", "red"), dot.scale = 15, group.by = "simplified_clusters") +
  RotatedAxis() +
  scale_y_discrete(limits=rev)
devoff_msg <- dev.off()

devoff_msg <- dev.off()
png(filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "dotplot2_pvm_bam_microglia.png")), width = 400, height = 1200)
DotPlot(sc10x.rna.seurat.combined, features = MARKERS_ALL, cols = c("blue", "orange", "red"), dot.scale = 15, group.by = "simplified_clusters") +
  RotatedAxis() +
  coord_flip()
devoff_msg <- dev.off()


# FeaturePlot for monitored genes
## @knitr heterogeneity_monitoredGenes_tsne_umap

#FeaturePlot(sc10x.rna.seurat.combined, features = "Itgam", cols = col_comet, order = TRUE)
for (i in MARKERS_ALL) {
  print(FeaturePlot(sc10x.rna.seurat.combined, features = i, cols = col_greyMagma, order = TRUE))
}

# Save FeaturePlot for markers of proliferation
devoff_msg <- dev.off()
png(filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "umap_prolif_Mki67.png")), width = 600, height = 600)
FeaturePlot(sc10x.rna.seurat.combined, features = "Mki67", cols = col_greyMagma, order = TRUE)
devoff_msg <- dev.off()

devoff_msg <- dev.off()
png(filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "umap_prolif_Birc5.png")), width = 600, height = 600)
FeaturePlot(sc10x.rna.seurat.combined, features = "Birc5", cols = col_greyMagma, order = TRUE)
devoff_msg <- dev.off()

devoff_msg <- dev.off()
png(filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "umap_prolif_Top2a.png")), width = 600, height = 600)
FeaturePlot(sc10x.rna.seurat.combined, features = "Top2a", cols = col_greyMagma, order = TRUE)
devoff_msg <- dev.off()

# Violin plot for monitored genes
## @knitr heterogeneity_monitoredGenes_vlnplot
for (i in MARKERS_ALL) {
  suppressMessages(suppressWarnings(print(VlnPlot(sc10x.rna.seurat.combined, features = i, pt.size = 0.1, cols = col_rainbow[1:nlevels( Idents( sc10x.rna.seurat.combined))]))) + 
                     theme( legend.position = "none"))
}

# Density plot for monitored genes
## @knitr heterogeneity_monitoredGenes_densityplot

for (i in MARKERS_ALL) {
  print(plot_density(object = sc10x.rna.seurat.combined, features = i, reduction = "umap", size = 0.5))
}




#######################
# Save combined seurat object
#######################

## @knitr save_data
saveRDS(sc10x.rna.seurat.combined, file = file.path(PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "sc10x.rna.seurat.integrated.RDS")))
print(paste0("Exporting combined seurat object to : ", file.path(PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "sc10x.rna.seurat.integrated.RDS"))))




