# #############################################
# In this script, we extract cell barcode and filtered 
# seurat object corresponding to Mrc1 positive cells
# #############################################



## @knitr extract_Lyve1_pos_cluster

#Extract cluster id with th hishest mean expression for Lyve1
cl_interest <- which.max(AverageExpression(sc10x.rna.seurat.combined, features = "Lyve1")$RNA) - 1

#Extract cell barcodes from cluster of interest
##Extract barcodes
cell_bc_interest <- names(Idents(sc10x.rna.seurat.combined)[which(Idents(sc10x.rna.seurat.combined) == cl_interest)])

##Save cell barcodes
filtered_cells_outfile = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "barcodes_cells_from_cl_interest_", cl_interest ,".txt"))
write.table( cell_bc_interest, file = filtered_cells_outfile, quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
cat("<BR><BR>Writing barcodes of cells from cluster", cl_interest, "to file: <BR>", filtered_cells_outfile)

##Subset seurat object
sc10x.rna.seurat.ad.subset <- subset(x = sc10x.rna.seurat.combined, idents = cl_interest)
#saveRDS( sc10x.rna.seurat.ad.subset, file = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "sc10x.rna.seurat.cl", cl_interest ,".RDS")))
cat("<BR><BR>Number of selected cells: <BR>", ncol(GetAssayData(sc10x.rna.seurat.ad.subset, slot = 'counts')))

