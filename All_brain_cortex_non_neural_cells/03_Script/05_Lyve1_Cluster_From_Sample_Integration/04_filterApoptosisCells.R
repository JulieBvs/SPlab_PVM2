# ########################################################
# This script aims to remove bad quality cells not filtered in QC
# ########################################################




# Identify apoptotic cells

## @knitr removeApoptoticCells

FeaturePlot(sc10x.rna.seurat, cols = col_greyMagma, reduction = 'umap', features = "Malat1")

VlnPlot(sc10x.rna.seurat, features = "percent.ribo")
VlnPlot(sc10x.rna.seurat, features = "percent.mito")



test <- sc10x.rna.seurat[[]][,c("percent.ribo", "percent.mito", "integrated_snn_res.0.5")]

for (i in unique(test[,"integrated_snn_res.0.5"])) {
  print(ggplot(test, aes(x = percent.ribo, y = percent.mito, color = integrated_snn_res.0.5)) +
          geom_point(size = 0.7) +
          scale_color_manual(values = c(rep("grey", i), "red", rep("grey", max(as.numeric(test[,"integrated_snn_res.0.5"])) - 1 ))) +
          #  scale_color_manual(values = col_rainbow) +
          theme_light()
  )
  
}


cl_apo <- suppressWarnings(which.min(AverageExpression(sc10x.rna.seurat, features = "Malat1")$RNA) - 1)
cell_bc_apo <- names(Idents(sc10x.rna.seurat)[which(Idents(sc10x.rna.seurat) == cl_apo)])

apo_cells_outfile = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "barcodes_apoptosis_cells.txt"))
write.table( cell_bc_apo, file = apo_cells_outfile, quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
cat("<BR><BR>Writing barcodes of apoptotic cells to file: <BR>", apo_cells_outfile)




# Identify contaminations

## @knitr removeContaminations

VlnPlot(sc10x.rna.seurat, features = "nFeature_RNA", group.by = "integrated_snn_res.0.5")
VlnPlot(sc10x.rna.seurat, features = "nCount_RNA", group.by = "integrated_snn_res.0.5")




cl_cont <- suppressWarnings(which.max(AverageExpression(sc10x.rna.seurat, features = "Ly6c1")$RNA) - 1)
cell_bc_cont <- names(Idents(sc10x.rna.seurat)[which(Idents(sc10x.rna.seurat) == cl_cont)])

cont_cells_outfile = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "barcodes_contamination_cells.txt"))
write.table( cell_bc_cont, file = cont_cells_outfile, quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
cat("<BR><BR>Writing barcodes of contamination cells to file: <BR>", cont_cells_outfile)








##Subset seurat object
cells_bc <- Cells(sc10x.rna.seurat)[!Cells(sc10x.rna.seurat) %in% cell_bc_apo]
cells_bc <- cells_bc[!cells_bc %in% cell_bc_cont]
sc10x.rna.seurat.subset <- subset(x = sc10x.rna.seurat, cells = cells_bc)
saveRDS( sc10x.rna.seurat.subset, file = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "withoutApoptoticAndContaminationCells_seurat_object.RDS")))

