
library(stringr)
library(reticulate)
library(limma)
library(pheatmap)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(geneplotter)
library(genefilter)
library(RColorBrewer)
library(lattice)
library(latticeExtra)
library(grid)
library(Rtsne)
library(xlsx)
library(gridExtra)
library(dendsort)
library(viridis)


get_plot_dims <- function(heat_map)
{
  plot_height <- sum(sapply(heat_map$gtable$heights, grid::convertHeight, "in"))
  plot_width  <- sum(sapply(heat_map$gtable$widths, grid::convertWidth, "in"))
  dev.off()
  return(list(height = plot_height, width = plot_width))
}


process_CM <- function(data){
	b <- apply(data,2, FUN=function(x) (x/sum(x)))
	b[is.nan(b)] <- NA
	b[b==0] <- NA
	b <- b[,colSums(is.na(b))<nrow(b)]
	b <- b[rowSums(is.na(b))<ncol(b),]
	return(b)
}

mean_acc <- function(data_mat){
	return(round(mean(diag(data_mat[!grepl('Unassigned', rownames(data_mat)), ])),3))
}

pd <- import("pandas")
myBreaks <- c(seq(0,  0.4, length.out= 20), seq(0.41, 0.79, length.out=10), seq(0.8, 1, length.out=20))
color_red_green <- colorRampPalette(c('#4E62CC','#D8DBE2' , '#BA463E'))(50)
color_red_green <- colorRampPalette(c('#cb5b4c','#D8DBE2', '#1aab2d'))(50)

draw_cfmt <- function(dataSet, path = NULL, out_path = NULL){
  dir.create(out_path, showWarnings = FALSE)
  fontsize = 8
  fontsize_number = 8
  print(dataSet)
  switch(dataSet,
         human_dataset_random = {
           dataSet_name = 'Human Hematopoiesis'},
         mouse_atlas_random = {
           dataSet_name = 'Mouse Atlas'},
         pancreas_01 = {
           dataSet_name = 'Pancreas Bar16-Mur16'},
         pancreas_02 = {
           dataSet_name = 'Pancreas Bar16-Seg16'},
         human_blood_01 = {
           dataSet_name = 'PBMC 10x_v3-10x_v5'},
         stop("Does Not Exist!")
  )
  
  annotation <- pd$read_pickle(file.path(path, dataSet, 'JIND', 'JIND_assignmentbrftune.pkl'))
  annotation$cell_names <- rownames(annotation)
  
  annotation_seurat <- pd$read_pickle(file.path(path, dataSet, 'seurat', 'seurat_assignment.pkl'))
  annotation_seurat$cell_names <- rownames(annotation_seurat)
  
  data <- as.data.frame.matrix(table(annotation$raw_predictions, annotation$labels))
  b <- process_CM(data)
  b[is.na(b)] <- 0
  
  jind_raw <- pheatmap(b, color=color_red_green, cluster_rows = F, cluster_cols=F, display_numbers=T, cellheight=25, cellwidth=25, na_col= '#cc5a4e', main=paste0('JIND+ ', dataSet_name, '\n Mean Accuracy (per Cell Type): ', mean_acc(b)), legend=F, breaks = myBreaks, fontsize = fontsize, fontsize_number = fontsize_number)
  file_xlsx <- file.path(out_path, paste0(dataSet, '_CM.xlsx' ))
  if (file.exists(file_xlsx)) {
    file.remove(file_xlsx)
  }
  b[is.na(b)] <-0
  write.xlsx(b, file=file_xlsx, sheetName='JIND_raw', row.names = TRUE, append=TRUE)
  
  data <- as.data.frame.matrix(table(annotation$predictions, annotation$labels))
  b <- process_CM(data)
  b[is.na(b)] <- 0
  
  jind <- pheatmap(b, color=color_red_green, cluster_rows = F, cluster_cols=F, display_numbers=T, cellheight=25, cellwidth=25, na_col= '#cc5a4e', main=paste0('JIND+ ', dataSet_name, '\n Mean Accuracy (per Cell Type): ', mean_acc(b)), legend=F, breaks = myBreaks, fontsize = fontsize, fontsize_number = fontsize_number)
  b[is.na(b)] <-0
  write.xlsx(b, file=file_xlsx, sheetName='JIND', row.names = TRUE, append=TRUE)
  
  data_seurat <- as.data.frame.matrix(table(annotation_seurat$raw_predictions, annotation_seurat$labels))
  b_seurat <- process_CM(data_seurat)
  b_seurat[is.na(b_seurat)] <-0
  
  seurat_raw <- pheatmap(b_seurat, color=color_red_green, cluster_rows = F, cluster_cols=F, display_numbers=T, cellheight=25, cellwidth=25, na_col= '#cc5a4e', main=paste0('Seurat-LT ', dataSet_name, '\n Mean Accuracy (per Cell Type): ', mean_acc(b_seurat)), legend=T, breaks = myBreaks, fontsize = fontsize, fontsize_number = fontsize_number)
  write.xlsx(b_seurat, file=file_xlsx, sheetName='seurat', row.names = TRUE, append=TRUE)
  
  data_seurat <- as.data.frame.matrix(table(annotation_seurat$predictions, annotation_seurat$labels))
  b_seurat <- process_CM(data_seurat)
  b_seurat[is.na(b_seurat)] <-0
  
  seurat <- pheatmap(b_seurat, color=color_red_green, cluster_rows = F, cluster_cols=F, display_numbers=T, cellheight=25, cellwidth=25, na_col= '#cc5a4e', main=paste0('Seurat-LT rej', dataSet_name, '\n Mean Accuracy (per Cell Type): ', mean_acc(b_seurat)), legend=T, breaks = myBreaks, fontsize = fontsize, fontsize_number = fontsize_number)
  write.xlsx(b_seurat, file=file_xlsx, sheetName='seurat_raw', row.names = TRUE, append=TRUE)
  
  plot_dims <- get_plot_dims(jind)
  
  pdf(file.path(out_path, paste0(dataSet, '_CM.pdf')), family="Times", height = plot_dims$height *2  , width = plot_dims$width*2)
  grid.arrange(grobs = list(jind_raw[[4]], jind[[4]], seurat_raw[[4]], seurat[[4]]), ncol=2)
  dev.off()
  
  pdf(file.path(out_path, paste0(dataSet, '_CM_onlySJ.pdf')), family="Times", height = plot_dims$height*2, width = plot_dims$width*2 + 1)
  grid.arrange(grobs = list(jind[[4]], jind_raw[[4]], seurat_raw[[4]]), ncol=2, nrow =2)
  dev.off()
}


path = "/home/mohit/mohit/seq-rna/Comparison/datasets"
out_path = "/home/mohit/mohit/seq-rna/Comparison/JIND_DE/Plots/MohitPlots"

draw_cfmt("pancreas_01", path = path, out_path = out_path)
draw_cfmt("pancreas_02", path = path, out_path = out_path)
draw_cfmt("human_blood_01", path = path, out_path = out_path)
draw_cfmt("human_dataset_random", path = path, out_path = out_path)
draw_cfmt("mouse_atlas_random", path = path, out_path = out_path)

