
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
library(dendsort)


"graphContrast" <- function(data, name, Bth, FCth, namecol) {
	hist(data$P.Value, main = name, xlab = "p-value", ylab = "Genes");
	hist(data$logFC, main = name, xlab="Log2(FoldChange)", ylab="Genes", 50);
    #hist(treatm_vs_ctrl$B, main = name, xlab="B", ylab="Genes", 50);
	volcanoCol(data, Bth, FCth, name, namecol)
}

"volcanoCol" <- function(res, Bth, FCth, title, namecol) {
	colVP <- rep("black", length(res$B))
	colVP[which(res$B>Bth & res$logFC>FCth)] <- "red"
	colVP[which(res$B>Bth & res$logFC<(FCth*(-1)))] <- "green"
	plot(res$logFC, (-1)*log(res$P.Value), pch = ".", col = colVP, main = title, xlab = "foldchange", ylab = "-log(pvalue)")
	abline(v = FCth)
	abline(v = (FCth*(-1)))
	abline(h = (-1)*log(max(res$P.Value[res$B>Bth])))
	selFC <- res$logFC[which(res$B>Bth & abs(res$logFC)>FCth)]
	colVP2 <- rep("red", length(selFC))
	colVP2[which(selFC<((-1)*FCth))] <- "green"
	if (length(res[which(res$B>Bth & abs(res$logFC)>FCth),namecol])>0)
		text(res$logFC[which(res$B>Bth & abs(res$logFC)>FCth)], (-1)*log(res$P.Value[which(res$B>Bth & abs(res$logFC)>FCth)]), res[which(res$B>Bth & abs(res$logFC)>FCth),namecol], pos = 3, cex=0.7, offset = 0.5, col = colVP2)
}

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


callback_reverse = function(hc, mat)
{
dendsort(hc, isReverse=T)

}

pd <- import("pandas")

gencode22 <- read.table('/home/sevastopol/data/gserranos/UTILS/gencode.v22.annotation', header=T)
gencode22$gene_id <- str_extract(gencode22$gene_id, 'ENSG[R]?[0-9]+')
data_path = '/home/sevastopol/data/gserranos/JIND_DE/Data/'         # Use forward slash at the end of the path!
plots_path = '/home/sevastopol/data/gserranos/JIND_DE/Plots/Dec20/' # Use forward slash at the end of the path!

DATASETS <- list.files(data_path)


color_red_green <- colorRampPalette(c('#4E62CC','#D8DBE2' , '#BA463E'))(50)
color_red_green <- colorRampPalette(c('#cb5b4c','#D8DBE2', '#1aab2d'))(50)
myBreaks <- c(seq(0,  0.2, length.out= 20), seq(0.21, 0.79, length.out=10), seq(0.8, 1, length.out=20))
for (dataSet in DATASETS){
	print(dataSet)
	switch(dataSet,
		HumanDatasetRandom = {
			dataSet_name = 'Human Hematopoiesis'},
		Mouse_atlas = {
			dataSet_name = 'Mouse Atlas'},
		Pancreas_01 = {
			dataSet_name = 'Pancreas Bar16-Mur16'},
		Pancreas_02 = {
			dataSet_name = 'Pancreas Bar16-Seg16'},
		PBMC = {
			dataSet_name = 'PBMC 10x_v3-10x_v5'}
	)

	annotation <- pd$read_pickle(paste0(data_path,dataSet,'/JIND_assignmentbrftune.pkl'))
	annotation$cell_names <- rownames(annotation)

	annotation_seurat <- pd$read_pickle(paste0(data_path,dataSet,'/seurat_assignment.pkl'))
	annotation_seurat$cell_names <- rownames(annotation_seurat)

	data <- as.data.frame.matrix(table(annotation$raw_predictions, annotation$labels))
	b <- process_CM(data)
	b[is.na(b)] <- 0

	jind_raw <- pheatmap(b, color=color_red_green, cluster_rows = F, cluster_cols=F, display_numbers=T, cellheight=25, cellwidth=25, na_col= '#cc5a4e', main=paste0('JIND raw ', dataSet_name, ' mean accuracy: ', mean_acc(b)), legend=F, breaks = myBreaks, fontsize = 7)
	file_xlsx <- paste0(plots_path,dataSet,'_CM.xlsx' )
	if (file.exists(file_xlsx)) {
  		file.remove(file_xlsx)
	}
	b[is.na(b)] <-0
	write.xlsx(b, file=file_xlsx, sheetName='JIND_raw', row.names = TRUE, append=TRUE)

	data <- as.data.frame.matrix(table(annotation$predictions, annotation$labels))
	b <- process_CM(data)
	b[is.na(b)] <- 0

	jind <- pheatmap(b, color=color_red_green, cluster_rows = F, cluster_cols=F, display_numbers=T, cellheight=25, cellwidth=25, na_col= '#cc5a4e', main=paste0('JIND ', dataSet_name, ' mean accuracy: ', mean_acc(b)), legend=F, breaks = myBreaks, fontsize = 7)
	b[is.na(b)] <-0
	write.xlsx(b, file=file_xlsx, sheetName='JIND', row.names = TRUE, append=TRUE)

	data_seurat <- as.data.frame.matrix(table(annotation_seurat$raw_predictions, annotation_seurat$labels))
	b_seurat <- process_CM(data_seurat)
	b_seurat[is.na(b_seurat)] <-0

	seurat_raw <- pheatmap(b_seurat, color=color_red_green, cluster_rows = F, cluster_cols=F, display_numbers=T, cellheight=25, cellwidth=25, na_col= '#cc5a4e', main=paste0('SEURAT ', dataSet_name, ' mean accuracy: ', mean_acc(b_seurat)), legend=T, breaks = myBreaks, fontsize = 7)
	write.xlsx(b_seurat, file=file_xlsx, sheetName='seurat', row.names = TRUE, append=TRUE)

	data_seurat <- as.data.frame.matrix(table(annotation_seurat$predictions, annotation_seurat$labels))
	b_seurat <- process_CM(data_seurat)
	b_seurat[is.na(b_seurat)] <-0

	seurat <- pheatmap(b_seurat, color=color_red_green, cluster_rows = F, cluster_cols=F, display_numbers=T, cellheight=25, cellwidth=25, na_col= '#cc5a4e', main=paste0('SEURAT rej', dataSet_name, ' mean accuracy: ', mean_acc(b_seurat)), legend=T, breaks = myBreaks, fontsize = 7)
	write.xlsx(b_seurat, file=file_xlsx, sheetName='seurat_raw', row.names = TRUE, append=TRUE)

	plot_dims <- get_plot_dims(jind)

	pdf(paste0(plots_path,dataSet,'_CM.pdf'), height = plot_dims$height *2  , width = plot_dims$width*2)
	grid.arrange(grobs = list(jind_raw[[4]], jind[[4]], seurat_raw[[4]], seurat[[4]]), ncol=2)
	dev.off()
	
	pdf(paste0(plots_path,dataSet,'_CM_onlySJ.pdf'), height = plot_dims$height *3, width = plot_dims$width + 5)
	grid.arrange(grobs = list(jind[[4]], jind_raw[[4]], seurat_raw[[4]]), ncol=1, nrow =3)
	dev.off()

}


DE_with_TSNE_train <- function(dataSet, target, obj, genes_displ, plot_selected_genes = NULL,data_path = data_path, plots_path = plots_path){

	df  <- pd$read_pickle(paste0(data_path,dataSet,'/train.pkl'))
	annotation <- data.frame(cell_names = rownames(df), labels = df$labels, stringsAsFactors = F)
	rownames(annotation) <- annotation$cell_names

	all_data <- t(df[, -which(colnames(df) %in% c('labels'))])

	stand_dev <- apply(all_data,1, sd)
	all_data <- all_data[!rownames(all_data) %in% names(which(stand_dev ==0)), ]

	# tsne_out <- Rtsne(as.matrix(t(all_data)))
	# plotter <-  as.data.frame(tsne_out$Y)

	annotation$labels <- gsub(' ', '_', annotation$labels)

	G1 <- annotation[annotation$labels == target , 'cell_names']
	G2 <- annotation[ annotation$labels == obj, 'cell_names']

	selection <- c(G1, G2)

	data_tmp <- all_data[ , colnames(all_data) %in%  selection]
	stand_dev <- apply(data_tmp,1, sd)
	data_tmp <- data_tmp[!rownames(data_tmp) %in% names(which(stand_dev ==0)), ]


	DESIGN <- data.frame(cells = colnames(data_tmp))
	labels <- c()
	for (cell in DESIGN$cells){
		ifelse(cell %in% G1, labels <- c(labels, 'G1'), labels <- c(labels, 'G2'))
	}
	DESIGN$labels <- labels

	design_tmp <- as.matrix(DESIGN[, 'labels'])
	design_tmp[design_tmp != 'G1'] <-1
	design_tmp[design_tmp == 'G1'] <-0
	design_tmp <- model.matrix(~0+as.factor(design_tmp[,1]))
	colnames(design_tmp) <- c('G1', 'G2')
	rownames(design_tmp) <- colnames(data_tmp)


	### DE
	# create the linear model
	fit_tmp <- lmFit(data_tmp, design_tmp)
	# model correction
	fit_tmp <- eBayes(fit_tmp)
	# results <- topTable(fit_tmp, n=Inf)
	x <- paste0('G1', '-', 'G2')
	contrast_mat_tmp <- makeContrasts(contrasts=x, levels= c('G1', 'G2'))
	fit2_tmp <- contrasts.fit(fit_tmp, contrast_mat_tmp)
	fit2_tmp <- eBayes(fit2_tmp)
	tmp   <- topTable(fit2_tmp, adjust="BH", n=Inf)
	tmp$gene_name <- rownames(tmp)
	# tmp <- merge(tmp, gencode22, by= 'gene_name')



	tmp[tmp$P.Value == 0, 'P.Value'] <- 1.445749e-281
	tmp[tmp$adj.P.Val == 0, 'P.Value'] <- 1.445749e-281

	# pdf('./Plots/12_CD14.Mono.2_VP.pdf')
	# 	graphContrast(tmp," (Ctrl B > 0, FC>0.5)", 0, 0.5, 1)
	# dev.off()

	tmp$logpval <- -log(tmp$P.Value)
	tmp2 <- tmp[ order(-tmp$logpval), c('gene_name', 'logFC', 'P.Value' ,'logpval')]


	data2heat <- data_tmp[rownames(data_tmp) %in%  tmp[tmp$adj.P.Val<0.001,'gene_name'],]
	# data2heat <- data_tmp[rownames(data_tmp) %in%  tmp[tmp$adj.P.Val<0.05 & abs(tmp$logFC)>0.5,'gene_name'],]
	data2heat[data2heat > 5] <- 5


	colors = c(rgb(31, 119, 180, max = 255), rgb(255, 127, 14, max =255),
				rgb(44, 160, 44, max = 255), rgb(214, 39, 40, max =255),
				rgb(148, 103, 189, max = 255), rgb(140, 86, 75, max =255),
				rgb(227, 119, 194, max = 255), rgb(127, 127, 127, max =255),
				rgb(188, 189, 34, max = 255), rgb(23, 190, 207, max =255))


	ann <- data.frame(Group = annotation[rownames(annotation) %in% colnames(data2heat), 'labels'])
	ann$Group <-  ifelse(ann$Group == target, 'G1', 'G2')
	rownames(ann) <- colnames(data2heat)
	Group <- colors[1:length(unique(annotation$labels))]
	names(Group) <-names(sort(table(as.character(annotation$labels) ), decreasing=T))
	# names(Group)  <- c(unique(ann$Group))
	anno_colors <- c(Group[target],Group[obj])
	names(anno_colors) <- c('G1', 'G2')
	# anno_colors   <- list(anno_colors = anno_colors)
	anno_colors <- list(Group = anno_colors)
	data2heat2 <- rbind(data2heat[ c('CD14' , 'FCGR3A'),], data2heat[ 1:23,])
	HM <- pheatmap( data2heat2, cluster_rows = F, treeheight_row = 0, annotation_col = ann,annotation_colors = anno_colors, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", cellheight= 4,cellwidth= 0.6, show_colnames = F, main = paste0('Heatmap between ',target,' classified as ',target,' (G1)\n and ',target,' classified as ',obj,' (G2)'), fontsize = 8,fontsize_row=4)
	plot_dims <- get_plot_dims(HM)
	pdf('Test.pdf',height = plot_dims$height, width = plot_dims$width)
	print(HM)
	dev.off()



}


DE_with_TSNE <- function(dataSet, target, obj, genes_displ, plot_selected_genes = NULL,data_path = data_path, plots_path = plots_path){

	df  <- pd$read_pickle(paste0(data_path,dataSet,'/test.pkl'))
	annotation <- pd$read_pickle(paste0(data_path, dataSet,'/JIND_assignmentbrftune.pkl'))
	annotation$cell_names <- rownames(annotation)

	all_data <- t(df[, -which(colnames(df) %in% c('labels'))])
	var_5k=FALSE
	if (var_5k==TRUE & !dataSet %in% c('PBMC', 'Pancreas_01', 'Pancreas_02')){
		var_ <- apply(all_data, 1, var)
		names(var_) <- rownames(all_data)
		genes_2_keep <- names(sort(var_, decreasing=T)[1:5000])
		all_data <- all_data[rownames(all_data) %in% genes_2_keep, ]
	}

	stand_dev <- apply(all_data,1, sd)
	all_data <- all_data[!rownames(all_data) %in% names(which(stand_dev ==0)), ]

	# tsne_out <- Rtsne(as.matrix(t(all_data)))
	# plotter <-  as.data.frame(tsne_out$Y)

	annotation$labels <- gsub(' ', '_', annotation$labels)
	annotation$predictions <- gsub(' ', '_', annotation$predictions)
	annotation$raw_predictions <- gsub(' ', '_', annotation$raw_predictions)

	G1 <- annotation[annotation$labels == target & annotation$prediction == target, 'cell_names']
	G2 <- annotation[annotation$labels == target & annotation$prediction == obj, 'cell_names']

	selection <- c(G1, G2)

	data_tmp <- all_data[ , colnames(all_data) %in%  selection]
	stand_dev <- apply(data_tmp,1, sd)
	data_tmp <- data_tmp[!rownames(data_tmp) %in% names(which(stand_dev ==0)), ]

	DESIGN <- data.frame(cells = colnames(data_tmp))
	labels <- c()
	for (cell in DESIGN$cells){
		ifelse(cell %in% G1, labels <- c(labels, 'G1'), labels <- c(labels, 'G2'))
	}
	DESIGN$labels <- labels

	design_tmp <- as.matrix(DESIGN[, 'labels'])
	design_tmp[design_tmp != 'G1'] <-1
	design_tmp[design_tmp == 'G1'] <-0
	design_tmp <- model.matrix(~0+as.factor(design_tmp[,1]))
	colnames(design_tmp) <- c('G1', 'G2')
	rownames(design_tmp) <- colnames(data_tmp)


	### DE
	# create the linear model
	fit_tmp <- lmFit(data_tmp, design_tmp)
	# model correction
	fit_tmp <- eBayes(fit_tmp)
	# results <- topTable(fit_tmp, n=Inf)
	x <- paste0('G1', '-', 'G2')
	contrast_mat_tmp <- makeContrasts(contrasts=x, levels= c('G1', 'G2'))
	fit2_tmp <- contrasts.fit(fit_tmp, contrast_mat_tmp)
	fit2_tmp <- eBayes(fit2_tmp)
	tmp   <- topTable(fit2_tmp, adjust="fdr", n=Inf)
	tmp$gene_name <- rownames(tmp)
	# tmp <- merge(tmp, gencode22, by= 'gene_name')


	tmp[tmp$P.Value == 0, 'P.Value'] <- 1.445749e-281
	tmp[tmp$adj.P.Val == 0, 'P.Value'] <- 1.445749e-281

	file_xlsx <- paste0(plots_path,dataSet,'_',target,'Vs',obj,'_DE.xlsx')
	if (file.exists(file_xlsx)) {
  		file.remove(file_xlsx)
	}
	write.xlsx( tmp[tmp$adj.P.Val<0.05,] , file=file_xlsx, sheetName='JIND+', row.names = TRUE, append=TRUE)
	# pdf('./Plots/12_CD14.Mono.2_VP.pdf')
	# 	graphContrast(tmp," (Ctrl B > 0, FC>0.5)", 0, 0.5, 1)
	# dev.off()

	tmp$logpval <- -log(tmp$P.Value)
	tmp2 <- tmp[ order(-tmp$logpval), c('gene_name', 'logFC', 'P.Value' ,'logpval', 'adj.P.Val')]


	data2heat <- data_tmp[rownames(data_tmp) %in%  tmp[tmp$adj.P.Val<0.05,'gene_name'],]
	# data2heat <- data_tmp[rownames(data_tmp) %in%  tmp[tmp$adj.P.Val<0.05 & abs(tmp$logFC)>0.5,'gene_name'],]
	data2heat[data2heat > 5] <- 5
	print(nrow(data2heat))
	colors = c(rgb(31, 119, 180, max = 255), rgb(255, 127, 14, max =255),
				rgb(44, 160, 44, max = 255), rgb(214, 39, 40, max =255),
				rgb(148, 103, 189, max = 255), rgb(140, 86, 75, max =255),
				rgb(227, 119, 194, max = 255), rgb(127, 127, 127, max =255),
				rgb(188, 189, 34, max = 255), rgb(23, 190, 207, max =255))


	ann <- data.frame(Group = annotation[rownames(annotation) %in% colnames(data2heat), 'predictions'])
	ann$Group <-  ifelse(ann$Group == target, 'G1', 'G2')
	rownames(ann) <- colnames(data2heat)
	Group <- colors[1:length(unique(annotation$labels))]
	names(Group) <-names(sort(table(as.character(annotation$labels) ), decreasing=T))
	# names(Group)  <- c(unique(ann$Group))
	anno_colors <- c(Group[target],Group[obj])
	names(anno_colors) <- c('G1', 'G2')
	# anno_colors   <- list(anno_colors = anno_colors)


	HM <- pheatmap( data2heat, cluster_rows = T, treeheight_row = 0, annotation_col = ann, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", cellheight= 4,cellwidth= 2, show_colnames = F, main = paste0('Heatmap between ',target,' classified as ',target,' (G1)\n and ',target,' classified as ',obj,' (G2)'), fontsize = 8,fontsize_row=4, callback = callback_reverse)

	plot_dims <- get_plot_dims(HM)


	
	cluster_names <- HM$tree_col$labels
	cluster_pos <- HM$tree_col$order

	cluster_seq <- which(cluster_names[cluster_pos] %in% G2)
	cluster_list <- split(cluster_seq, cumsum(c(1, diff(cluster_seq) != 1)))

	counter <- 1
	shape <- data.frame(cell_id = colnames(data2heat), Shape = 'G1', stringsAsFactors=FALSE)
	for (cl in cluster_list){
		shape[shape$cell_id %in% cluster_names[cluster_pos][cl] , 'Shape'] <- paste0('G2_', counter)
		counter <- counter +1
	}

	ann <- data.frame(cell_id = colnames(data2heat), Group = annotation[rownames(annotation) %in% colnames(data2heat), 'predictions'])
	ann$Group <-  ifelse(ann$Group == target, 'G1', 'G2')
	ann <- merge(ann, shape, by='cell_id')
	ann$cell_id <- NULL
	names(ann) <- c('Group', 'Cluster')
	rownames(ann) <- colnames(data2heat)
	# Group_         <- c(Group[target], Group[obj])
	
	Cluster         <- colors[which(!colors %in% c(anno_colors))][1:counter]
	# names(Group)  <- c(levels(ann$Group))
	# names(Group_) <- c(unique(ann$Group))
	names(Cluster) <- c(unique(ann$Cluster))
	# Cluster <- 
	Cluster <- setNames(c(Cluster, rgb(127, 127, 127, alpha = 50, max=255)), c(names(Cluster), 'G0'))


	anno_colors2  <- list(Group = anno_colors, Cluster=Cluster)
	anno_colors2$Cluster['G1'] <- anno_colors['G1']
	HM <- pheatmap( data2heat, cluster_rows = T, treeheight_row = 0, annotation_col = ann, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", cellheight= 4,cellwidth= 2, annotation_colors = anno_colors2, show_colnames = F, main = paste0('Heatmap between ',target,' classified as ',target,' (G1)\n and ',target,' classified as ',obj,' (G2)'), fontsize = 8,fontsize_row=4, callback = callback_reverse)

	plot_dims <- get_plot_dims(HM)

	pdf(paste0(plots_path,dataSet,'_',target,'Vs',obj,'_HM_All.pdf'), height = plot_dims$height, width = plot_dims$width)
	print(HM)
	dev.off()

	if(dataSet %in% c('Pancreas_01', 'PBMC')){
		switch(dataSet,
		Pancreas_01 = { golden_boys <- c('KRT19', 'PDX1', 'SOX9', 'UEA1', 'GP2', 'CD142', 'PRSS1', 'CTRC', 'CPA1', 'AMY2A', 'SYCN', 'RBPJL', 'MIST1', 'HNF1B', 'PTF1A', 'CA19.9', 'PARM1', 'GP2', 'CD142', 'RBPJ', 'MYC')},
		PBMC        = { golden_boys <- c('CD14', 'FCGR3A')}
		)
		golden_present <- tmp2[tmp2$gene_name %in% golden_boys,'gene_name']
		
		krtdata <- as.data.frame(data2heat[rownames(data2heat) %in% golden_present,])
		# rownames(krtdata) <- 'KRT19'
		data2heat2 <- data2heat[!rownames(data2heat) %in% golden_present,]
		data2heat2 <- rbind(krtdata, data2heat2)
		ann$Cluster <- NULL

		HM_sel <- pheatmap( data2heat2[1:genes_displ,],cluster_cols = HM$tree_col, cluster_rows = T, treeheight_row = 0, annotation_col = ann, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",cellwidth= 2, annotation_colors = anno_colors2, show_colnames = F, cellheight= 7, main = paste0('Heatmap between ',target,' classified as ',target,' (G1)\n and ',target,' classified as ',obj,' (G2)'), fontsize = 8,fontsize_row=8)
		plot_dims <- get_plot_dims(HM_sel)
		pdf(paste0(plots_path,dataSet,'_',target,'Vs',obj,'_HM_Sel.pdf'), height = plot_dims$height, width = plot_dims$width)
			print(HM_sel)
		dev.off()
	}

	##
	# Plot only the marker genes
	##
	if (!is.null(plot_selected_genes)){
		data2heat2 <- data_tmp[rownames(data_tmp) %in% plot_selected_genes,]
		ann <- data.frame(cell_id = colnames(data2heat2), Group = annotation[rownames(annotation) %in% colnames(data2heat2), 'predictions'])
		ann$Group <-  ifelse(ann$Group == target, 'G1', 'G2')
		ann$cell_id <- NULL
		rownames(ann) <- colnames(data2heat2)
		data2heat2 <- data2heat2[,order(-data2heat2[plot_selected_genes[1],])]
		HM_1 <- pheatmap( data2heat2[plot_selected_genes[1],,drop=F], cluster_rows = F, cluster_cols = F, treeheight_row = 0, annotation_col = ann, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", cellheight= 4,cellwidth= 2, annotation_colors = anno_colors, show_colnames = F, main = paste0('Heatmap between ',target,' classified as ',target,' (G1)\n and ',target,' classified as ',obj,' (G2)'), fontsize = 8,fontsize_row=4)
		data2heat2 <- data2heat2[,order(-data2heat2[plot_selected_genes[2],])]
		HM_2 <- pheatmap( data2heat2[plot_selected_genes[2],,drop=F], cluster_rows = F, cluster_cols = F, treeheight_row = 0, annotation_col = ann, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", cellheight= 4,cellwidth= 2, annotation_colors = anno_colors, show_colnames = F, main = paste0('Heatmap between ',target,' classified as ',target,' (G1)\n and ',target,' classified as ',obj,' (G2)'), fontsize = 8,fontsize_row=4)

		plot_dims <- get_plot_dims(HM_1)

		pdf(paste0(plots_path,dataSet,'_',target,'Vs',obj,'_HM_2genes.pdf'), height = plot_dims$height*3, width = plot_dims$width)
		grid.arrange(grobs = list(HM_1[[4]], HM_2[[4]]), nrow=2)
		dev.off()
	}
	#############
	# TSNE
	#############

	G3 <- annotation[annotation$labels == target & annotation$prediction == target, 'cell_names']
	G4 <- annotation[annotation$labels == obj & annotation$prediction == obj, 'cell_names']

	all_data2 <- all_data[, colnames(all_data) %in%  c(G3, G4, G2)]
	# all_data2 <- all_data[, colnames(all_data) %in% annotation[annotation$labels %in% c(target, obj), 'cell_names']]
	# set.seed(123)
	tsne_out <- Rtsne(as.matrix(t(all_data2)))
	plotter <-  as.data.frame(tsne_out$Y)

	plotter$labels <- annotation[rownames(annotation) %in% colnames(all_data2), 'labels']
	plotter$predictions <- annotation[rownames(annotation) %in% colnames(all_data2), 'predictions']
	plotter$cell_id <- rownames(annotation[rownames(annotation) %in% colnames(all_data2),])

	labels <- data.frame(cell_id = colnames(all_data2) )
	labels$Cluster <- 'G0'
	labels[labels$cell_id %in% G1,'Cluster'] <- 'G1'
	counter <- 1
	for (cl in cluster_list){
		labels[labels$cell_id %in% cluster_names[cluster_pos][cl], 'Cluster'] <- paste0('G2_', counter)
		counter <- counter +1
	}
	
	plotter <- merge(plotter, labels, by ='cell_id')


	euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

	ratio_dist <- function(a, b, c) {
		return(as.character(euc.dist(unlist(a),unlist(b))/euc.dist(unlist(a),unlist(c))))
	}

	mean_dist <- list()
	for (gr in unique(plotter$Cluster)){
		mean_dist[[gr]] <- c(c(mean(plotter[plotter$Cluster %in% c(gr), 'V1']), mean(plotter[plotter$Cluster %in% c(gr), 'V2'])))
	}
    # print(mean_dist)
	title <- c()
	for (i in seq_along(mean_dist)){
		print(i)
		tmp <- paste0(names(mean_dist[i]), ' ratio dist.  ', paste(ratio_dist(mean_dist[i],  mean_dist['G0'], mean_dist['G1'])))
		title <- paste(title, tmp, sep=' \n')
	}
	title <- paste(title, paste0('Ratio:  ',obj,'/',target,''), sep='\n')

	pdf(paste0(plots_path,dataSet,'_',target,'Vs',obj,'_TSNE.pdf'), 15,15)
	print(ggplot(plotter, aes(x=V1, y=V2, color = Cluster, label = labels)) + theme_classic() + geom_text()  +scale_color_manual(values = anno_colors2$Cluster) + ggtitle(title)+ theme(legend.text=element_text(size=15)))
	print(ggplot(plotter, aes(x=V1, y=V2, color = Cluster)) + theme_classic() + geom_point(alpha = 0.8,size = 8)  +scale_color_manual(values = anno_colors2$Cluster) + ggtitle(title)  + 
	theme(legend.text=element_text(size=15)) +
	theme(axis.title.x=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank()) )
	dev.off()

	probs <- data.frame(cell_id = NULL, prob = NULL)
	for (cell_name in plotter$cell_id){
		tmp <- data.frame(cell_id = cell_name, prob = annotation[annotation$cell_names == cell_name, as.character(annotation[annotation$cell_names == cell_name, 'labels'])])
		probs <- rbind(probs,tmp)
	}

	plotter <- merge(plotter, probs, by='cell_id')
	library(wesanderson)
	pal <- wes_palette("Zissou1", 100, type = "continuous")

	pdf(paste0(plots_path,dataSet,'_',target,'Vs',obj,'_TSNE_probs.pdf'), 15,15)
	print(ggplot(plotter, aes(x=V1, y=V2, color = prob, label = labels)) + theme_bw() + geom_text()  +   scale_color_gradientn(colours = pal))
	print(ggplot(plotter, aes(x=V1, y=V2, fill = prob)) + theme_bw() + geom_point(aes(fill=prob),colour="black", shape=21, size = 5) +scale_fill_gradientn(colours = pal) )# scale_fill_viridis(option="inferno"))
	dev.off()
	return(setNames(list(cluster_list, cluster_names), c('cluster_list', 'cluster_names')))
}



TSNE_for_all <- function(dataSet){

	df  <- pd$read_pickle(paste0('/home/sevastopol/data/gserranos/JIND_DE/Data/',dataSet,'/test.pkl'))

	annotation <- pd$read_pickle(paste0('/home/sevastopol/data/gserranos/JIND_DE/Data/',dataSet,'/JIND_assignmentbrftune.pkl'))
	annotation$cell_names <- rownames(annotation)

	all_data <- t(df[, -which(colnames(df) %in% c('labels'))])
	tsne_out <- Rtsne(as.matrix(t(all_data)))
	plotter <-  as.data.frame(tsne_out$Y)

	plotter$labels <- annotation[rownames(annotation) %in% colnames(all_data), 'labels']
	plotter$predictions <- annotation[rownames(annotation) %in% colnames(all_data), 'predictions']
	plotter$cell_id <- rownames(annotation[rownames(annotation) %in% colnames(all_data),])


	probs <- data.frame(cell_id = NULL, prob = NULL)
	for (cell_name in plotter$cell_id){
		tmp <- data.frame(cell_id = cell_name, prob = annotation[annotation$cell_names == cell_name, as.character(annotation[annotation$cell_names == cell_name, 'labels'])])
		probs <- rbind(probs,tmp)
	}

	plotter <- merge(plotter, probs, by='cell_id')
	library(wesanderson)
	pal <- wes_palette("Zissou1", 100, type = "continuous")

	colors_labels = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
	colors_labels_ = sample(colors_labels, length(unique(plotter$labels)))
	pdf(paste0(plots_path,dataSet,'_TSNE_probs.pdf'), 15,15)
	print(ggplot(plotter, aes(x=V1, y=V2, fill = labels)) + theme_bw() + geom_point(aes(fill=labels) , colour="black", shape=21, size = 5, alpha = 0.7) + scale_fill_manual(values=colors_labels_) )
	
	print(ggplot(plotter, aes(x=V1, y=V2, color = prob, label = labels)) + theme_bw() + geom_text(alpha = 0.6)  +   scale_color_gradientn(colours = pal))
	
	print(ggplot(plotter, aes(x=V1, y=V2, color = prob, label = predictions)) + theme_bw() + geom_text(alpha = 0.6)  +   scale_color_gradientn(colours = pal))
	
	print(ggplot(plotter, aes(x=V1, y=V2, fill = prob)) + theme_bw() + geom_point(aes(fill=prob),colour="black", shape=21, size = 5, alpha = 0.6) +scale_fill_gradientn(colours = pal) )
	dev.off()


}


NegativeControl_DE <- function(dataSet, target, plots_path){

	df  <- pd$read_pickle(paste0(data_path,dataSet,'/test.pkl'))
	annotation <- pd$read_pickle(paste0(data_path, dataSet,'/JIND_assignmentbrftune.pkl'))
	annotation$cell_names <- rownames(annotation)

	all_data <- t(df[, -which(colnames(df) %in% c('labels'))])

	G1 <- annotation[annotation$labels == target & annotation$prediction == target, 'cell_names']
	G2 <- sample(G1, round(length(G1)/2))
	G1 <- setdiff(G1, G2)

	selection <- c(G1, G2)

	data_tmp <- all_data[ , colnames(all_data) %in%  selection]

	DESIGN <- data.frame(cells = colnames(data_tmp))
	labels <- c()
	for (cell in DESIGN$cells){
		ifelse(cell %in% G1, labels <- c(labels, 'G1'), labels <- c(labels, 'G2'))
	}
	DESIGN$labels <- labels

	design_tmp <- as.matrix(DESIGN[, 'labels'])
	design_tmp[design_tmp != 'G1'] <-1
	design_tmp[design_tmp == 'G1'] <-0
	design_tmp <- model.matrix(~0+as.factor(design_tmp[,1]))
	colnames(design_tmp) <- c('G1', 'G2')
	rownames(design_tmp) <- colnames(data_tmp)


	### DE
	# create the linear model
	fit_tmp <- lmFit(data_tmp, design_tmp)
	# model correction
	fit_tmp <- eBayes(fit_tmp)
	# results <- topTable(fit_tmp, n=Inf)
	x <- paste0('G1', '-', 'G2')
	contrast_mat_tmp <- makeContrasts(contrasts=x, levels= c('G1', 'G2'))
	fit2_tmp <- contrasts.fit(fit_tmp, contrast_mat_tmp)
	fit2_tmp <- eBayes(fit2_tmp)
	tmp   <- topTable(fit2_tmp, adjust="fdr", n=Inf)
	tmp$gene_name <- rownames(tmp)



	tmp[tmp$adj.P.Val == 0, 'P.Value'] <- 1.445749e-281

	# pdf('./Plots/Final_PBMC_NC_classMonCD_VP.pdf')
	# 	graphContrast(tmp," (Ctrl B > 0, FC>0.5)", 0, 0.5, 1)
	# dev.off()

	tmp$logpval <- -log(tmp$P.Value)
	tmp2 <- tmp[ order(-tmp$logpval), c('gene_name', 'logFC', 'P.Value' ,'logpval', 'adj.P.Val')]
	tmp2 <- tmp[tmp$P.Value < 0.05,]

	file_xlsx <- paste0(plots_path,dataSet,'_',target,'_NC.xlsx')
	if (file.exists(file_xlsx)) {
			file.remove(file_xlsx)
	}

	write.xlsx(tmp2 , file=file_xlsx, sheetName='P.Value Lt 0.05', row.names = TRUE, append=TRUE)

	data2heat <- data_tmp[rownames(data_tmp) %in%  tmp2$gene_name,]
	data2heat[data2heat > 5] <- 5

	colors = c(rgb(31, 119, 180, max = 255), rgb(255, 127, 14, max =255),
					rgb(44, 160, 44, max = 255), rgb(214, 39, 40, max =255),
					rgb(148, 103, 189, max = 255), rgb(140, 86, 75, max =255),
					rgb(227, 119, 194, max = 255), rgb(127, 127, 127, max =255),
					rgb(188, 189, 34, max = 255), rgb(23, 190, 207, max =255))
	ann <- data.frame(Var1 = annotation[rownames(annotation) %in% colnames(data2heat), 'predictions'])
	rownames(ann) <- colnames(data2heat)
	Var1        <- c(colors[1], colors[2])
	names(Var1) <- c(levels(ann$Var1))
	anno_colors <- list(Var1 = Var1)


	ann <- data.frame(Group = c(G1, G2))
	ann$Group <-  ifelse(ann$Group %in% G1, 'G1', 'G2')
	rownames(ann) <- colnames(data2heat)
	Group         <- c(Var1)
	names(Group) <- c('G1', 'G2')
	anno_colors   <- list(Group = Group)

	p <- pheatmap( data2heat, cluster_rows = T, treeheight_row = 0, annotation_col = ann, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", cellheight= 10,annotation_colors = anno_colors, show_colnames = F, main = paste0('Heatmap between ',target,' classified as ',target,' (G1)\n  and ',target,' (G2) \nNEGATIVE CONTROL'), fontsize = 8,fontsize_row=10)

	plot_dims <- get_plot_dims(p)

	pdf(paste0(plots_path,dataSet,'_',target,'_NC_HM.pdf'), height = plot_dims$height, width = plot_dims$width)
	print(p)
	dev.off()

}


# dataSet <- 'PBMC'
# target <- 'Monocyte_FCGR3A'
# obj    <- 'Monocyte_CD14'

# dataSet <- 'Pancreas_01'
# target <- 'ductal'
# obj    <- 'acinar'

a = DE_with_TSNE('Pancreas_01', 'ductal', 'acinar', 25, ,data_path = data_path, plots_path = plots_path)
a = DE_with_TSNE('PBMC', 'Monocyte_FCGR3A', 'Monocyte_CD14', 25,data_path = data_path, plots_path = plots_path)
DE_with_TSNE("HumanDatasetRandom", )


TSNE_for_all('PBMC')
TSNE_for_all('Pancreas_01')

NegativeControl_DE('Pancreas_01', 'ductal', plots_path = plots_path)
NegativeControl_DE('PBMC', 'Monocyte_FCGR3A', plots_path = plots_path)


a$cluster_names[a$cluster_list[[3]]]
dataSet <- 'Pancreas_01'
target <- 'acinar'
obj <- 'custom'
df  <- pd$read_pickle(paste0('/home/sevastopol/data/gserranos/JIND_DE/Data/',dataSet,'/test.pkl'))
	annotation <- pd$read_pickle(paste0('/home/sevastopol/data/gserranos/JIND_DE/Data/',dataSet,'/JIND_assignmentbrftune.pkl'))
	annotation$cell_names <- rownames(annotation)

	all_data <- t(df[, -which(colnames(df) %in% c('labels'))])
	# all_data2 <- t(df2[, -which(colnames(df2) %in% c('labels'))])
	var_5k=FALSE
	if (var_5k==TRUE & !dataSet %in% c('PBMC', 'Pancreas_01', 'Pancreas_02')){
		var_ <- apply(all_data, 1, var)
		names(var_) <- rownames(all_data)
		genes_2_keep <- names(sort(var_, decreasing=T)[1:5000])
		all_data <- all_data[rownames(all_data) %in% genes_2_keep, ]
	}

	stand_dev <- apply(all_data,1, sd)
	all_data <- all_data[!rownames(all_data) %in% names(which(stand_dev ==0)), ]

	# tsne_out <- Rtsne(as.matrix(t(all_data)))
	# plotter <-  as.data.frame(tsne_out$Y)

	annotation$labels <- gsub(' ', '_', annotation$labels)
	annotation$predictions <- gsub(' ', '_', annotation$predictions)
	annotation$raw_predictions <- gsub(' ', '_', annotation$raw_predictions)

	G1 <- annotation[annotation$labels == target & annotation$prediction == target, 'cell_names']
	G2 <- a$cluster_names[a$cluster_list[[3]]]

	selection <- c(G1, G2)

	data_tmp <- all_data[ , colnames(all_data) %in%  selection]
	stand_dev <- apply(data_tmp,1, sd)
	data_tmp <- data_tmp[!rownames(data_tmp) %in% names(which(stand_dev ==0)), ]

	DESIGN <- data.frame(cells = colnames(data_tmp))
	labels <- c()
	for (cell in DESIGN$cells){
		ifelse(cell %in% G1, labels <- c(labels, 'G1'), labels <- c(labels, 'G2'))
	}
	DESIGN$labels <- labels

	design_tmp <- as.matrix(DESIGN[, 'labels'])
	design_tmp[design_tmp != 'G1'] <-1
	design_tmp[design_tmp == 'G1'] <-0
	design_tmp <- model.matrix(~0+as.factor(design_tmp[,1]))
	colnames(design_tmp) <- c('G1', 'G2')
	rownames(design_tmp) <- colnames(data_tmp)


	### DE
	# create the linear model
	fit_tmp <- lmFit(data_tmp, design_tmp)
	# model correction
	fit_tmp <- eBayes(fit_tmp)
	# results <- topTable(fit_tmp, n=Inf)
	x <- paste0('G1', '-', 'G2')
	contrast_mat_tmp <- makeContrasts(contrasts=x, levels= c('G1', 'G2'))
	fit2_tmp <- contrasts.fit(fit_tmp, contrast_mat_tmp)
	fit2_tmp <- eBayes(fit2_tmp)
	tmp   <- topTable(fit2_tmp, adjust="BH", n=Inf)
	tmp$gene_name <- rownames(tmp)
	tmp <- merge(tmp, gencode22, by= 'gene_name')



	tmp[tmp$P.Value == 0, 'P.Value'] <- 1.445749e-281
	tmp[tmp$adj.P.Val == 0, 'P.Value'] <- 1.445749e-281

	# pdf('./Plots/12_CD14.Mono.2_VP.pdf')
	# 	graphContrast(tmp," (Ctrl B > 0, FC>0.5)", 0, 0.5, 1)
	# dev.off()

	tmp$logpval <- -log(tmp$P.Value)
	tmp2 <- tmp[ order(-tmp$logpval), c('gene_name', 'logFC', 'P.Value' ,'logpval')]


	data2heat <- data_tmp[rownames(data_tmp) %in%  tmp[tmp$adj.P.Val<0.05,'gene_name'],]

	colors = c(
	"G2"  = rgb(105, 191, 109, max=255),
	"G2_3"  = rgb(140, 86, 75, max=255), #BROWN
	"G2_4" = rgb(255, 127, 14, max=255), #ORANGE
	"G1" =  rgb(214, 38, 40, max=255), #RED
	"G2_5" =rgb(44, 160, 44, max=255), #GREEN
	"G2_2" = rgb(31, 119, 180, max=255), #BLUE
	"18_Plasma" = rgb(227, 119, 194, max=255), # PINK
	"G2" = rgb(148, 103, 189, max=255), #PURPLE
	"G0" = rgb(127, 127, 127, max=255), #GREY
	"G2_1"  = rgb(188, 189, 34, max=255)) # LIME

	ann <- data.frame(Group = annotation[rownames(annotation) %in% colnames(data2heat), 'predictions'])
	ann$Group <-  ifelse(ann$Group == target, 'G1', 'G2')
	rownames(ann) <- colnames(data2heat)
	Group         <- c(colors['G1'], colors['G2'])
	# names(Group)  <- c(levels(ann$Group))
	names(Group) <- c('G1', 'G2')
	anno_colors   <- list(Group = Group)

	HM <- pheatmap( data2heat, cluster_rows = T, treeheight_row = 0, annotation_col = ann, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", cellheight= 4,cellwidth= 2, annotation_colors = anno_colors, show_colnames = F, main = paste0('Heatmap between ',target,' classified as ',target,' (G1)\n and ',target,' classified as ',obj,' (G2)'), fontsize = 8,fontsize_row=4)

	plot_dims <- get_plot_dims(HM)

	pdf(paste0(plots_path,dataSet,'_',target,'Vs',obj,'_HM.pdf'), height = plot_dims$height, width = plot_dims$width)
	# print(HM$tree_col)
	print(HM)
	dev.off()