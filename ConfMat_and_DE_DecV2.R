
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


gencode22 <- read.table('/home/sevastopol/data/gserranos/UTILS/gencode.v22.annotation', header=T)
gencode22$gene_id <- str_extract(gencode22$gene_id, 'ENSG[R]?[0-9]+')


process_CM <- function(data){
	b <- apply(data,2, FUN=function(x) (x/sum(x)))
	b[is.nan(b)] <- NA
	b[b==0] <- NA
	b <- b[,colSums(is.na(b))<nrow(b)]
	b <- b[rowSums(is.na(b))<ncol(b),]
	return(b)
}

pd <- import("pandas")

DATASETS <- list.files('/home/sevastopol/data/gserranos/JIND_DE/Data/')

color_red_green <- colorRampPalette(c('#cc5a4e', '#1dcc4b'))(50)
for (dataSet in DATASETS){
	print(dataSet)
	# df2 <- pd$read_pickle(paste0('/home/sevastopol/data/gserranos/JIND_DE/Data/',dataSet,'/train.pkl'))

	annotation <- pd$read_pickle(paste0('/home/sevastopol/data/gserranos/JIND_DE/Data/',dataSet,'/JIND_assignmentbrftune.pkl'))
	annotation$cell_names <- rownames(annotation)

	annotation_seurat <- pd$read_pickle(paste0('/home/sevastopol/data/gserranos/JIND_DE/Data/',dataSet,'/seurat_assignment.pkl'))
	annotation_seurat$cell_names <- rownames(annotation_seurat)

	data <- as.data.frame.matrix(table(annotation$raw_predictions, annotation$labels))
	b <- process_CM(data)
	jind_raw <- pheatmap(b, color=color_red_green, cluster_rows = F, cluster_cols=F, display_numbers=T, cellheight=25, cellwidth=25, na_col= '#61312c', main=paste0('JIND_RAW_', dataSet), legend=F)
	file_xlsx <- paste0('/home/sevastopol/data/gserranos/JIND_DE/Plots/',dataSet,'_CM.xlsx' )
	if (file.exists(file_xlsx)) {
  		file.remove(file_xlsx)
	}
	b[is.na(b)] <-0
	write.xlsx(b, file=file_xlsx, sheetName='JIND_raw', row.names = TRUE, append=TRUE)


	data <- as.data.frame.matrix(table(annotation$predictions, annotation$labels))
	b <- process_CM(data)
	jind <- pheatmap(b, color=color_red_green, cluster_rows = F, cluster_cols=F, display_numbers=T, cellheight=25, cellwidth=25, na_col= '#61312c', main=paste0('JIND_', dataSet), legend=F)
	b[is.na(b)] <-0
	write.xlsx(b, file=file_xlsx, sheetName='JIND', row.names = TRUE, append=TRUE)

	data_seurat <- as.data.frame.matrix(table(annotation_seurat$raw_predictions, annotation_seurat$labels))
	b_seurat <- process_CM(data_seurat)
	seurat_raw <- pheatmap(b_seurat, color=color_red_green, cluster_rows = F, cluster_cols=F, display_numbers=T, cellheight=25, cellwidth=25, na_col= '#61312c', main=paste0('SEURAT_RAW_', dataSet), legend=F)
	b_seurat[is.na(b_seurat)] <-0
	write.xlsx(b_seurat, file=file_xlsx, sheetName='seurat', row.names = TRUE, append=TRUE)


	data_seurat <- as.data.frame.matrix(table(annotation_seurat$predictions, annotation_seurat$labels))
	b_seurat <- process_CM(data_seurat)
	seurat <- pheatmap(b_seurat, color=color_red_green, cluster_rows = F, cluster_cols=F, display_numbers=T, cellheight=25, cellwidth=25, na_col= '#61312c', main=paste0('SEURAT_', dataSet), legend=T)
	b_seurat[is.na(b_seurat)] <-0
	write.xlsx(b_seurat, file=file_xlsx, sheetName='seurat_raw', row.names = TRUE, append=TRUE)


	plot_dims <- get_plot_dims(seurat)
	dev.off()
	pdf(paste0('./Plots/',dataSet,'_CM.pdf'), height = plot_dims$height *2  , width = plot_dims$width*2)
	grid.arrange(grobs = list(jind_raw[[4]], jind[[4]], seurat_raw[[4]], seurat[[4]]), ncol=2)
	dev.off()

}




dataSet <- 'Pancreas_01'

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

colors = c(
"acinar"  = rgb(140, 86, 75, max=255), #BROWN
"endothelial" = rgb(255, 127, 14, max=255), #ORANGE
"Monocyte_FCGR3A" =  rgb(214, 39, 40, max=255), #RED
"ductal" =rgb(44, 160, 44, max=255), #GREEN
"CD4_T_cell" = rgb(31, 119, 180, max=255), #BLUE
"18_Plasma" = rgb(227, 119, 194, max=255), # PINK
"Monocyte_CD14" = rgb(148, 103, 189, max=255), #PURPLE
"Megakaryocyte" = rgb(127, 127, 127, max=255), #GREY
"Unassigned"    = rgb(188, 189, 34, max=255)) # LIME




target <- 'ductal'
obj    <- 'acinar'

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
# data2heat <- data_tmp[rownames(data_tmp) %in%  tmp[tmp$adj.P.Val<0.05 & abs(tmp$logFC)>0.5,'gene_name'],]
data2heat[data2heat > 5] <- 5

ann <- data.frame(Group = annotation[rownames(annotation) %in% colnames(data2heat), 'predictions'])
ann$Group <-  ifelse(ann$Group == target, 'G1', 'G2')
rownames(ann) <- colnames(data2heat)
Group         <- c(colors[target], colors[obj])
# names(Group)  <- c(levels(ann$Group))
names(Group) <- c('G1', 'G2')
anno_colors   <- list(Group = Group)


HM <- pheatmap( data2heat, cluster_rows = T, treeheight_row = 0, annotation_col = ann, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", cellheight= 4,cellwidth= 2, annotation_colors = anno_colors, show_colnames = F, main = paste0('Heatmap between ',target,' classified as ',target,' (G1)\n and ',target,' classified as ',obj,' (G2)'), fontsize = 8,fontsize_row=4)

plot_dims <- get_plot_dims(HM)

pdf(paste0('./Plots/',dataSet,'_',target,'Vs',obj,'_HM.pdf'), height = plot_dims$height, width = plot_dims$width)
	print(HM)
dev.off()


library(dendsort)
callback = function(hc, mat)
{
dendsort(hc, isReverse=T)

}
HM <- pheatmap( data2heat, cluster_rows = T, treeheight_row = 0, annotation_col = ann, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", cellheight= 4,cellwidth= 2, annotation_colors = anno_colors, show_colnames = F, main = paste0('Heatmap between ',target,' classified as ',target,' (G1)\n and classified as ',obj,' (G2)'), fontsize = 8,fontsize_row=4, clustering_callback=callback)

plot_dims <- get_plot_dims(HM)

pdf(paste0('./Plots/',dataSet,'_',target,'Vs',obj,'_HM2.pdf'), height = plot_dims$height, width = plot_dims$width)
	print(HM)
dev.off()




# Acinar stuff

pdf('./Plots/Ductal.pdf', width=20)
plot(HM$tree_col)
dev.off()



clusters <- cutree(HM$tree_col, k = 2)
all_data2 <- all_data[, colnames(all_data) %in% annotation[annotation$labels %in% c('acinar', 'ductal'), 'cell_names']]
tsne_out <- Rtsne(as.matrix(t(all_data2)))
plotter <-  as.data.frame(tsne_out$Y)

plotter$labels <- annotation[rownames(annotation) %in% colnames(all_data2), 'labels']


labels <- data.frame(cell_id = colnames(all_data2) )
labels$group <- 'G0'
labels[labels$cell_id %in% G1,'group'] <- 'G1'
labels[labels$cell_id %in% G2,'group'] <- 'G2_1'
labels[labels$cell_id %in% names(which(clusters == 2)),'group'] <- 'G2_2'
plotter$shape <- labels$group

colors_2_pop <- 
c("alpha" = 'grey',
"ductal" = rgb(44, 160, 44, max=255),
"endothelial" = 'grey',
"delta" = 'grey',
"acinar" = rgb(140, 86, 75, max=255),
"beta" = 'grey',
"gamma" = 'grey',
"epsilon" = 'grey')



euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

ration_dist <- function(a, b, c) {
	return(as.character(euc.dist(a,b)/euc.dist(a,c)))
}


G2_2_mean  <- c(mean(plotter[plotter$shape %in% c('G2_2'), 'V1']), mean(plotter[plotter$shape %in% c('G2_2'), 'V2']))
G2_1_mean  <- c(mean(plotter[plotter$shape %in% c('G2_1'), 'V1']), mean(plotter[plotter$shape %in% c('G2_1'), 'V2']))

G0_mean  <- c(mean(plotter[plotter$shape %in% c('G0'), 'V1']), mean(plotter[plotter$shape %in% c('G0'), 'V2']))
G1_mean  <- c(mean(plotter[plotter$shape %in% c('G1'), 'V1']), mean(plotter[plotter$shape %in% c('G1'), 'V2']))

ductal_mean <- c(mean(plotter[plotter$shape %in% c('G2_2', 'G2_1', 'G1'), 'V1']), mean(plotter[plotter$shape %in% c('G2_2', 'G2_1', 'G1'), 'V2']))
acinar_mean <- c(mean(plotter[plotter$shape %in% c('G0'), 'V1']), mean(plotter[plotter$shape %in% c('G0'), 'V2']))

ration_dist(G2_2_mean,  acinar_mean, ductal_mean)
ration_dist(G2_1_mean,  acinar_mean, ductal_mean)

ration_dist(G0_mean,  acinar_mean, ductal_mean)
ration_dist(G1_mean,  acinar_mean, ductal_mean)




colors_2_pop <- 
c("G0" = rgb(147, 147, 147, alpha=50, max=255),
"G1" = rgb(44, 160, 44, max=255),
"G2_1" = rgb(140, 86, 75, max=255),
"G2_2" = rgb(31, 119, 180, max=255))

pdf('./Plots/Ductal_Vs_Acinar_TSNE_2Pops5.pdf')
ggplot(plotter, aes(x=V1, y=V2, color = shape, label = labels)) + theme_bw() + geom_text()  +scale_color_manual(values = colors_2_pop) + ggtitle(
	paste0('G0 ratio dist.', ration_dist(G0_mean,  acinar_mean, ductal_mean), '\n',
	'G1 ratio dist.', ration_dist(G1_mean,  acinar_mean, ductal_mean), '\n',
	'G2_1 ratio dist. ' , ration_dist(G2_1_mean,  acinar_mean, ductal_mean), '\n ', 
	'G2_2 ratio dist. ' , ration_dist(G2_2_mean,  acinar_mean, ductal_mean), '\n ',
	'Ratio:  Acinar/Ductal'))
dev.off()






# Monocyte stuff

pdf('./Plots/Monocyte_FCGR3A.pdf', width=20)
plot(HM$tree_col)
dev.off()



clusters <- cutree(HM$tree_col, k = 2)

cluster_names <- HM$tree_col$labels
cluster_pos <- HM$tree_col$order

#####################
# lost code 
cluster_names[cluster_pos]
which(cluster_names[cluster_pos] %in% G2)
split(dat, cumsum(c(1, diff(dat) != 1)))
#####################


all_data2 <- all_data[, colnames(all_data) %in% annotation[annotation$labels %in% c('Monocyte_CD14', 'Monocyte_FCGR3A'), 'cell_names']]
tsne_out <- Rtsne(as.matrix(t(all_data2)))
plotter <-  as.data.frame(tsne_out$Y)

plotter$labels <- annotation[rownames(annotation) %in% colnames(all_data2), 'labels']


labels <- data.frame(cell_id = colnames(all_data2) )
labels$group <- 'G0'
labels[labels$cell_id %in% G1,'group'] <- 'G1'
labels[labels$cell_id %in% G2,'group'] <- 'G2_1'
labels[labels$cell_id %in% names(which(clusters == 2)),'group'] <- 'G2_2'
plotter$shape <- labels$group



euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

ration_dist <- function(a, b, c) {
	return(as.character(euc.dist(a,b)/euc.dist(a,c)))
}


G2_2_mean  <- c(mean(plotter[plotter$shape %in% c('G2_2'), 'V1']), mean(plotter[plotter$shape %in% c('G2_2'), 'V2']))
G2_1_mean  <- c(mean(plotter[plotter$shape %in% c('G2_1'), 'V1']), mean(plotter[plotter$shape %in% c('G2_1'), 'V2']))

G0_mean  <- c(mean(plotter[plotter$shape %in% c('G0'), 'V1']), mean(plotter[plotter$shape %in% c('G0'), 'V2']))
G1_mean  <- c(mean(plotter[plotter$shape %in% c('G1'), 'V1']), mean(plotter[plotter$shape %in% c('G1'), 'V2']))

ductal_mean <- c(mean(plotter[plotter$shape %in% c('G2_2', 'G2_1', 'G1'), 'V1']), mean(plotter[plotter$shape %in% c('G2_2', 'G2_1', 'G1'), 'V2']))
acinar_mean <- c(mean(plotter[plotter$shape %in% c('G0'), 'V1']), mean(plotter[plotter$shape %in% c('G0'), 'V2']))

ration_dist(G2_2_mean,  acinar_mean, ductal_mean)
ration_dist(G2_1_mean,  acinar_mean, ductal_mean)

ration_dist(G0_mean,  acinar_mean, ductal_mean)
ration_dist(G1_mean,  acinar_mean, ductal_mean)




colors_2_pop <- 
c("G0" = rgb(147, 147, 147, alpha=50, max=255),
"G1" = rgb(44, 160, 44, max=255),
"G2_1" = rgb(140, 86, 75, max=255),
"G2_2" = rgb(31, 119, 180, max=255))

pdf('./Plots/Ductal_Vs_Acinar_TSNE_2Pops4.pdf', 20, 20)
ggplot(plotter, aes(x=V1, y=V2, color = shape, label = labels)) + theme_bw() + geom_text()  +scale_color_manual(values = colors_2_pop) + ggtitle(
	paste0('G2_1 ratio dist. ' , ration_dist(G2_1_mean,  acinar_mean, ductal_mean), '  // ', 
	'G2_2 ratio dist. ' , ration_dist(G2_2_mean,  acinar_mean, ductal_mean)))
dev.off()
