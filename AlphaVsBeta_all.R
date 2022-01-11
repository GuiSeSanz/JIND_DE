library(reticulate)
use_condaenv(condaenv = "jind", conda = "auto", required = TRUE)
library(argparse)
library(stringr)
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
library(wesanderson)

callback = function(hc, mat)
{
  dendsort(hc, isReverse=T)
  
}


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

get_plot_dims <- function(heat_map){
	plot_height <- sum(sapply(heat_map$gtable$heights, grid::convertHeight, "in"))
	plot_width  <- sum(sapply(heat_map$gtable$widths, grid::convertWidth, "in"))
	dev.off()
	return(list(height = plot_height, width = plot_width))
}

get_topk_features <- function(data, k = 5000){
	print(paste0("Reudcing to top", k, " features"))
	var_ <- apply(data, 1, var)
	names(var_) <- rownames(data)
	genes_2_keep <- names(sort(var_, decreasing=T)[1:5000])
	red_data <- data[rownames(data) %in% genes_2_keep, ]
	return(red_data)
}

perform_DE <- function(data, group){
  ### DE
  # create the linear model
  fit_tmp <- lmFit(data, group)
  # model correction
  fit_tmp <- eBayes(fit_tmp)
  # results <- topTable(fit_tmp, n=Inf)
  x <- paste0('G1', '-', 'G2')
  contrast_mat_tmp <- makeContrasts(contrasts=x, levels= c('G1', 'G2'))
  fit2_tmp <- contrasts.fit(fit_tmp, contrast_mat_tmp)
  fit2_tmp <- eBayes(fit2_tmp)
  tmp   <- topTable(fit2_tmp, adjust="fdr", n=Inf)
  tmp$gene_name <- rownames(tmp)
  
  tmp[tmp$P.Value == 0, 'P.Value'] <- 1.445749e-281
  tmp[tmp$adj.P.Val == 0, 'adj.P.Val'] <- 1.445749e-281
  
  return(tmp)
}



parser <- ArgumentParser(description='alpha beta DE')
parser$add_argument('--test_path', default="./Data/test_indexUnified.pkl", type="character",
                    help='path to train data frame with labels')
parser$add_argument('--assignment', default="./Data/Pancreas_01/JIND_assignmentbrftune.pkl", type="character",
                    help='path to assignment data frame with labels')
parser$add_argument('--method', default="JIND", type="character",
                    help='method name')

args <- parser$parse_args()




pd <- import("pandas")


df <- pd$read_pickle(args$test_path)

jind_labels <- pd$read_pickle(args$assignment)

jind_labels$cell_names <- rownames(df)
all_data <- t(df[, -which(colnames(df) %in% c('labels'))])
stand_dev <- apply(all_data,1, sd)
data_tmp <- all_data[!rownames(all_data) %in% names(which(stand_dev ==0)), ]

df$labels <- as.character(df$labels)
jind_labels$labels <- as.character(df$labels)

table(jind_labels$labels, jind_labels$predictions)
table(jind_labels[jind_labels$predictions != jind_labels$labels, 'labels'])
table(jind_labels[jind_labels$labels == 'microglia', 'predictions'])

G1 <- jind_labels[jind_labels$labels == 'alpha' &  jind_labels$predictions == 'alpha', 'cell_names']
# G1 <- jind_labels[jind_labels$labels == 'alpha' &  jind_labels$predictions == 'alpha', 'cellname']
G2 <- jind_labels[jind_labels$labels == 'alpha' &  jind_labels$raw_predictions == 'beta', 'cell_names']
# G2 <- jind_labels[jind_labels$labels == 'alpha' &  jind_labels$raw_predictions == 'beta', 'cellname']


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
tmp <- perform_DE(data_tmp, design_tmp)

tmp$logpval <- -log(tmp$P.Value)
tmp2 <- tmp[ order(-tmp$logpval), c('gene_name', 'logFC', 'P.Value' ,'logpval', 'adj.P.Val')]

data2heat <- data_tmp[tmp2[tmp2$adj.P.Val<0.05 & abs(tmp2$logFC)>2,'gene_name'],]
dim(data2heat)
rownames(data2heat) <- toupper(rownames(data2heat))
ann <- data.frame(Group = colnames(data2heat))
ann$Group <-  ifelse(ann$Group %in% G1, 'G1', 'G2')
rownames(ann) <- colnames(data2heat)
Group <- ggthemes::tableau_color_pal('Classic 10')(2)
names(Group)  <- c(unique(ann$Group))
data2xlsx <- data2heat
ann_tmp <- as.data.frame(ann)
data2xlsx <- t(merge(ann_tmp, t(data2xlsx), by=0))
xlsx::write.xlsx(tmp2, file = paste0('./Plots/', args$method, '_DE_results.xlsx'), row.names=FALSE, sheetName = 'JIND+')
xlsx::write.xlsx(tmp2[tmp2$adj.P.Val<0.05 & abs(tmp2$logFC)>2,], file = paste0('./Plots/', args$method, '_DE_results.xlsx'), row.names=FALSE, sheetName = 'JIND+ (FDR < 0.05 & logFC>2)', append = TRUE)

anno_colors <- list(Group = Group)
HM <- pheatmap( data2heat, scale='row', color=viridis::viridis(10), 
				cluster_rows = T, treeheight_row = 0, annotation_col = ann,  
				annotation_colors = anno_colors, clustering_distance_rows = "euclidean", 
				clustering_distance_cols = "euclidean", show_colnames = F, 
				main = paste0('Heatmap between ', 'alpha',' (G1)\n and ','unassigned alpha',' (G2)'), 
				fontsize = 8,fontsize_row=6 , silent=TRUE, legend=FALSE) #, callback = callback)
plot_dims <- get_plot_dims(HM)

genes_2_plot <- reshape2::melt(data2heat)
genes_2_plot <- merge(genes_2_plot, ann, by.x='Var2', by.y =0)



alpha <- jind_labels[jind_labels$labels == 'alpha' , 'cell_names']
beta <- jind_labels[jind_labels$labels == 'beta', 'cell_names']
selection_AB <- c(alpha, beta)
data_tmp <- all_data[ , colnames(all_data) %in%  selection_AB]
stand_dev <- apply(data_tmp,1, sd)
data_tmp <- data_tmp[!rownames(data_tmp) %in% names(which(stand_dev ==0)), ]
DESIGN <- data.frame(cells = colnames(data_tmp))
labels <- c()
for (cell in DESIGN$cells){
	ifelse(cell %in% G1, labels <- c(labels, 'alpha'), labels <- c(labels, 'beta'))
}
DESIGN$labels <- labels
design_tmp <- as.matrix(DESIGN[, 'labels'])
design_tmp[design_tmp != 'beta'] <-1
design_tmp[design_tmp == 'alpha'] <-0
design_tmp <- model.matrix(~0+as.factor(design_tmp[,1]))
colnames(design_tmp) <- c('alpha', 'beta')
rownames(design_tmp) <- colnames(data_tmp)


### DE
tmp <- perform_DE(data_tmp, design_tmp)

tmp$logpval <- -log(tmp$P.Value)
tmp2 <- tmp[ order(-tmp$logpval), c('gene_name', 'logFC', 'P.Value' ,'logpval', 'adj.P.Val')]


data2heat_AB <- data_tmp[tmp2[tmp2$adj.P.Val<0.05 & abs(tmp2$logFC)>2,'gene_name'],]
data2heat_AB <- data_tmp[tmp2[1:100, 'gene_name'],]

ann2 <- data.frame(Group = ifelse(colnames(data2heat_AB) %in% alpha, 'Alpha', 'Beta'))
ann2[['Prediction']] <- ifelse(!colnames(data2heat_AB) %in% selection, 'Beta', ifelse(colnames(data2heat_AB) %in% G1, 'G1', 'G2'))
rownames(ann2) <- colnames(data2heat_AB)

anno_colors2 <- list(Group=setNames(ggthemes::tableau_color_pal('Classic 10')(4)[3:4], c('Beta', 'Alpha')),
					 Prediction=setNames(ggthemes::tableau_color_pal('Classic 10')(3), c('G1', 'G2', 'Beta')))

HM_all <- pheatmap( data2heat_AB, scale='row', color=viridis::viridis(10), cluster_rows = T, treeheight_row = 0,treeheight_col = 5, annotation_col = ann2,  annotation_colors = anno_colors2,clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", show_colnames = F, show_rownames = F, main = paste0('Heatmap between ', 'alpha',' (G1)\n and ','beta',' (G2)'), fontsize = 8,fontsize_row=6 ,silent=TRUE,legend=FALSE)
# pdf('./Plots/Test.pdf', 12, 10)
# print(HM_all)
# dev.off()





tsne_out <- Rtsne(as.matrix(t(all_data)))
plotter <-  as.data.frame(tsne_out$Y)
plotter$cell_names <- colnames(all_data)
plotter <- merge(plotter, jind_labels[, c('cell_names', 'labels', 'raw_predictions')], by='cell_names')
plotter$Cluster <- ifelse(plotter$cell_names %in% jind_labels[jind_labels$labels == 'beta', 'cell_names'],'Beta',
				   ifelse(plotter$cell_names %in% G1, 'G1' , 
				   ifelse(plotter$cell_names %in% G2, 'G2' , 'Other')))
plotter$transparency <- ifelse(plotter$Cluster %in% c('G1', 'G2', 'Beta'), TRUE, FALSE)

pdf(paste0('./Plots/Test_AlphaVsBeta_',args$method, '.pdf'), 12, 10)
cowplot::plot_grid(
	cowplot::plot_grid(
		cowplot::plot_grid(
			HM$gtable,
			ggplot(plotter, aes(x=V1, y=V2, color = Cluster, label = labels, alpha=transparency)) + theme_classic() + geom_text()  +scale_color_manual(values = c('Other' ='808080',setNames(ggthemes::tableau_color_pal('Classic 10')(3), c('G1', 'G2', 'Beta')))) + theme(legend.text=element_text(size=10), axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank()) + scale_alpha_discrete(range = c(0.2, 0.9), guide='none'),
		nrow=2, rel_heights=c(2,3)),

	ggplot(genes_2_plot, aes(x=Group, y=value, fill = Group)) + geom_violin(alpha=0.7) + scale_fill_manual(values = Group) + ggtitle('title') + theme(legend.text=element_text(size=15), axis.ticks=element_blank(),axis.title.x=element_blank(),axis.text=element_blank())+ ggprism::theme_prism()  + facet_wrap(~Var1) + ggsignif::geom_signif(comparisons = list(c("G1", "G2")), map_signif_level = TRUE, vjust=0.5),
	nrow=2, rel_heights=c(3,2)),
	HM_all$gtable,
ncol=2, rel_heights=c(3,2))
dev.off()


