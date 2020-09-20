library(SummarizedExperiment)
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



"multi" = function(type, x, xlim, title1, title2, title3, ...) {

    if (type == "density")
      multidensity(x,
                   xlim = xlim,
                   main = "",
                   xlab = "",
                   ylab = "", ...)
    if (type == "ecdf")
      multiecdf(x,
                xlim = xlim,
                main = "",
                xlab = "",
                ylab = "", ...)
    mtext(title1, side = 2, adj = 0.5, padj = -4 , cex = 0.7)
    mtext(title2, side = 1, adj = 0.5, padj = 4 , cex = 0.7)
    mtext(title3, side = 3, adj = 0.5, padj = -1 , cex = 0.7)
     
}

"detOutliers" <- function(data, typeSamples, filenamePdf, W, H) { 

	colourRange = rgb(seq(0,1,l=256),seq(0,1,l=256),seq(1,0,l=256))
	
	sArray = apply(data, 1, median)
	dat <- data
	for (i in 1:ncol(data))
		dat[,i] <- dat[,i]-sArray
		
	outM = as.dist(dist2(na.omit(dat)))
		
	d.row = as.dendrogram(hclust(outM))
	od.row = order.dendrogram(d.row)
	m = as.matrix(outM)
	namesCol <- colnames(data)
	colnames(m) = namesCol
	rownames(m) = namesCol
   
	covar = typeSamples
	lev = levels(as.factor(covar))
	corres = matrix(0,nrow=length(lev),ncol=2)
	colourCov = brewer.pal(12,"Set3")

	pdf(file = filenamePdf, width = W, heigh = H, colormode = "rgb")
		print(levelplot(m[od.row,od.row],
			scales=list(x=list(rot=90)),
			legend=list(
			top=list(fun=dendrogramGrob,args=list(x=d.row,side="top")),
			right=list(fun=dendrogramGrob,args=list(x=d.row,side="right", size.add= 1, add = list(rect = list(col = "transparent", fill = colourCov[as.factor(covar)])), type = "rectangle"))),
			colorkey = list(space ="left"),
			xlab="",ylab="",
			col.regions=colourRange))
		
		x=0.06
		y=0.98
		
		for(i in 1:length(lev))
		{
		corres[i,] = c(unique(covar[covar == lev[i]]),colourCov[i])
		grid.text(lev[i], x=x, y=y, just="left")
		grid.rect(gp=gpar(fill=corres[i,2],col="transparent"), x=x-0.02, y=y, width=0.02, height=0.02)
		y=y-0.03
		}
	dev.off()

	return(outM)
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

gencode22 <- read.table('/home/sevastopol/data/gserranos/UTILS/gencode.v22.annotation', header=T)
gencode22$gene_id <- str_extract(gencode22$gene_id, 'ENSG[R]?[0-9]+')



# Load Data


pd <- import("pandas")
df <- pd$read_pickle('test.pkl')

annotation <- pd$read_pickle('JIND_assignmentbrftune.pkl')
annotation$cell_names <- rownames(annotation)

all_data <- t(df[, -which(colnames(df) %in% c('labels'))])
stand_dev <- apply(all_data,1, sd)
all_data <- all_data[!rownames(all_data) %in% names(which(stand_dev ==0)), ]

tsne_out <- Rtsne(as.matrix(t(all_data)))
plotter <-  as.data.frame(tsne_out$Y)

annotation$labels <- gsub(' ', '_', annotation$labels)
annotation$predictions <- gsub(' ', '_', annotation$predictions)
annotation$raw_predictions <- gsub(' ', '_', annotation$raw_predictions)

# All CD4 Vs All CD8

target <- 'CD4_T_cell'
obj    <- 'CD8_T_cell'

G1 <- annotation[annotation$labels == target, 'cell_names']
G2 <- annotation[annotation$labels == obj, 'cell_names']

selection <- c(G1, G2)
data_tmp <- all_data[ , colnames(all_data) %in%  selection]


DESIGN <- data.frame(cells = rownames(df), label = gsub(' ', '_',df$labels))

design_tmp <- as.matrix(DESIGN[DESIGN$cells %in% selection, 'label'])

design_tmp[design_tmp != target] <-1
design_tmp[design_tmp == target] <-0
design_tmp <- model.matrix(~0+as.factor(design_tmp[,1]))
colnames(design_tmp) <- c(target, obj)
rownames(design_tmp) <- DESIGN[DESIGN$cells %in% selection, 'cells']


### DE
# create the linear model
fit_tmp <- lmFit(data_tmp, design_tmp)
# model correction
fit_tmp <- eBayes(fit_tmp)
# results <- topTable(fit_tmp, n=Inf)
x <- paste0(target, '-', obj)
contrast_mat_tmp <- makeContrasts(contrasts=x, levels= c(target, obj))
fit2_tmp <- contrasts.fit(fit_tmp, contrast_mat_tmp)
fit2_tmp <- eBayes(fit2_tmp)
tmp   <- topTable(fit2_tmp, adjust="BH", n=Inf)
tmp$gene_name <- rownames(tmp)
tmp <- merge(tmp, gencode22, by= 'gene_name')

pdf('./Plots/Test_gene.pdf')
plt <- data.frame(Value = all_data['CD8A', colnames(all_data) %in% G2], Cell = 'G2') 
plt <- rbind(plt, data.frame(Value = all_data['CD8A', colnames(all_data) %in% G1], Cell = 'G1'))
ggplot(plt,  aes(x=Cell, y=Value)) + geom_boxplot()
dev.off()


tmp[tmp$P.Value == 0, 'P.Value'] <- 1.445749e-281
tmp[tmp$adj.P.Val == 0, 'P.Value'] <- 1.445749e-281

pdf('./Plots/TT_cd4Vscd8_VP.pdf')
	graphContrast(tmp," (Ctrl B > 0, FC>0.5)", 0, 0.5, 1)
dev.off()

tmp$logpval <- -log(tmp$P.Value)
tmp2 <- tmp[ order(-tmp$logpval), c('gene_name', 'logFC', 'P.Value' ,'logpval')]

saveRDS(tmp2, './Data/TT_cd4Vscd8.rds')

data2heat <- data_tmp[rownames(data_tmp) %in%  tmp2[1:50, 'gene_name'],]
data2heat[data2heat > 5] <- 5

ann <- data.frame(Var1 = annotation[rownames(annotation) %in% colnames(data2heat), 'labels'])
rownames(ann) <- colnames(data2heat)
Var1        <- c("#b577a8", "#979858")
names(Var1) <- c(levels(ann$Var1))
anno_colors <- list(Var1 = Var1)

pdf('./Plots/TT_cd4Vscd8_HM.pdf')
pheatmap( data2heat, cluster_rows = T, annotation = ann,  annotation_colors = anno_colors, show_colnames = F)
dev.off()



# unique(annotation$labels)
GREY <- '#74797a'
colors <- c(
"NK_cell"  = GREY,
"Monocyte_CD14" = GREY,
"CD8_T_cell" = '#59a1bd',
"B_cell" = GREY,
"CD4_T_cell" = '#a7993d',
"Plasmacytoid_dendritic_cell" = GREY,
"Monocyte_FCGR3A" = GREY,
"Megakaryocyte" = GREY,
"Hematopoietic_stem_cell"    = GREY
)


plotter$labels <- annotation[rownames(annotation) %in% colnames(all_data), 'labels']

pdf('./Plots/TT_cd4Vscd8_TSNE.pdf')
ggplot(plotter, aes(x=V1, y=V2, color = labels)) + theme_bw() +geom_point(size = 1.5) +scale_color_manual(values = colors)
dev.off()




# Real CD4 Vs CD4 predicted as CD8

target <- 'CD4_T_cell'
obj    <- 'CD8_T_cell'

G1 <- annotation[annotation$labels == target & annotation$prediction == target, 'cell_names']
G2 <- annotation[annotation$labels == target & annotation$prediction == obj, 'cell_names']

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
tmp   <- topTable(fit2_tmp, adjust="BH", n=Inf)
tmp$gene_name <- rownames(tmp)
tmp <- merge(tmp, gencode22, by= 'gene_name')



tmp[tmp$P.Value == 0, 'P.Value'] <- 1.445749e-281
tmp[tmp$adj.P.Val == 0, 'P.Value'] <- 1.445749e-281


pdf('./Plots/TT_cd4realVscd4_classCd8_VP.pdf')
	graphContrast(tmp," (Ctrl B > 0, FC>0.5)", 0, 0.5, 1)
dev.off()

tmp$logpval <- -log(tmp$P.Value)
tmp2 <- tmp[ order(-tmp$logpval), c('gene_name', 'logFC', 'P.Value' ,'logpval')]

saveRDS(tmp2, './Data/TT_cd4realVscd4_classCd8.rds')


data2heat <- data_tmp[rownames(data_tmp) %in%  tmp2[1:20, 'gene_name'],]
data2heat[data2heat > 5] <- 5

ann <- data.frame(Var1 = annotation[rownames(annotation) %in% colnames(data2heat), 'predictions'])
rownames(ann) <- colnames(data2heat)
Var1        <- c("#b577a8", "#979858")
names(Var1) <- c(levels(ann$Var1))
anno_colors <- list(Var1 = Var1)

pdf('./Plots/TT_cd4realVscd4_classCd8_HM.pdf')
pheatmap( data2heat, cluster_rows = T, annotation = ann, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", annotation_colors = anno_colors, show_colnames = F)
dev.off()



# unique(annotation$labels)
GREY <- '#74797a'
TARGET <- '#cc5f43'
colors <- c(
"TARGET" = TARGET,
"NK_cell"  = GREY,
"Monocyte_CD14" = GREY,
"CD8_T_cell" = '#59a1bd',
"B_cell" = GREY,
"CD4_T_cell" = '#a7993d',
"Plasmacytoid_dendritic_cell" = GREY,
"Monocyte_FCGR3A" = GREY,
"Megakaryocyte" = GREY,
"Hematopoietic_stem_cell"    = GREY
)

ann_color <- annotation[rownames(annotation) %in% colnames(all_data), c('cell_names', 'labels')]
ann_color[ann_color$cell_names %in% G2, 'labels'] <- 'TARGET'


plotter$labels <- ann_color$labels

pdf('./Plots/TT_cd4realVscd4_classCd8_TSNE.pdf')
ggplot(plotter, aes(x=V1, y=V2, color = labels)) + theme_bw() +geom_point(size = 1.5) +scale_color_manual(values = colors)
dev.off()






# Real CD4 Vs CD8 predicted as CD4



target <- 'CD4_T_cell'
obj    <- 'CD8_T_cell'

G1 <- annotation[annotation$labels == target & annotation$prediction == target, 'cell_names']
G2 <- annotation[annotation$labels == obj & annotation$prediction == target, 'cell_names']

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
tmp   <- topTable(fit2_tmp, adjust="BH", n=Inf)
tmp$gene_name <- rownames(tmp)
tmp <- merge(tmp, gencode22, by= 'gene_name')


tmp[tmp$P.Value == 0, 'P.Value'] <- 1.445749e-281
tmp[tmp$adj.P.Val == 0, 'P.Value'] <- 1.445749e-281


pdf('./Plots/TT_cd4realVscd8_classCd4_VP.pdf')
	graphContrast(tmp," (Ctrl B > 0, FC>0.5)", 0, 0.5, 1)
dev.off()

tmp$logpval <- -log(tmp$P.Value)
tmp2 <- tmp[ order(-tmp$logpval), c('gene_name', 'logFC', 'P.Value' ,'logpval')]

saveRDS(tmp2, './Data/TT_cd4realVscd8_classCd4.rds')


data2heat <- data_tmp[rownames(data_tmp) %in%  tmp2[1:20, 'gene_name'],]
data2heat[data2heat > 5] <- 5

ann <- data.frame(Var1 = annotation[rownames(annotation) %in% colnames(data2heat), 'labels'])
rownames(ann) <- colnames(data2heat)
Var1        <- c("#b577a8", "#979858")
names(Var1) <- c(levels(ann$Var1))
anno_colors <- list(Var1 = Var1)

pdf('./Plots/TT_cd4realVscd8_classCd4_HM.pdf')
pheatmap( data2heat, cluster_rows = T, annotation = ann, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", annotation_colors = anno_colors, show_colnames = F)
dev.off()



# unique(annotation$labels)
GREY <- '#74797a'
TARGET <- '#cc5f43'
colors <- c(
"TARGET" = TARGET,
"NK_cell"  = GREY,
"Monocyte_CD14" = GREY,
"CD8_T_cell" = '#59a1bd',
"B_cell" = GREY,
"CD4_T_cell" = '#a7993d',
"Plasmacytoid_dendritic_cell" = GREY,
"Monocyte_FCGR3A" = GREY,
"Megakaryocyte" = GREY,
"Hematopoietic_stem_cell"    = GREY
)

ann_color <- annotation[rownames(annotation) %in% colnames(all_data), c('cell_names', 'labels')]
ann_color[ann_color$cell_names %in% G2, 'labels'] <- 'TARGET'


plotter$labels <- ann_color$labels

pdf('./Plots/TT_cd4realVscd8_classCd4_TSNE.pdf')
ggplot(plotter, aes(x=V1, y=V2, color = labels)) + theme_bw() +geom_point(size = 1.5) +scale_color_manual(values = colors)
dev.off()



# Real CD84 Vs CD8 predicted as CD4



target <- 'CD4_T_cell'
obj    <- 'CD8_T_cell'

G1 <- annotation[annotation$labels == obj & annotation$prediction == obj, 'cell_names']
G2 <- annotation[annotation$labels == obj & annotation$prediction == target, 'cell_names']

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
tmp   <- topTable(fit2_tmp, adjust="BH", n=Inf)
tmp$gene_name <- rownames(tmp)
tmp <- merge(tmp, gencode22, by= 'gene_name')


tmp[tmp$P.Value == 0, 'P.Value'] <- 1.445749e-281
tmp[tmp$adj.P.Val == 0, 'P.Value'] <- 1.445749e-281


pdf('./Plots/TT_cd48realVscd8_classCd4_VP.pdf')
	graphContrast(tmp," (Ctrl B > 0, FC>0.5)", 0, 0.5, 1)
dev.off()

tmp$logpval <- -log(tmp$P.Value)
tmp2 <- tmp[ order(-tmp$logpval), c('gene_name', 'logFC', 'P.Value' ,'logpval')]

saveRDS(tmp2, './Data/TT_cd48realVscd8_classCd4.rds')


data2heat <- data_tmp[rownames(data_tmp) %in%  tmp2[1:20, 'gene_name'],]
data2heat[data2heat > 5] <- 5

ann <- data.frame(Var1 = annotation[rownames(annotation) %in% colnames(data2heat), 'labels'])
rownames(ann) <- colnames(data2heat)
Var1        <- c("#b577a8", "#979858")
names(Var1) <- c(levels(ann$Var1))
anno_colors <- list(Var1 = Var1)

pdf('./Plots/TT_cd48realVscd8_classCd4_HM.pdf')
pheatmap( data2heat, cluster_rows = T, annotation = ann, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", annotation_colors = anno_colors, show_colnames = F)
dev.off()



# unique(annotation$labels)
GREY <- '#74797a'
TARGET <- '#cc5f43'
colors <- c(
"TARGET" = TARGET,
"NK_cell"  = GREY,
"Monocyte_CD14" = GREY,
"CD8_T_cell" = '#59a1bd',
"B_cell" = GREY,
"CD4_T_cell" = '#a7993d',
"Plasmacytoid_dendritic_cell" = GREY,
"Monocyte_FCGR3A" = GREY,
"Megakaryocyte" = GREY,
"Hematopoietic_stem_cell"    = GREY
)

ann_color <- annotation[rownames(annotation) %in% colnames(all_data), c('cell_names', 'labels')]
ann_color[ann_color$cell_names %in% G2, 'labels'] <- 'TARGET'


plotter$labels <- ann_color$labels

pdf('./Plots/TT_cd48realVscd8_classCd4_TSNE.pdf')
ggplot(plotter, aes(x=V1, y=V2, color = labels)) + theme_bw() +geom_point(size = 1.5) +scale_color_manual(values = colors)
dev.off()




# CD4 Vs CD4 labelled as Unassigned



target <- 'CD4_T_cell'
obj    <- 'CD8_T_cell'

G1 <- annotation[annotation$labels == target & annotation$prediction == target, 'cell_names']
G2 <- annotation[annotation$labels == target & annotation$raw_predictions == target & annotation$prediction == 'Unassigned', 'cell_names']

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
tmp   <- topTable(fit2_tmp, adjust="BH", n=Inf)
tmp$gene_name <- rownames(tmp)
tmp <- merge(tmp, gencode22, by= 'gene_name')


tmp[tmp$P.Value == 0, 'P.Value'] <- 1.445749e-281
tmp[tmp$adj.P.Val == 0, 'P.Value'] <- 1.445749e-281


pdf('./Plots/TT_cd4realVscd4_Unassigned_VP.pdf')
	graphContrast(tmp," (Ctrl B > 0, FC>0.5)", 0, 0.5, 1)
dev.off()

tmp$logpval <- -log(tmp$P.Value)
tmp2 <- tmp[ order(-tmp$logpval), c('gene_name', 'logFC', 'P.Value' ,'logpval')]

saveRDS(tmp2, './Data/TT_cd4realVscd4_Unassigned.rds')


data2heat <- data_tmp[rownames(data_tmp) %in%  tmp2[1:20, 'gene_name'],]
data2heat[data2heat > 5] <- 5

ann <- data.frame(Var1 = annotation[rownames(annotation) %in% colnames(data2heat), 'predictions'])
rownames(ann) <- colnames(data2heat)
Var1        <- c("#b577a8", "#979858")
names(Var1) <- c(levels(ann$Var1))
anno_colors <- list(Var1 = Var1)

pdf('./Plots/TT_cd4realVscd4_Unassigned_HM.pdf')
pheatmap( data2heat, cluster_rows = T, annotation = ann, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", annotation_colors = anno_colors, show_colnames = F)
dev.off()



# unique(annotation$labels)

GREY <- '#74797a'
TARGET <- '#cc5f43'
colors <- c(
"TARGET" = TARGET,
"NK_cell"  = GREY,
"Monocyte_CD14" = GREY,
"CD8_T_cell" = '#59a1bd',
"B_cell" = GREY,
"CD4_T_cell" = '#a7993d',
"Plasmacytoid_dendritic_cell" = GREY,
"Monocyte_FCGR3A" = GREY,
"Megakaryocyte" = GREY,
"Hematopoietic_stem_cell"    = GREY
)

ann_color <- annotation[rownames(annotation) %in% colnames(all_data), c('cell_names', 'labels')]
ann_color[ann_color$cell_names %in% G2, 'labels'] <- 'TARGET'


plotter$labels <- ann_color$labels

pdf('./Plots/TT_cd4realVscd4_Unassigned_TSNE.pdf')
ggplot(plotter, aes(x=V1, y=V2, color = labels)) + theme_bw() +geom_point(size = 1.5) +scale_color_manual(values = colors)
dev.off()


# pdf_files <- list.files('./Plots/', pattern = '^TT')

# convert TT_cd4realVscd4_classCd8_HM.pdf TT_cd4realVscd4_classCd8_HM.png
# convert TT_cd4realVscd4_classCd8_TSNE.pdf TT_cd4realVscd4_classCd8_TSNE.png
# convert TT_cd4realVscd4_classCd8_VP.pdf TT_cd4realVscd4_classCd8_VP.png
# convert TT_cd4realVscd4_Unassigned_HM.pdf TT_cd4realVscd4_Unassigned_HM.png
# convert TT_cd4realVscd4_Unassigned_TSNE.pdf TT_cd4realVscd4_Unassigned_TSNE.png
# convert TT_cd4realVscd4_Unassigned_VP.pdf TT_cd4realVscd4_Unassigned_VP.png
# convert TT_cd4realVscd8_classCd4_HM.pdf TT_cd4realVscd8_classCd4_HM.png
# convert TT_cd4realVscd8_classCd4_TSNE.pdf TT_cd4realVscd8_classCd4_TSNE.png
# convert TT_cd4realVscd8_classCd4_VP.pdf TT_cd4realVscd8_classCd4_VP.png
# convert TT_cd4Vscd8_HM.pdf TT_cd4Vscd8_HM.png
# convert TT_cd4Vscd8_TSNE.pdf TT_cd4Vscd8_TSNE.png
# convert TT_cd4Vscd8_VP.pdf TT_cd4Vscd8_VP.png


# MONOCYTES

target <- 'Monocyte_FCGR3A'
obj    <- 'Monocyte_CD14'

G1 <- annotation[annotation$labels == target, 'cell_names']
G2 <- annotation[annotation$labels == obj, 'cell_names']

selection <- c(G1, G2)

data_tmp <- all_data[ , colnames(all_data) %in%  selection]


DESIGN <- data.frame(cells = rownames(df), label = gsub(' ', '_',df$labels))
annotation$labels <- gsub(' ', '_', annotation$labels)
annotation$predictions <- gsub(' ', '_', annotation$predictions)
annotation$raw_predictions <- gsub(' ', '_', annotation$raw_predictions)
design_tmp <- as.matrix(DESIGN[DESIGN$cells %in% selection, 'label'])

design_tmp[design_tmp != target] <-1
design_tmp[design_tmp == target] <-0
design_tmp <- model.matrix(~0+as.factor(design_tmp[,1]))
colnames(design_tmp) <- c(target, obj)
rownames(design_tmp) <- DESIGN[DESIGN$cells %in% selection, 'cells']


### DE
# create the linear model
fit_tmp <- lmFit(data_tmp, design_tmp)
# model correction
fit_tmp <- eBayes(fit_tmp)
# results <- topTable(fit_tmp, n=Inf)
x <- paste0(target, '-', obj)
contrast_mat_tmp <- makeContrasts(contrasts=x, levels= c(target, obj))
fit2_tmp <- contrasts.fit(fit_tmp, contrast_mat_tmp)
fit2_tmp <- eBayes(fit2_tmp)
tmp   <- topTable(fit2_tmp, adjust="BH", n=Inf)
tmp$gene_name <- rownames(tmp)
tmp <- merge(tmp, gencode22, by= 'gene_name')


tmp[tmp$P.Value == 0, 'P.Value'] <- 1.445749e-281
tmp[tmp$adj.P.Val == 0, 'P.Value'] <- 1.445749e-281


pdf('./Plots/TT_MonFcVsMonCD_VP.pdf')
	graphContrast(tmp," (Ctrl B > 0, FC>0.5)", 0, 0.5, 1)
dev.off()

tmp$logpval <- -log(tmp$P.Value)
tmp2 <- tmp[ order(-tmp$logpval), c('gene_name', 'logFC', 'P.Value' ,'logpval')]

saveRDS(tmp2, './Data/TT_MonFcVsMonCD.rds')

data2heat <- data_tmp[rownames(data_tmp) %in%  tmp2[1:20, 'gene_name'],]
data2heat[data2heat > 5] <- 5

ann <- data.frame(Var1 = annotation[rownames(annotation) %in% colnames(data2heat), 'labels'])
rownames(ann) <- colnames(data2heat)
Var1        <- c("#b577a8", "#979858")
names(Var1) <- c(levels(ann$Var1))
anno_colors <- list(Var1 = Var1)

pdf('./Plots/TT_MonFcVsMonCD_HM.pdf')
pheatmap( data2heat, cluster_rows = T, annotation = ann,  annotation_colors = anno_colors, show_colnames = F)
dev.off()



# unique(annotation$labels)

GREY <- '#74797a'
TARGET <- '#cc5f43'
colors <- c(
"NK_cell"  = GREY,
"Monocyte_CD14" = '#59a1bd',
"CD8_T_cell" = GREY,
"B_cell" = GREY,
"CD4_T_cell" = GREY,
"Plasmacytoid_dendritic_cell" = GREY,
"Monocyte_FCGR3A" = '#a7993d',
"Megakaryocyte" = GREY,
"Hematopoietic_stem_cell"    = GREY
)


plotter$labels <- annotation[rownames(annotation) %in% colnames(all_data), 'labels']

pdf('./Plots/TT_MonFcVsMonCD_TSNE.pdf')
ggplot(plotter, aes(x=V1, y=V2, color = labels)) + theme_bw() +geom_point(size = 1.5) +scale_color_manual(values = colors)
dev.off()




# MonFC Vs MonFC assigned as MonCD



target <- 'Monocyte_FCGR3A'
obj    <- 'Monocyte_CD14'

G1 <- annotation[annotation$labels == target & annotation$prediction == target, 'cell_names']
G2 <- annotation[annotation$labels == target & annotation$prediction == obj, 'cell_names']

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
tmp   <- topTable(fit2_tmp, adjust="BH", n=Inf)
tmp$gene_name <- rownames(tmp)
tmp <- merge(tmp, gencode22, by= 'gene_name')



tmp[tmp$P.Value == 0, 'P.Value'] <- 1.445749e-281
tmp[tmp$adj.P.Val == 0, 'P.Value'] <- 1.445749e-281

pdf('./Plots/TT_MonFCVsMonFC_classMonCD_VP.pdf')
	graphContrast(tmp," (Ctrl B > 0, FC>0.5)", 0, 0.5, 1)
dev.off()

tmp$logpval <- -log(tmp$P.Value)
tmp2 <- tmp[ order(-tmp$logpval), c('gene_name', 'logFC', 'P.Value' ,'logpval')]

saveRDS(tmp2, './Data/TT_MonFCVsMonFC_classMonCD.rds')


data2heat <- data_tmp[rownames(data_tmp) %in%  tmp2[1:20, 'gene_name'],]
data2heat[data2heat > 5] <- 5

ann <- data.frame(Var1 = annotation[rownames(annotation) %in% colnames(data2heat), 'predictions'])
rownames(ann) <- colnames(data2heat)
Var1        <- c("#b577a8", "#979858")
names(Var1) <- c(levels(ann$Var1))
anno_colors <- list(Var1 = Var1)

pdf('./Plots/TT_MonFCVsMonFC_classMonCD_HM.pdf')
pheatmap( data2heat, cluster_rows = T, annotation = ann, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", annotation_colors = anno_colors, show_colnames = F)
dev.off()



# unique(annotation$labels)
GREY <- '#74797a'
TARGET <- '#cc5f43'
colors <- c(
"TARGET" = TARGET,
"NK_cell"  = GREY,
"Monocyte_CD14" = '#59a1bd',
"CD8_T_cell" = GREY,
"B_cell" = GREY,
"CD4_T_cell" = GREY,
"Plasmacytoid_dendritic_cell" = GREY,
"Monocyte_FCGR3A" = '#a7993d',
"Megakaryocyte" = GREY,
"Hematopoietic_stem_cell"    = GREY
)


ann_color <- annotation[rownames(annotation) %in% colnames(all_data), c('cell_names', 'labels')]
ann_color[ann_color$cell_names %in% G2, 'labels'] <- 'TARGET'


plotter$labels <- ann_color$labels

pdf('./Plots/TT_MonFCVsMonFC_classMonCD_TSNE.pdf')
ggplot(plotter, aes(x=V1, y=V2, color = labels)) + theme_bw() +geom_point(size = 1.5) +scale_color_manual(values = colors)
dev.off()




################################
# MonFC Vs MonCD assigned as MonFC

target <- 'Monocyte_FCGR3A'
obj    <- 'Monocyte_CD14'

G1 <- annotation[annotation$labels == target & annotation$prediction == target, 'cell_names']
G2 <- annotation[annotation$labels == obj & annotation$prediction == target, 'cell_names']
# Not enough samples

################################
################################

# MonFC Vs MonFC as Unassigned

target <- 'Monocyte_FCGR3A'
obj    <- 'Monocyte_CD14'

G1 <- annotation[annotation$labels == target & annotation$prediction == target, 'cell_names']
G2 <- annotation[annotation$labels == target & annotation$raw_predictions == target & annotation$prediction == 'Unassigned', 'cell_names']

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
tmp   <- topTable(fit2_tmp, adjust="BH", n=Inf)
tmp$gene_name <- rownames(tmp)
tmp <- merge(tmp, gencode22, by= 'gene_name')


tmp[tmp$P.Value == 0, 'P.Value'] <- 1.445749e-281
tmp[tmp$adj.P.Val == 0, 'P.Value'] <- 1.445749e-281


pdf('./Plots/TT_MonFCVsMonFC_Unassigned_VP.pdf')
	graphContrast(tmp," (Ctrl B > 0, FC>0.5)", 0, 0.5, 1)
dev.off()

tmp$logpval <- -log(tmp$P.Value)
tmp2 <- tmp[ order(-tmp$logpval), c('gene_name', 'logFC', 'P.Value' ,'logpval')]

saveRDS(tmp2, './Data/TT_MonFCVsMonFC_Unassigned.rds')


data2heat <- data_tmp[rownames(data_tmp) %in%  tmp2[1:20, 'gene_name'],]
data2heat[data2heat > 5] <- 5

ann <- data.frame(Var1 = annotation[rownames(annotation) %in% colnames(data2heat), 'predictions'])
rownames(ann) <- colnames(data2heat)
Var1        <- c("#b577a8", "#979858")
names(Var1) <- c(levels(ann$Var1))
anno_colors <- list(Var1 = Var1)

pdf('./Plots/TT_MonFCVsMonFC_Unassigned_HM.pdf')
pheatmap( data2heat, cluster_rows = T, annotation = ann, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", annotation_colors = anno_colors, show_colnames = F)
dev.off()



# unique(annotation$labels)

GREY <- '#74797a'
TARGET <- '#cc5f43'
colors <- c(
"TARGET" = TARGET,
"NK_cell"  = GREY,
"Monocyte_CD14" = '#59a1bd',
"CD8_T_cell" = GREY,
"B_cell" = GREY,
"CD4_T_cell" = GREY,
"Plasmacytoid_dendritic_cell" = GREY,
"Monocyte_FCGR3A" = '#a7993d',
"Megakaryocyte" = GREY,
"Hematopoietic_stem_cell"    = GREY
)


ann_color <- annotation[rownames(annotation) %in% colnames(all_data), c('cell_names', 'labels')]
ann_color[ann_color$cell_names %in% G2, 'labels'] <- 'TARGET'


plotter$labels <- ann_color$labels

pdf('./Plots/TT_MonFCVsMonFC_Unassigned_TSNE.pdf')
ggplot(plotter, aes(x=V1, y=V2, color = labels)) + theme_bw() +geom_point(size = 1.5) +scale_color_manual(values = colors)
dev.off()




# Check the confused cell samples

library(VennDiagram)
pd <- import("pandas")


annotation_all <- pd$read_pickle('/home/sevastopol/data/gserranos/JIND_DE/Data/JIND_assignmentbrftune_allLabels.pkl')
annotation_all$cell_names <- rownames(annotation_all)

annotation_cd4cd8 <- pd$read_pickle('/home/sevastopol/data/gserranos/JIND_DE/Data/JIND_assignmentbrftune_cd4_cd8.pkl')
annotation_cd4cd8$cell_names <- rownames(annotation_cd4cd8)

cd4 <- 'CD4 T cell'
cd8 <- 'CD8 T cell'

all    <- annotation_all[annotation_all$labels == cd4 & annotation_all$predictions != cd4, 'cell_names']
cd4cd8 <- annotation_cd4cd8[annotation_cd4cd8$labels == cd4 & annotation_cd4cd8$predictions != cd4, 'cell_names']


venn.diagram(list(all, cd4cd8), category.names=c('all', 'cd4cd8'), filename='/home/sevastopol/data/gserranos/JIND_DE/Plots/Test.png', imagetype='png')

all_missclass <- all[which(!(all %in% cd4cd8))]
annotation_all[annotation_all$cell_names %in% all_missclass, 'labels']




# pdf_files <- list.files('./Plots/', pattern = '^TT')

# convert TT_cd4realVscd4_classCd8_HM.pdf TT_cd4realVscd4_classCd8_HM.png
# convert TT_cd4realVscd4_classCd8_TSNE.pdf TT_cd4realVscd4_classCd8_TSNE.png
# convert TT_cd4realVscd4_classCd8_VP.pdf TT_cd4realVscd4_classCd8_VP.png
# convert TT_cd4realVscd4_Unassigned_HM.pdf TT_cd4realVscd4_Unassigned_HM.png
# convert TT_cd4realVscd4_Unassigned_TSNE.pdf TT_cd4realVscd4_Unassigned_TSNE.png
# convert TT_cd4realVscd4_Unassigned_VP.pdf TT_cd4realVscd4_Unassigned_VP.png
# convert TT_cd4realVscd8_classCd4_HM.pdf TT_cd4realVscd8_classCd4_HM.png
# convert TT_cd4realVscd8_classCd4_TSNE.pdf TT_cd4realVscd8_classCd4_TSNE.png
# convert TT_cd4realVscd8_classCd4_VP.pdf TT_cd4realVscd8_classCd4_VP.png
# convert TT_cd4Vscd8_HM.pdf TT_cd4Vscd8_HM.png
# convert TT_cd4Vscd8_TSNE.pdf TT_cd4Vscd8_TSNE.png
# convert TT_cd4Vscd8_VP.pdf TT_cd4Vscd8_VP.png
# convert TT_MonFcVsMonCD_HM.pdf TT_MonFcVsMonCD_HM.png
# convert TT_MonFcVsMonCD_TSNE.pdf TT_MonFcVsMonCD_TSNE.png
# convert TT_MonFcVsMonCD_VP.pdf TT_MonFcVsMonCD_VP.png
# convert TT_MonFCVsMonFC_classMonCD_HM.pdf TT_MonFCVsMonFC_classMonCD_HM.png
# convert TT_MonFCVsMonFC_classMonCD_TSNE.pdf TT_MonFCVsMonFC_classMonCD_TSNE.png
# convert TT_MonFCVsMonFC_classMonCD_VP.pdf TT_MonFCVsMonFC_classMonCD_VP.png
# convert TT_MonFCVsMonFC_Unassigned_HM.pdf TT_MonFCVsMonFC_Unassigned_HM.png
# convert TT_MonFCVsMonFC_Unassigned_TSNE.pdf TT_MonFCVsMonFC_Unassigned_TSNE.png
# convert TT_MonFCVsMonFC_Unassigned_VP.pdf TT_MonFCVsMonFC_Unassigned_VP.png



#######################################
#######################################
# FINAL FIGURE 
#######################################
#######################################

annotation$Final <- ifelse(annotation$label == annotation$predictions, annotation$label , 'Miss')

colors <- c(
"NK_cell"  = "#cb5698",
"Monocyte_CD14" = "#b07052",
"CD8_T_cell" = "#92cbb2",
"B_cell" = "#4d623c",
"CD4_T_cell" = "#4b2e48",
"Plasmacytoid_dendritic_cell" = "#ccbb51",
"Monocyte_FCGR3A" = '#8049be',
"Megakaryocyte" = "#79cd57",
"Hematopoietic_stem_cell"    = '#5a628c',
"Miss" = '#ea5148'
)



ann_color <- annotation[rownames(annotation) %in% colnames(all_data), c('cell_names', 'labels')]
ann_color[ann_color$cell_names %in% G2, 'labels'] <- 'TARGET'


plotter$labels <- annotation$Final

pdf('./Plots/Final_TSNE.pdf', width=10, height=10)
ggplot(plotter, aes(x=V1, y=V2, color = labels)) + theme_classic() +geom_point(size = 1.5, alpha = 0.7) +scale_color_manual(values = colors)
dev.off()




target <- 'Monocyte_FCGR3A'
obj    <- 'Monocyte_CD14'

G1 <- annotation[annotation$labels == target & annotation$prediction == target, 'cell_names']
G2 <- annotation[annotation$labels == target & annotation$prediction == obj, 'cell_names']

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
tmp   <- topTable(fit2_tmp, adjust="BH", n=Inf)
tmp$gene_name <- rownames(tmp)
tmp <- merge(tmp, gencode22, by= 'gene_name')



tmp[tmp$P.Value == 0, 'P.Value'] <- 1.445749e-281
tmp[tmp$adj.P.Val == 0, 'P.Value'] <- 1.445749e-281

pdf('./Plots/Final_MonFCVsMonFC_classMonCD_VP.pdf')
	graphContrast(tmp," (Ctrl B > 0, FC>0.5)", 0, 0.5, 1)
dev.off()

tmp$logpval <- -log(tmp$P.Value)
tmp2 <- tmp[ order(-tmp$logpval), c('gene_name', 'logFC', 'P.Value' ,'logpval')]

saveRDS(tmp2, './Data/TT_MonFCVsMonFC_classMonCD.rds')


data2heat <- data_tmp[rownames(data_tmp) %in%  tmp2[1:20, 'gene_name'],]
data2heat[data2heat > 5] <- 5

ann <- data.frame(Var1 = annotation[rownames(annotation) %in% colnames(data2heat), 'predictions'])
rownames(ann) <- colnames(data2heat)
Var1        <- c(colors['Monocyte_CD14'], colors['Monocyte_FCGR3A'])
names(Var1) <- c(levels(ann$Var1))
anno_colors <- list(Var1 = Var1)

pdf('./Plots/Final_MonFCVsMonFC_classMonCD_HM.pdf')
pheatmap( data2heat, cluster_rows = T, annotation_col = ann, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", cellheight= 10,annotation_colors = anno_colors, show_colnames = F)
dev.off()



# unique(annotation$labels)
GREY <- '#74797a'
TARGET <- '#cc5f43'
colors <- c(
"TARGET" = TARGET,
"NK_cell"  = GREY,
"Monocyte_CD14" = '#59a1bd',
"CD8_T_cell" = GREY,
"B_cell" = GREY,
"CD4_T_cell" = GREY,
"Plasmacytoid_dendritic_cell" = GREY,
"Monocyte_FCGR3A" = '#a7993d',
"Megakaryocyte" = GREY,
"Hematopoietic_stem_cell"    = GREY
)


ann_color <- annotation[rownames(annotation) %in% colnames(all_data), c('cell_names', 'labels')]
ann_color[ann_color$cell_names %in% G2, 'labels'] <- 'TARGET'


plotter$labels <- ann_color$labels

pdf('./Plots/Final_MonFCVsMonFC_classMonCD_TSNE.pdf')
ggplot(plotter, aes(x=V1, y=V2, color = labels)) + theme_bw() +geom_point(size = 1.5) +scale_color_manual(values = colors)
dev.off()

