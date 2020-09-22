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




#################################
#################################
# Load Data  ||  HUMAN BLOOD ||
#################################
#################################


pd <- import("pandas")
df <- pd$read_pickle('/home/sevastopol/data/gserranos/JIND_DE/Data/Pancreas_02/test.pkl')

annotation <- pd$read_pickle('/home/sevastopol/data/gserranos/JIND_DE/Data/Pancreas_02/JIND_assignmentbrftune.pkl')
annotation$cell_names <- rownames(annotation)

all_data <- t(df[, -which(colnames(df) %in% c('labels'))])
stand_dev <- apply(all_data,1, sd)
all_data <- all_data[!rownames(all_data) %in% names(which(stand_dev ==0)), ]

tsne_out <- Rtsne(as.matrix(t(all_data)))
plotter <-  as.data.frame(tsne_out$Y)

annotation$labels <- gsub(' ', '_', annotation$labels)
annotation$predictions <- gsub(' ', '_', annotation$predictions)
annotation$raw_predictions <- gsub(' ', '_', annotation$raw_predictions)





#######
# TSNE
#######





colors = c(
"alpha" = rgb(31, 119, 180, max=255), #BLUE
"ductal" = rgb(255, 127, 14, max=255), #ORANGE
"mast" = rgb(127, 127, 127, max=255), #GREY
"beta" =rgb(44, 160, 44, max=255), #GREEN
"gamma" =  rgb(214, 39, 40, max=255), #RED
"acinar" = rgb(148, 103, 189, max=255), #PURPLE
"delta"  = rgb(140, 86, 75, max=255), #BROWN
"epsilon"    = rgb(188, 189, 34, max=255), # LIME
"endothelial" = rgb(227, 119, 194, max=255)) # PINK



ann_color <- annotation[rownames(annotation) %in% colnames(all_data), c('cell_names', 'labels')]


names(plotter) <- c('tSNE1', 'tSNE2')
plotter$predictions <- annotation$labels
rownames(plotter) <- annotation$cell_names
plotter$shapes <- ifelse(annotation$predictions == 'Unassigned', 18, 20)

pdf('./Plots/Final_TSNE_Pancreas02.pdf', width=7, height=5)
ggplot(plotter, aes(x=tSNE1, y=tSNE2, color = predictions)) + 
geom_point(size = 2, alpha = 1) + 
scale_shape_identity() + 
scale_color_manual(values = colors) + 
theme_classic() + 
theme(  legend.position="top", 
		legend.title=element_blank(),
		axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
		axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), 
		plot.title = element_text(hjust = 0.5)) +
		guides(colour=guide_legend(nrow=2, override.aes = list(size=5)))
dev.off()




target <- 'ductal'
obj    <- 'acinar'

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

pdf('./Plots/Final_AcinarVsDuctal_classMonCD_VP.pdf')
	graphContrast(tmp," (Ctrl B > 0, FC>0.5)", 0, 0.5, 1)
dev.off()

tmp$logpval <- -log(tmp$P.Value)
tmp2 <- tmp[ order(-tmp$logpval), c('gene_name', 'logFC', 'P.Value' ,'logpval')]


data2heat <- data_tmp[rownames(data_tmp) %in%  tmp2[1:100, 'gene_name'],]
data2heat[data2heat > 5] <- 5

ann <- data.frame(Var1 = annotation[rownames(annotation) %in% colnames(data2heat), 'predictions'])
rownames(ann) <- colnames(data2heat)
Var1        <- c(colors['acinar'], colors['ductal'])
names(Var1) <- c(levels(ann$Var1))
anno_colors <- list(Var1 = Var1)

ann <- data.frame(Group = annotation[rownames(annotation) %in% colnames(data2heat), 'predictions'])
ann$Group <-  ifelse(ann$Group == 'ductal', 'G1', 'G2')
rownames(ann) <- colnames(data2heat)
Group         <- c(colors['ductal'], colors['acinar'])
# names(Group)  <- c(levels(ann$Group))
names(Group) <- c('G1', 'G2')
anno_colors   <- list(Group = Group)

pdf('./Plots/Final_Pancreas02_HM.pdf', heigh = 30)
pheatmap( data2heat, cluster_rows = T, treeheight_row = 0, annotation_col = ann, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", cellheight= 10,annotation_colors = anno_colors, show_colnames = F, main = 'Heatmap between ductal classified as ductal (G1)\n and ductal classified as acinar (G2)', fontsize = 8,fontsize_row=10)
dev.off()

# annotation2 <- annotation[, c('labels', 'predictions')]
# annotation3 <- melt(annotation2)
# annotation3 <- as.data.frame(table(annotation3))
# cols10 <- c("#578d66","#6b49c4","#65b648","#cc54c3","#b99c49","#64396c","#d55936","#7a90c6","#7d4a31","#ce5f80")

# pdf('./Test.pdf', heigh=15, width=15)
# ggplot(annotation3,
#        aes(y = Freq, axis1 = labels, axis2 = predictions)) +
#   geom_alluvium(aes(fill = predictions, color = predictions), width = 2/5, alpha = alpha, knot.pos = 0.4, , discern = TRUE) +
#   geom_stratum(width = 2/5, color = "grey", discern=TRUE) +
#   geom_text(stat = "stratum", discern = TRUE, aes(label = after_stat(stratum))) +
#   scale_x_continuous(breaks = 1:2, labels = c("Labels", "Predictions"))     +
#   scale_fill_manual(values  = cols10) +
#   scale_color_manual(values = cols10) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(size = 12, face = "bold") , legend.position = 'none')
# dev.off()
  


target <- 'ductal'
obj    <- 'acinar'

G1 <- annotation[annotation$labels == target & annotation$prediction == target, 'cell_names']
G2 <- annotation[annotation$labels == target & annotation$raw_prediction == obj, 'cell_names']

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

pdf('./Plots/Final_AcinarVsDuctal_classMonCD_VP.pdf')
	graphContrast(tmp," (Ctrl B > 0, FC>0.5)", 0, 0.5, 1)
dev.off()

tmp$logpval <- -log(tmp$P.Value)
tmp2 <- tmp[ order(-tmp$logpval), c('gene_name', 'logFC', 'P.Value' ,'logpval')]


data2heat <- data_tmp[rownames(data_tmp) %in%  tmp2[1:100, 'gene_name'],]
data2heat[data2heat > 5] <- 5

ann <- data.frame(Var1 = annotation[rownames(annotation) %in% colnames(data2heat), 'predictions'])
rownames(ann) <- colnames(data2heat)
Var1        <- c(colors['acinar'], colors['ductal'])
names(Var1) <- c(levels(ann$Var1))
anno_colors <- list(Var1 = Var1)

ann <- data.frame(Group = annotation[rownames(annotation) %in% colnames(data2heat), 'predictions'])
ann$Group <-  ifelse(ann$Group == 'ductal', 'G1', 'G2')
rownames(ann) <- colnames(data2heat)
Group         <- c(colors['ductal'], colors['acinar'])
# names(Group)  <- c(levels(ann$Group))
names(Group) <- c('G1', 'G2')
anno_colors   <- list(Group = Group)

pdf('./Plots/Final_Pancreas02RAW_HM.pdf', heigh = 30)
pheatmap( data2heat, cluster_rows = T, treeheight_row = 0, annotation_col = ann, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", cellheight= 10,annotation_colors = anno_colors, show_colnames = F, main = 'Heatmap between ductal classified as ductal (G1)\n and ductal classified as acinar (G2)', fontsize = 8,fontsize_row=10)
dev.off()
