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



#### QC ####


pd <- import("pandas")
df <- pd$read_pickle('test.pkl')

df2 <- pd$read_pickle('JIND_assignmentbrftune.pkl')
cd8_cd4_samples <- df2[df2$labels %in% c('CD8 T cell','CD4 T cell'),]

# Working with everithing

data <- t(df[, -which(colnames(df) %in% c('labels'))])
stand_dev <- apply(data,1, sd)
data <- data[!rownames(data) %in% names(which(stand_dev ==0)), ]


pdf(file="./Plots/QC_JIND_ALL_Aug20_Bxplt.pdf", width = 40, height = 20)
boxplot(data[,1:100], names = colnames(data)[1:100],cex.axis=0.7, las=2)
dev.off()

xlim_r = c(min(na.omit(data)),max(na.omit(data)))

pdf(file="./Plots/QC_JIND_ALL_Aug20_Multi.pdf", width = 40, height = 20)
multi("density", data, xlim_r, "Density","","")
dev.off()

clust.euclid.average <- hclust(dist(t(data)),method="average")
pdf("./Plots/QC_JIND_ALL_Aug20_clust.pdf", width = 40, height = 20)
plot(clust.euclid.average, main="Hierarchical clustering of normalized samples",  hang=-1)
dev.off()



x="PC1"
y="PC2"
fit <- prcomp(t(na.omit(as.matrix(data))), scale=T)
pcData <- data.frame(fit$x)
ggp <- ggplot(pcData, aes_string(x=x, y=y)) + geom_point(aes(colour=df2$labels, shape=df2$predictions), size = 5) + theme_bw() + theme(legend.title=element_blank())+ scale_shape_manual(values=seq(0,15))
pdf("./Plots/PCA_JIND_ALL_Aug20.pdf", width = 15, height = 15, colormodel = "rgb")
ggp
dev.off()

###########################
# DE for all data -- 1 vs 1
###########################

data <- t(df[, -which(colnames(df) %in% c('labels'))])
stand_dev <- apply(data,1, sd)
data <- data[!rownames(data) %in% names(which(stand_dev ==0)), ]

contrasts_results1Vs1 <- list()
DESIGN <- data.frame(cells = rownames(df), label = gsub(' ', '_',df$labels))

for (target in levels(DESIGN$label)){
	for (obj in levels(DESIGN$label)){
		print(paste0(target, '_Vs_', obj))
		if (target == obj){next}
		cells_2_keep <- DESIGN[DESIGN$label %in% c(target, obj), 'cells']
		data_tmp <- data[, colnames(data) %in% cells_2_keep]

		design_tmp <- as.matrix(DESIGN[DESIGN$cells %in% cells_2_keep, 'label'])
    	design_tmp[design_tmp != target] <-1
    	design_tmp[design_tmp == target] <-0
    	design_tmp <- model.matrix(~0+as.factor(design_tmp[,1]))
    	colnames(design_tmp) <- c(target, obj)
    	rownames(design_tmp) <- DESIGN[DESIGN$cells %in% cells_2_keep, 'cells']


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
    	contrasts_results1Vs1[[paste0(target, '_vs_', obj)]] <- tmp
	}
}

saveRDS(contrasts_results1Vs1, '/home/sevastopol/data/gserranos/JIND/CD4_CD8_1Vs1.rds')



pdf('./Plots/1Vs1_B10_FC_1.pdf')
for (contrast in names(contrasts_results1Vs1)){
	print(contrast)
	graphContrast(contrasts_results1Vs1[[contrast]], paste0( contrast ," (Ctrl B > 10, FC>1)"), 10, 1, 1)
}
dev.off()



##########################
# Cd4 real Vs cd4 labelled
##########################


pd <- import("pandas")
df <- pd$read_pickle('test.pkl')

df2 <- pd$read_pickle('JIND_assignmentbrftune.pkl')
cd8_cd4_samples <- df2[df2$labels %in% c('CD8 T cell','CD4 T cell'),]

# Working with everithing

data <- t(df[, -which(colnames(df) %in% c('labels'))])
stand_dev <- apply(data,1, sd)
data <- data[!rownames(data) %in% names(which(stand_dev ==0)), ]

data <- t(df[, -which(colnames(df) %in% c('labels'))])
stand_dev <- apply(data,1, sd)
data <- data[!rownames(data) %in% names(which(stand_dev ==0)), ]

gencode22 <- read.table('/home/sevastopol/data/gserranos/UTILS/gencode.v22.annotation', header=T)
gencode22$gene_id <- str_extract(gencode22$gene_id, 'ENSG[R]?[0-9]+')

cd8_cd4_samples$cell_names <- rownames(cd8_cd4_samples)

true_cd4 <- cd8_cd4_samples[cd8_cd4_samples$labels == 'CD4 T cell' & cd8_cd4_samples$predictions == 'CD4 T cell', 'cell_names']
false_cd4 <- cd8_cd4_samples[ cd8_cd4_samples$labels == 'CD4 T cell' & cd8_cd4_samples$predictions == 'CD8 T cell','cell_names']

selection <- c(true_cd4, false_cd4)
data_cd4_labelVspred <- data[ , colnames(data) %in%  selection]


DESIGN <- data.frame(cells = true_cd4, labels = 'CD4')
DESIGN <- rbind(DESIGN, data.frame(cells= false_cd4, labels = 'CD8'))

DESIGN <- DESIGN[ order(colnames(data_cd4_labelVspred)), ]


design_tmp <- as.matrix(DESIGN[, 'labels'])
design_tmp[design_tmp != 'CD4'] <-1
design_tmp[design_tmp == 'CD8'] <-0
design_tmp <- model.matrix(~0+as.factor(design_tmp[,1]))
colnames(design_tmp) <- c('CD4', 'CD8')
rownames(design_tmp) <- colnames(data_cd4_labelVspred)


### DE
# create the linear model
fit_tmp <- lmFit(data_cd4_labelVspred, design_tmp)
# model correction
fit_tmp <- eBayes(fit_tmp)
# results <- topTable(fit_tmp, n=Inf)
x <- paste0('CD4', '-', 'CD8')
contrast_mat_tmp <- makeContrasts(contrasts=x, levels= c('CD4', 'CD8'))
fit2_tmp <- contrasts.fit(fit_tmp, contrast_mat_tmp)
fit2_tmp <- eBayes(fit2_tmp)
tmp   <- topTable(fit2_tmp, adjust="BH", n=Inf)
tmp$gene_name <- rownames(tmp)
tmp <- merge(tmp, gencode22, by= 'gene_name')



saveRDS(tmp, '/home/sevastopol/data/gserranos/JIND/CD4_realVsPred.rds')



pdf('./Plots/CD4_realVsPred_B10_FC_1.pdf')
	graphContrast(tmp, paste0( 'CD4_realVsPred' ," (Ctrl B > 10, FC>1)"), 10, 1, 1)
dev.off()


tmp$logpval <- -log(tmp$P.Value)
tmp2 <- tmp[ order(-tmp$logpval), c('gene_name', 'logFC', 'logpval')]

# gene_name      logFC  logpval
#    LCN8 0.12159118 21.64040
#   PTAFR 0.04547547 21.63423
#  MARCKS 0.02390872 21.61546
# FAM129C 0.02159865 21.60963
#   PILRA 0.06759608 19.78322


oneVsone <- readRDS('/home/sevastopol/data/gserranos/JIND/CD4_CD8_1Vs1.rds')
CD4vsFalseCD4 <- readRDS('/home/sevastopol/data/gserranos/JIND/CD4_realVsPred.rds')

cd4Vscd8 <- oneVsone[['CD4_T_cell_vs_CD8_T_cell']]
cd4Vscd8$logpval <- -log(cd4Vscd8$P.Value)
genes <- c('LCN8','PTAFR','MARCKS','FAM129C','PILRA')

cd4Vscd8[cd4Vscd8$gene_name %in% genes, c('gene_name', 'logFC','logpval')]

#     gene_name         logFC    logpval
# 287   FAM129C -8.581239e-03 4.88012333
# 606      LCN8 -1.876204e-04 0.03018702
# 664    MARCKS -1.254598e-03 0.88690657
# 777     PILRA  4.196452e-05 0.01105023
# 813     PTAFR  1.422263e-03 0.94833303




# CD4 Vs CD4_unassigned

gencode22 <- read.table('/home/sevastopol/data/gserranos/UTILS/gencode.v22.annotation', header=T)
gencode22$gene_id <- str_extract(gencode22$gene_id, 'ENSG[R]?[0-9]+')


pd <- import("pandas")
df <- pd$read_pickle('test.pkl')

df2 <- pd$read_pickle('JIND_assignmentbrftune.pkl')
cd8_cd4_samples <- df2[df2$labels %in% c('CD8 T cell','CD4 T cell'),]

# Working with everithing

data <- t(df[, -which(colnames(df) %in% c('labels'))])
stand_dev <- apply(data,1, sd)
data <- data[!rownames(data) %in% names(which(stand_dev ==0)), ]


stats <- as.data.frame(table(df2[df2$predictions == 'Unassigned','labels']))

df2$cell_names <- rownames(df2)
true_cd4 <- df2[df2$labels == 'CD4 T cell' & df2$predictions == 'CD4 T cell', 'cell_names']
cd4_unassigned <- df2[ df2$labels == 'CD4 T cell' & df2$predictions == 'Unassigned' & df2$raw_predictions == 'CD4 T cell','cell_names']

selection <- c(true_cd4, cd4_unassigned)
data_cd4VsUn <- data[ , colnames(data) %in%  selection]


DESIGN <- data.frame(cells = true_cd4, labels = 'CD4')
DESIGN <- rbind(DESIGN, data.frame(cells= cd4_unassigned, labels = 'Unassigned'))

DESIGN <- DESIGN[ order(colnames(data_cd4VsUn)), ]


design_tmp <- as.matrix(DESIGN[, 'labels'])
design_tmp[design_tmp != 'CD4'] <-1
design_tmp[design_tmp == 'Unassigned'] <-0
design_tmp <- model.matrix(~0+as.factor(design_tmp[,1]))
colnames(design_tmp) <- c('CD4', 'CD8')
rownames(design_tmp) <- colnames(data_cd4VsUn)


### DE
# create the linear model
fit_tmp <- lmFit(data_cd4VsUn, design_tmp)
# model correction
fit_tmp <- eBayes(fit_tmp)
# results <- topTable(fit_tmp, n=Inf)
x <- paste0('CD4', '-', 'Unassigned')
contrast_mat_tmp <- makeContrasts(contrasts=x, levels= c('CD4', 'Unassigned'))
fit2_tmp <- contrasts.fit(fit_tmp, contrast_mat_tmp)
fit2_tmp <- eBayes(fit2_tmp)
tmp   <- topTable(fit2_tmp, adjust="BH", n=Inf)
tmp$gene_name <- rownames(tmp)
tmp <- merge(tmp, gencode22, by= 'gene_name')

saveRDS(tmp, 'CD4VsUnassignedCD4.rds')

pdf('./Plots/CD4_VsUnassigned_B10_FC_1.pdf')
	graphContrast(tmp, paste0( 'CD4_VsUnassigned' ," (Ctrl B > 10, FC>1)"), 10, 1, 1)
dev.off()




# Mohit Live



pd <- import("pandas")
df <- pd$read_pickle('test.pkl')

df2 <- pd$read_pickle('JIND_assignmentbrftune.pkl')


data <- t(df[, -which(colnames(df) %in% c('labels'))])
stand_dev <- apply(data,1, sd)
data <- data[!rownames(data) %in% names(which(stand_dev ==0)), ]




df2$cell_names <- rownames(df2)
G4 <- df2[df2$labels == 'CD8 T cell' & df2$predictions == 'Unassigned' ,'cell_names']
G3 <- df2[df2$labels == 'CD4 T cell' & df2$predictions == 'Unassigned' ,'cell_names']


G2 <- df2[df2$labels == 'CD8 T cell' & df2$raw_predictions == 'CD4 T cell', 'cell_names']
real_cd8 <- df2[df2$labels == 'CD8 T cell' , 'cell_names']

real_cd8 <-  real_cd8[-which(real_cd8 %in% G2)]



selection <- c(G2, real_cd8)
data_CD8VsG2 <- data[ , colnames(data) %in%  selection]


DESIGN <- data.frame(cells = G2, labels = 'G2')
DESIGN <- rbind(DESIGN, data.frame(cells= real_cd8, labels = 'CD8'))

DESIGN <- DESIGN[ order(colnames(data_CD8VsG2)), ]


design_tmp <- as.matrix(DESIGN[, 'labels'])
design_tmp[design_tmp == 'G2'] <-1
design_tmp[design_tmp == 'CD8'] <-0
design_tmp <- model.matrix(~0+as.factor(design_tmp[,1]))
colnames(design_tmp) <- c('CD8', 'G2')
rownames(design_tmp) <- colnames(data_CD8VsG2)



### DE
# create the linear model
fit_tmp <- lmFit(data_CD8VsG2, design_tmp)
# model correction
fit_tmp <- eBayes(fit_tmp)
# results <- topTable(fit_tmp, n=Inf)
x <- paste0('CD8', '-', 'G2')
contrast_mat_tmp <- makeContrasts(contrasts=x, levels= c('CD8', 'G2'))
fit2_tmp <- contrasts.fit(fit_tmp, contrast_mat_tmp)
fit2_tmp <- eBayes(fit2_tmp)
tmp   <- topTable(fit2_tmp, adjust="BH", n=Inf)
tmp$gene_name <- rownames(tmp)


pdf('./Plots/CD8_VsG2_B10_FC_1.pdf')
	graphContrast(tmp, paste0( 'CD4_VsUnassigned' ," (Ctrl B > 5, FC>1)"), 5, 1, 1)
dev.off()