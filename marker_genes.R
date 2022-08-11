expression <- latent.sample[,1:4]
expression$new_labels <- latent.sample$new_labels

#DEG <- read.table("AML_data/AML707B/all/scCODE_707B_results.csv", sep=",", header=TRUE)

#select.genes <- DEG$DE_results.Gene_name
#select.genes <- latent.sample[,colnames(latent.sample) %in% select.genes]
#select.genes <- scale(select.genes)
select.genes <- latent.sample[5:(ncol(latent.sample)-1)]

#select.genes <- latent.sample[,6:(ncol(latent.sample)-1)]
select.genes <- scale(select.genes)
#select.genes <- as.data.frame(select.genes)

expression <- cbind(expression, select.genes)

expression <- expression[order(expression$new_labels),]
select.genes <- expression[6:ncol(expression)]

# see expression of some marker genes
ggplot(data = expression) + geom_point(mapping = aes(x = X0, y = X1, color = HBB), size = 1)
ggplot(data = expression) + geom_point(mapping = aes(x = X0, y = X1, color = MS4A1), size = 1)
ggplot(data = expression) + geom_point(mapping = aes(x = X0, y = X1, color = CCR7), size = 1)

mean.ery <- apply(expression[which(expression$new_labels == "Erythrocytes"), 6:ncol(expression)], 2, mean)
mean.gmp <- apply(expression[which(expression$new_labels == "GMPs"), 6:ncol(expression)], 2, mean)
mean.mono <- apply(expression[which(expression$new_labels == "Monocytes"), 6:ncol(expression)], 2, mean)
mean.nk <- apply(expression[which(expression$new_labels == "NK cells"), 6:ncol(expression)], 2, mean)
mean.b <- apply(expression[which(expression$new_labels == "Plasma and B cells"), 6:ncol(expression)], 2, mean)
mean.t <- apply(expression[which(expression$new_labels == "T cells"), 6:ncol(expression)], 2, mean)

DEG.matrix <- as.data.frame(matrix(data = 0, nrow = 6, ncol = ncol(expression)-5))
DEG.matrix[1,] <- mean.ery
DEG.matrix[2,] <- mean.gmp
DEG.matrix[3,] <- mean.mono
DEG.matrix[4,] <- mean.nk
DEG.matrix[5,] <- mean.b
DEG.matrix[6,] <- mean.t
colnames(DEG.matrix) <- colnames(expression)[6:ncol(expression)]

rownames(DEG.matrix) <- c("Erythrocytes", "GMPs", "Monocytes", "NK cells", "Plasma and B cells", "T cells")

DEG.matrix.sorted <- DEG.matrix[,order(DEG.matrix[1,], decreasing = TRUE)]
marker.ery <- colnames(DEG.matrix.sorted)[1:20]

DEG.matrix.sorted <- DEG.matrix[,order(DEG.matrix[2,], decreasing = TRUE)]
marker.gmp <- colnames(DEG.matrix.sorted)[1:20]
marker.gmp <- c(marker.gmp, "CD34")

DEG.matrix.sorted <- DEG.matrix[,order(DEG.matrix[3,], decreasing = TRUE)]
marker.mono <- colnames(DEG.matrix.sorted)[1:20]

DEG.matrix.sorted <- DEG.matrix[,order(DEG.matrix[4,], decreasing = TRUE)]
marker.nk <- colnames(DEG.matrix.sorted)[20:40]
marker.nk <- c(marker.nk, "NKG7")

DEG.matrix.sorted <- DEG.matrix[,order(DEG.matrix[5,], decreasing = TRUE)]
marker.b <- colnames(DEG.matrix.sorted)[145:160]
marker.b <- c("MS4A1", marker.b, "FCRLA")

DEG.matrix.sorted <- DEG.matrix[,order(DEG.matrix[6,], decreasing = TRUE)]
marker.t <- colnames(DEG.matrix.sorted)[1:20]

marker.all <- c(marker.ery, marker.gmp, marker.mono, marker.nk, marker.b, marker.t)
marker.all <- unique(marker.all)

select.DEG <- DEG.matrix[,marker.all]

#DEG <- DEG[order(DEG$DE_results.P_adjust),]
#gene.order <- DEG[DEG$sig == "up", "DE_results.Gene_name"]

library(dplyr)
library(reshape2)
#rownames(select.genes) <- c(1:nrow(select.genes))
select.DEG <- select.DEG %>% mutate(B=rownames(select.DEG)) %>% melt()
colnames(select.DEG)[3] <- "z_score"
p <- ggplot(select.DEG,aes(x=B,y=variable,fill=z_score))
p <- p + guides(color=guide_legend(title = "Scaled expression")) + theme_bw()
p <- p + theme(axis.text=element_text(size=6.5), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + geom_raster() + scale_fill_gradient2(low="blue", high="red", mid="white")
p <- p + xlab(NULL) + ylab(NULL)
p

marker.all.reverse <- rev(marker.all)
select.DEG.reverse <- DEG.matrix[,marker.all.reverse]
select.DEG.reverse <- select.DEG.reverse %>% mutate(B=rownames(select.DEG.reverse)) %>% melt()
colnames(select.DEG.reverse)[3] <- "z_score"
used.marker <- c("HBB", "AHSP", "CA2", "HBD", "CLEC11A", "PRSS57", "EGFL7", "SOX4", "CD34", "LYZ", "CD14", "GNLY", "GZMB", "GZMA", "NKG7", "MS4A1", "FCRLA", "FCRL5", "MZB1", "IL7R", "CD3D", "CCR7")
marker.all.reverse[!marker.all.reverse %in% used.marker] <- " "
p <- ggplot(select.DEG.reverse,aes(x=B,y=variable,fill=z_score))
p <- p + guides(color=guide_legend(title = "Scaled expression")) + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text=element_text(size=6.5), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + geom_raster() + scale_fill_gradient2(low="blue", high="red", mid="white")
p <- p + xlab(NULL) + ylab(NULL)
p
p <- ggplot(select.DEG.reverse,aes(x=B,y=variable,fill=z_score))
p <- p + guides(color=guide_legend(title = "Scaled expression")) + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text=element_text(size=12), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + geom_raster() + scale_fill_gradient2(low="blue", high="red", mid="white")
p <- p + scale_y_discrete(labels = marker.all.reverse)
p <- p + xlab(NULL) + ylab(NULL)
p


p <- ggplot(data = expression, mapping = aes(x = new_labels, y = CCR7, fill = new_labels)) + geom_violin()
p <- p + guides(color=guide_legend(title = "Cells")) + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text=element_text(size=15), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + labs(x = NULL)
p
