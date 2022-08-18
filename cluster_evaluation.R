latent <- read.table("AML_data/AML707B/all/AML707B-all_clusters.csv", sep=",", header=TRUE)
predict.labels <- read.table("AML_data/AML707B/all/AML707B-all_predict.csv", sep=",", header=TRUE)
truth.labels <- read.table("AML_data/AML707B/all/AML707B-all_celltype_merge.csv", sep="\t", header=TRUE)
truth.labels <- as.data.frame(t(truth.labels))
truth.labels <- truth.labels[2:nrow(truth.labels),]

latent$predict <- as.character(predict.labels$X0)
latent$truth <- as.character(truth.labels)
latent <- cbind(latent, expression)

latent$truth[latent$truth == "TRUE"] <- "T"

latent.sample <- latent[latent$truth == "B" | latent$truth == "lateEry" | latent$truth == "T" | latent$truth == "GMP" | latent$truth == "Mono" | latent$truth == "NK" | latent$truth == "Plasma",]

latent.sample <- latent.sample[order(latent.sample$new_labels),]
dist.mat <- dist(latent.sample[,1:2])
#dist.mat <- data.frame(dist.mat)

library(ggplot2)
library(reshape2)
library(dplyr)

# heatmap
dist.plot <- as.matrix(dist.mat)
dist.plot <- as.data.frame(dist.plot)
#dist.plot <- dist.mat
colnames(dist.plot) <- c(1:ncol(dist.plot))
rownames(dist.plot) <- c(1:nrow(dist.plot))
dist.plot <- dist.plot %>% mutate(B=colnames(dist.plot)) %>% melt()
dist.plot$B <- as.integer(dist.plot$B)
dist.plot$variable <- as.integer(dist.plot$variable)
p1 <- ggplot(dist.plot,aes(x=variable,y=B,fill=value)) 
p2 <- p1 + geom_raster() + scale_fill_gradient2(low="yellow", high="blue", mid="white")
p2
p2 + xlab("Cells") + ylab("Cells")

#latent.sample <- latent.sample[order(latent.sample$new_labels),]
#latent.sample.select <- seq(1, nrow(latent.sample), by=10)
#latent.sample <- latent.sample[latent.sample.select,]

clusters <- latent.sample$new_labels
dist.mat <- dist(latent.sample[,1:2])
dist.plot <- as.matrix(dist.mat)
dist.plot <- as.data.frame(dist.plot)
colnames(dist.plot) <- c(1:ncol(dist.plot))
rownames(dist.plot) <- c(1:nrow(dist.plot))
dist.plot <- dist.plot %>% mutate(B=colnames(dist.plot)) %>% melt()
dist.plot$B <- as.integer(dist.plot$B)
dist.plot$variable <- as.integer(dist.plot$variable)
#rownames(clusters) <- c(1:nrow(clusters))

# histogram
ggplot(data = dist.plot, mapping = aes(x = value)) + geom_histogram(binwidth = 0.5)

dist.plot$celltype = c(1:nrow(dist.plot))
dist.plot$celltype[dist.plot$B %in% which(clusters == "Erythrocytes") & dist.plot$variable %in% which(clusters == "Erythrocytes")] <- "Erythrocytes"
dist.plot$celltype[dist.plot$B %in% which(clusters == "Erythrocytes") & dist.plot$variable %in% which(clusters == "GMPs")] <- "GMPs"
dist.plot$celltype[dist.plot$B %in% which(clusters == "Erythrocytes") & dist.plot$variable %in% which(clusters == "Monocytes")] <- "Monocytes"
dist.plot$celltype[dist.plot$B %in% which(clusters == "Erythrocytes") & dist.plot$variable %in% which(clusters == "NK cells")] <- "NK cells"
dist.plot$celltype[dist.plot$B %in% which(clusters == "Erythrocytes") & dist.plot$variable %in% which(clusters == "Plasma and B cells")] <- "Plasma and B cells"
dist.plot$celltype[dist.plot$B %in% which(clusters == "Erythrocytes") & dist.plot$variable %in% which(clusters == "T cells")] <- "T cells"

dist.plot <- dist.plot[dist.plot$B %in% which(clusters == "Erythrocytes"),]

p <- ggplot(data = dist.plot, mapping = aes(x = value, y=..density.., color = celltype)) + geom_freqpoly(binwidth = 2.5)
p <- p + guides(color=guide_legend(title = "Cells")) + theme_bw()
p <- p + theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + labs(x = "Distance", y = "Density")
p
