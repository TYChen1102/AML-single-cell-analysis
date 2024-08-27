library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)
library(Rtsne)
library(uwot)

d <- read.table("AML_res/BM4/BM4_RNN_after.csv", sep=",", header=TRUE)
id <- read.table("AML_res/BM4/BM4_predict.csv", sep=",", header=TRUE)
nclusters <- max(id$X0) + 1
id$X0 <- as.character(id$X0)

# t-SNE plot

tsne <- Rtsne(d, pca_scale = FALSE)
df.tsne <- data.frame(tsne$Y)
colnames(df.tsne) <- c("tSNE1","tSNE2")
rownames(df.tsne) <- rownames(d)
df.tsne <- cbind(df.tsne, id)
colnames(df.tsne)[ncol(df.tsne)] <- "id"

ggplot(df.tsne, aes(x=tSNE1,y=tSNE2, color=id)) + geom_point()

p <- ggplot(data = df.tsne) + geom_point(mapping = aes(x = tSNE1,y = tSNE2, color = id), size = 1)
p <- p + guides(color=guide_legend(title = "Clusters")) + theme_bw()
p <- p + theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + labs(x = "t-SNE1", y = "t-SNE2")
p

# ground truth

raw <- read.table("AML_data/BM4/BM4_data.csv", sep="\t", header=TRUE)
cells <- read.table("AML_data/BM4/BM4_celltype.csv", sep="\t", header=TRUE)
cells <- cells[2:ncol(cells)]
cells <- t(cells)
raw <- cbind(raw, cells)
colnames(raw)[ncol(raw)] <- "celltype"

count1 <- raw[cells == "lateEry",]
count2 <- raw[cells == "TRUE",]
count3 <- raw[cells == "B",]
count4 <- raw[cells == "Mono",]
count5 <- raw[cells == "NK",]
count6 <- raw[cells == "Plasma",]
count7 <- raw[cells == "GMP",]
count8 <- raw[cells == "Prog",]
raw <- rbind(count1,count2,count3,count4,count5,count6,count7,count8)

cells1 <- cells[cells == "lateEry"]
cells2 <- cells[cells == "TRUE"]
cells3 <- cells[cells == "B"]
cells4 <- cells[cells == "Mono"]
cells5 <- cells[cells == "NK"]
cells6 <- cells[cells == "Plasma"]
cells7 <- cells[cells == "GMP"]
cells8 <- cells[cells == "Prog"]
cells <- c(cells1,cells2,cells3,cells4,cells5,cells6,cells7,cells8)

se <- seq(1,nrow(raw),2)
raw <- raw[se,]
cells <- cells[se]
d <- cbind(raw, df.tsne)

ggplot(d, aes(x=tSNE1,y=tSNE2, color=celltype)) + geom_point()

d$celltype[d$celltype==TRUE] <- "T"

p <- ggplot(data = d) + geom_point(mapping = aes(x = tSNE1,y = tSNE2, color = celltype), size = 1)
p <- p + guides(color=guide_legend(title = "Clusters")) + theme_bw()
p <- p + theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + labs(x = "t-SNE1", y = "t-SNE2")
p

# DEG

d[,1:(ncol(d)-4)] <- scale(d[,1:(ncol(d)-4)])
exp.mean <- matrix(data = 0, nrow = nclusters, ncol = ncol(d)-4)
exp.mean <- as.data.frame(exp.mean)
colnames(exp.mean) <- colnames(d)[2:(ncol(d)-4)]
for (i in 1:nclusters){
  m <- apply(d[d$id == as.character(i-1), 1:(ncol(d)-4)], 2, mean)
  exp.mean[i,] <- m
}

genelist <- c()
for (i in 1:nclusters){
  m <- exp.mean[i,]
  x <- as.numeric(m[1,])
  x <- order(x, decreasing = T)
  x <- x[1:50]
  g <- colnames(m)[x]
  genelist <- c(genelist, g)
}
genelist <- unique(genelist)

d.plot <- exp.mean[,genelist]
d.plot <- d.plot %>% mutate(B=rownames(d.plot)) %>% melt()
colnames(d.plot)[3] <- "z_score"
p <- ggplot(d.plot,aes(x=B,y=variable,fill=z_score))
p <- p + guides(color=guide_legend(title = "Scaled expression")) + theme_bw()
p <- p + theme(axis.text=element_text(size=6.5), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + geom_raster() + scale_fill_gradient2(low="blue", high="red", mid="white")
p <- p + xlab(NULL) + ylab(NULL)
p

c("HBB", "AHSP", "CA2", "HBD") %in% genelist
which(genelist == "HBB")
c("CLEC11A", "PRSS57", "EGFL7", "SOX4", "CD34") %in% genelist
c("LYZ", "CD14") %in% genelist
c("GNLY", "GZMB", "GZMA", "NKG7") %in% genelist
c("MS4A1", "FCRLA", "FCRL5", "MZB1") %in% genelist
c("IL7R", "CD3D", "CCR7") %in% genelist

used.marker <- c("HBB", "AHSP", "CA2", "HBD", "CLEC11A", "PRSS57", "EGFL7", "SOX4", "CD34", "LYZ", "CD14", "GNLY", "GZMB", "GZMA", "NKG7", "MS4A1", "FCRLA", "FCRL5", "MZB1", "IL7R", "CD3D", "CCR7")
marker <- genelist
marker[!marker %in% used.marker] <- " "

p <- ggplot(d.plot,aes(x=B,y=variable,fill=z_score))
p <- p + guides(color=guide_legend(title = "Scaled expression")) + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text=element_text(size=12), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + geom_raster() + scale_fill_gradient2(low="blue", high="red", mid="white")
p <- p + scale_y_discrete(labels = marker)
p <- p + xlab(NULL) + ylab(NULL)
p

# annotation

d$pre <- d$id
d$pre[d$id=="0"] <- "T"
d$pre[d$id=="1"] <- "Mono"
d$pre[d$id=="2"] <- "Prog & GMP"
d$pre[d$id=="3"] <- "Prog & GMP"
d$pre[d$id=="4"] <- "Prog & GMP"
d$pre[d$id=="5"] <- "NK"
d$pre[d$id=="6"] <- "Ery"
d$pre[d$id=="7"] <- "Plasma & B"
d$pre[d$id=="8"] <- "Plasma & B"

ggplot(d, aes(x=tSNE1,y=tSNE2, color=pre)) + geom_point()

p <- ggplot(data = d) + geom_point(mapping = aes(x = tSNE1,y = tSNE2, color = pre), size = 1)
p <- p + guides(color=guide_legend(title = "Cells")) + theme_bw()
p <- p + theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + labs(x = "t-SNE1", y = "t-SNE2")
p

exp.mean.new <- matrix(data = 0, nrow = 6, ncol = ncol(d)-5)
exp.mean.new <- as.data.frame(exp.mean.new)
celltype <- unique(d$pre)
colnames(exp.mean.new) <- colnames(d)[1:(ncol(d)-5)]
for (i in 1:6){
  m <- apply(d[d$pre == celltype[i], 1:(ncol(d)-5)], 2, mean)
  exp.mean.new[i,] <- m
}
rownames(exp.mean.new) <- celltype

genelist.new <- c()
for (i in 1:6){
  m <- exp.mean.new[i,]
  x <- as.numeric(m[1,])
  x <- order(x, decreasing = T)
  x <- x[1:50]
  g <- colnames(m)[x]
  genelist.new <- c(genelist.new, g)
}
genelist.new <- unique(genelist.new)

d.plot <- exp.mean.new[,genelist.new]
d.plot <- d.plot %>% mutate(B=rownames(d.plot)) %>% melt()
colnames(d.plot)[3] <- "z_score"
d.plot <-within(d.plot,{B <-factor(B,levels=c("Prog & GMP","Ery","NK","T","Plasma & B","Mono"))})

used.marker <- c("HBB", "CA2", "HBD", "CLEC11A", "PRSS57", "EGFL7", "SOX4", "CD34", "LYZ", "CD14", "GNLY", "GZMB", "NKG7", "MS4A1", "IL7R", "CCR7")
marker <- genelist.new
marker[!marker %in% used.marker] <- " "

p <- ggplot(d.plot,aes(x=B,y=variable,fill=z_score))
p <- p + guides(color=guide_legend(title = "Scaled expression")) + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 10), axis.text=element_text(size=12), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + geom_raster() + scale_fill_gradient2(low="blue", high="red", mid="white")
p <- p + scale_y_discrete(labels = marker)
p <- p + xlab(NULL) + ylab(NULL)
p

# cf matrix

d$celltype[d$celltype == "TRUE"] <- "T"
d <-within(d,{pre <-factor(pre,levels=c("Prog & GMP","Ery","NK","T","Plasma & B","Mono"))})
d <-within(d,{celltype <-factor(celltype,levels=c("Prog","GMP","lateEry","NK","T","Plasma","B","Mono"))})

p <-ggplot(data = d) + geom_count(mapping = aes(x = pre, y = celltype))
p <- p + guides(color=guide_legend(title = "Count")) + theme_bw()
p <- p + theme(axis.text=element_text(size=7.5), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + ylab("Labels from paper") + xlab("Our labels")
p

p <-ggplot(data = d) + geom_count(mapping = aes(x = pre, y = celltype))
p <- p + guides(color=guide_legend(title = "Count")) + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12), axis.text.y = element_text(size=12),
               axis.text=element_text(size=7.5), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + ylab("Labels from paper") + xlab("Our labels")
p

# distance

d.after <- read.table("AML_res/BM4/BM4_RNN_after.csv", sep=",", header=TRUE)
d.ori <- read.table("AML_res/BM4/BM4_RNN_before.csv", sep=",", header=TRUE)

for (i in 1:32){
  if(max(d.ori[,i]) != 0){
    d.ori[,i] <- 100*d.ori[,i]/max(d.ori[,i])
  }
  if(max(d.after[,i]) != 0){
    d.after[,i] <- 100*d.after[,i]/max(d.after[,i])
  }
}

pre <- id
pre[id=="0"] <- "T"
pre[id=="1"] <- "Mono"
pre[id=="2"] <- "Prog & GMP"
pre[id=="3"] <- "Prog & GMP"
pre[id=="4"] <- "Prog & GMP"
pre[id=="5"] <- "NK"
pre[id=="6"] <- "Ery"
pre[id=="7"] <- "Plasma & B"
pre[id=="8"] <- "Plasma & B"

id.order <- order(pre)
pre <- pre[id.order,]

d <- d.ori
d <- d[id.order,]
d.dist <- dist(d)
d.dist <- as.matrix(d.dist)
d.dist <- as.data.frame(d.dist)
colnames(d.dist) <- c(1:ncol(d.dist))
rownames(d.dist) <- c(1:nrow(d.dist))
d.plot <- d.dist %>% mutate(B=colnames(d.dist)) %>% melt()
d.plot$B <- as.integer(d.plot$B)
d.plot$variable <- as.integer(d.plot$variable)
p <- ggplot(data = d.plot, aes(x=value)) + geom_freqpoly(bins = 250)
p
d.plot$type <- "Before"
d.all <- d.plot

d <- d.after
d <- d[id.order,]
d.dist <- dist(d)
d.dist <- as.matrix(d.dist)
d.dist <- as.data.frame(d.dist)
colnames(d.dist) <- c(1:ncol(d.dist))
rownames(d.dist) <- c(1:nrow(d.dist))
d.plot <- d.dist %>% mutate(B=colnames(d.dist)) %>% melt()
d.plot$B <- as.integer(d.plot$B)
d.plot$variable <- as.integer(d.plot$variable)
p <- ggplot(data = d.plot, aes(x=value)) + geom_freqpoly(bins = 250)
p
d.plot$type <- "After"
d.all <- rbind(d.all, d.plot)

d.all$group <- "within"
d.all[pre[d.all$B] != pre[d.all$variable],5] <- "between"

p <- ggplot(data = d.all, mapping = aes(x = value, y=..density.., color = type, linetype = group)) + geom_freqpoly(bins = 250)
p <- p + theme_bw() + theme(axis.text.x = element_text(size = 15),
                            axis.text.y = element_text(size = 15),
                            axis.title.x = element_text(size = 18),
                            axis.title.y = element_text(size = 18),
                            panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + xlab("Distance") + ylab("Density")
p

p <- ggplot(data = d.all[d.all$type=="After",], aes(x=B,y=variable)) + geom_tile(aes(fill = value))
p <- p + scale_fill_continuous(low = "blue", high = "yellow")
p <- p + guides(color=guide_legend(title = "Distance"))
p <- p + theme_bw() + theme(axis.title.x = element_text(size = 18),
                            axis.title.y = element_text(size = 18),
                            panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + xlab("Cells") + ylab("Cells")
p








