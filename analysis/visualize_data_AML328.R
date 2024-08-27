library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)
library(Rtsne)
library(uwot)

d <- read.table("AML_res/AML328/AML328_RNN_after_2.csv", sep=",", header=TRUE)
id <- read.table("AML_res/AML328/AML328_predict.csv", sep=",", header=TRUE)
nclusters <- max(id$X0) + 1
id$X0 <- as.character(id$X0)

# t-SNE plot

tsne <- Rtsne(d, pca_scale = FALSE)
df.tsne <- data.frame(tsne$Y)
#df.tsne <- read.table("AML_data/AML328/AML328_tsne_res.csv", sep="\t", header=TRUE)
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

raw <- read.table("AML_data/AML328/AML328_all_data.csv", sep="\t", header=TRUE)
cells <- read.table("AML_data/AML328/AML328_idents.csv", sep="\t", header=TRUE)
#days <- read.table("AML_data/AML328/AML328-all_days_1.csv", sep="\t", header=TRUE)
cells <- cells[2:ncol(cells)]
cells <- t(cells)
#days <- days[2:ncol(days)]
#days <- t(days)
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

#days1 <- days[cells == "lateEry"]
#days2 <- days[cells == "TRUE"]
#days3 <- days[cells == "B"]
#days4 <- days[cells == "Mono"]
#days5 <- days[cells == "NK"]
#days6 <- days[cells == "Plasma"]
#days7 <- days[cells == "GMP"]
#days8 <- days[cells == "Prog"]
#days <- c(days1,days2,days3,days4,days5,days6,days7,days8)

se <- seq(1,nrow(raw),2)
raw <- raw[se,]
cells <- cells[se]
d <- cbind(raw, df.tsne)

#days <- days[se]

ggplot(d, aes(x=tSNE1,y=tSNE2, color=celltype)) + geom_point()

p <- ggplot(data = d) + geom_point(mapping = aes(x = tSNE1,y = tSNE2, color = celltype), size = 1)
p <- p + guides(color=guide_legend(title = "Clusters")) + theme_bw()
p <- p + theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + labs(x = "t-SNE1", y = "t-SNE2")
p

# DEG

d[,1:(ncol(d)-4)] <- scale(d[,1:(ncol(d)-4)])
exp.mean <- matrix(data = 0, nrow = nclusters, ncol = ncol(d)-4)
exp.mean <- as.data.frame(exp.mean)
colnames(exp.mean) <- colnames(d)[1:(ncol(d)-4)]
for (i in 1:nclusters){
  m <- apply(d[d$id == as.character(i-1), 1:(ncol(d)-4)], 2, mean)
  exp.mean[i,] <- m
}

genelist <- c()
for (i in 1:nclusters){
  m <- exp.mean[i,]
  #m <- sort(m[1,], decreasing = TRUE)
  #g <- colnames(m)[1:50]
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
d$pre[d$id=="1"] <- "T"
d$pre[d$id=="2"] <- "NK"
d$pre[d$id=="3"] <- "Prog & GMP"
d$pre[d$id=="4"] <- "Prog & GMP"
d$pre[d$id=="5"] <- "Mono"
d$pre[d$id=="6"] <- "Prog & GMP"
d$pre[d$id=="7"] <- "T"
d$pre[d$id=="8"] <- "Ery"
d$pre[d$id=="9"] <- "Plasma & B"

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

used.marker <- c("HBB", "AHSP", "CA2", "HBD", "CLEC11A", "PRSS57", "EGFL7", "SOX4", "CD34", "LYZ", "CD14", "GNLY", "GZMB", "GZMA", "NKG7", "MS4A1", "FCRLA", "FCRL5", "MZB1", "IL7R", "CD3D", "CCR7")
used.marker <- c("HBB", "AHSP", "CA2", "HBD", "CLEC11A", "EGFL7", "CD34", "LYZ", "CD14", "NKG7", "GZMB", "MS4A1", "FCRL5", "MZB1", "CD3D", "CCR7")
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

# change

#days <- c(rep("D0",1904),rep("D29",1880),rep("D113",2029),rep("D171",1402))
days <- read.table("AML_data/AML328/AML328-all_days_1.csv", sep="\t", header=TRUE)
days <- days[2:ncol(days)]
days <- t(days)
cells <- read.table("AML_data/AML328/AML328_idents.csv", sep="\t", header=TRUE)
cells <- cells[2:ncol(cells)]
cells <- t(cells)

days1 <- days[cells == "lateEry"]
days2 <- days[cells == "TRUE"]
days3 <- days[cells == "B"]
days4 <- days[cells == "Mono"]
days5 <- days[cells == "NK"]
days6 <- days[cells == "Plasma"]
days7 <- days[cells == "GMP"]
days8 <- days[cells == "Prog"]
days <- c(days1,days2,days3,days4,days5,days6,days7,days8)
days <- days[se]

d$days <- days
d$d0 <- d$days == "D0"
d$d29 <- d$days == "D29"
d$d113 <- d$days == "D113"
d$d171 <- d$days == "D171"

sum(d$d0)
sum(d$d29)
sum(d$d113)
sum(d$d171)

p <- ggplot(data = d) + geom_point(mapping = aes(x = tSNE1,y = tSNE2, color = d0), size = 1)
p <- p + guides(color=guide_legend(title = NULL)) + theme_bw()
p <- p + theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + labs(x = "t-SNE1", y = "t-SNE2")
p <- p + scale_colour_manual(values=c("red", "gray"),breaks = c("TRUE","FALSE"), labels = c("D0 (398 cells)","2611 cells"))
p

p <- ggplot(data = d) + geom_point(mapping = aes(x = tSNE1,y = tSNE2, color = d29), size = 1)
p <- p + guides(color=guide_legend(title = NULL)) + theme_bw()
p <- p + theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + labs(x = "t-SNE1", y = "t-SNE2")
p <- p + scale_colour_manual(values=c("red", "gray"),breaks = c("TRUE","FALSE"), labels = c("D29 (650 cells)","2611 cells"))
p

p <- ggplot(data = d) + geom_point(mapping = aes(x = tSNE1,y = tSNE2, color = d113), size = 1)
p <- p + guides(color=guide_legend(title = NULL)) + theme_bw()
p <- p + theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + labs(x = "t-SNE1", y = "t-SNE2")
p <- p + scale_colour_manual(values=c("red", "gray"),breaks = c("TRUE","FALSE"), labels = c("D113 (941 cells)","2611 cells"))
p

p <- ggplot(data = d) + geom_point(mapping = aes(x = tSNE1,y = tSNE2, color = d171), size = 1)
p <- p + guides(color=guide_legend(title = NULL)) + theme_bw()
p <- p + theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + labs(x = "t-SNE1", y = "t-SNE2")
p <- p + scale_colour_manual(values=c("red", "gray"),breaks = c("TRUE","FALSE"), labels = c("D171 (622 cells)","2611 cells"))
p

a <- which(colnames(d) == "pre")
b <- which(colnames(d) == "days")

days.matrix <- table(d[,c(a,b)])
days.matrix <- as.data.frame(days.matrix)
days.matrix <-within(days.matrix,{days <-factor(days,levels=c("D0","D29","D113","D171"))})
days.matrix$Freq[days.matrix$days == "D0"] <- days.matrix$Freq[days.matrix$days == "D0"]/sum(d$d0)
days.matrix$Freq[days.matrix$days == "D29"] <- days.matrix$Freq[days.matrix$days == "D29"]/sum(d$d29)
days.matrix$Freq[days.matrix$days == "D113"] <- days.matrix$Freq[days.matrix$days == "D113"]/sum(d$d113)
days.matrix$Freq[days.matrix$days == "D171"] <- days.matrix$Freq[days.matrix$days == "D171"]/sum(d$d171)

p <- ggplot(data = days.matrix, aes(x = days, y = Freq,  color = pre, group = pre))
p <- p + geom_point(size = 3)
p <- p + geom_line(size = 1)
p <- p + guides(color=guide_legend(title = "Cells")) + theme_bw()
p <- p + theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15))
p <- p + labs(x = "Days", y = "Proportion (%)")
p

# mutations

library(stringr)

test.data.0.anno <- read.table("AML_data/AML328/GSM3587932_AML328-D0.anno.txt", sep="\t", header=TRUE)
test.data.29.anno <- read.table("AML_data/AML328/GSM3587938_AML328-D29.anno.txt", sep="\t", header=TRUE)
test.data.113.anno <- read.table("AML_data/AML328/GSM3587934_AML328-D113.anno.txt", sep="\t", header=TRUE)
test.data.171.anno <- read.table("AML_data/AML328/GSM3587936_AML328-D171.anno.txt", sep="\t", header=TRUE)

test.data.0.anno <- test.data.0.anno$CellType
test.data.29.anno <- test.data.29.anno$CellType
test.data.113.anno <- test.data.113.anno$CellType
test.data.171.anno <- test.data.171.anno$CellType

test.data.all.anno <- c(test.data.0.anno, test.data.29.anno, test.data.113.anno, test.data.171.anno)

anno1 <- test.data.all.anno[cells == "lateEry"]
anno2 <- test.data.all.anno[cells == "TRUE"]
anno3 <- test.data.all.anno[cells == "B"]
anno4 <- test.data.all.anno[cells == "Mono"]
anno5 <- test.data.all.anno[cells == "NK"]
anno6 <- test.data.all.anno[cells == "Plasma"]
anno7 <- test.data.all.anno[cells == "GMP"]
anno8 <- test.data.all.anno[cells == "Prog"]
anno <- c(anno1,anno2,anno3,anno4,anno5,anno6,anno7,anno8)

anno <- anno[se]
d <- cbind(d, anno)

mut <- anno

mut[1:length(anno)] <- "Normal"
mut[str_detect(anno, "like")] <- "Malignant"
d <- cbind(d, mut)

p <- ggplot(data = d) + geom_point(mapping = aes(x = tSNE1,y = tSNE2, color = mut), size = 1)
p <- p + guides(color=guide_legend(title = "Cells")) + theme_bw()
p <- p + theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + labs(x = "t-SNE1", y = "t-SNE2")
p

a <- which(colnames(d) == "days")
b <- which(colnames(d) == "pre")
c <- which(colnames(d) == "mut")

d.mut <- table(d[,c(a,b,c)])
d.mut <- as.data.frame(d.mut)
d.mut <-within(d.mut,{days <-factor(days,levels=c("D0","D29","D113","D171"))})
d.mut <-within(d.mut,{pre <-factor(pre,levels=c("Prog & GMP","Ery","NK","T","Plasma & B","Mono"))})
d.mut$count <- d.mut$Freq

for (i in 1:nrow(d.mut)){
  a <- which(d.mut$days == d.mut$days[i] & d.mut$pre == d.mut$pre[i])
  c <- d.mut[a[1],4] + d.mut[a[2],4]
  b <- d.mut[i,4]/c
  d.mut[i,5] <- c
  if (c != 0){
    d.mut[i,4] <- b
  }
}

d.mut <- d.mut[d.mut$mut == "Malignant",]

p <-ggplot(data = d.mut) + geom_point(mapping = aes(x = days, y = pre, color = Freq, size = count))
p <- p + guides(color=guide_legend(title = "Count")) + theme_bw()
p <- p + theme(axis.text=element_text(size=12), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + ylab("Cell types") + xlab("Samples")
p

p <- ggplot(data = d) + geom_bar(aes(x = pre, fill = mut), stat = "count", position = "stack")
p <- p + guides(color=guide_legend(title = "Cells")) + theme_bw()
p <- p + theme(axis.text.x = element_text(size = 10),
               axis.text.y = element_text(size=12), 
               axis.title=element_text(size=15), 
               legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + ylab("Number of cells") + xlab("Cell types")
p

d <-within(d,{days <-factor(days,levels=c("D0","D29","D113","D171"))})

p <- ggplot(data = d) + geom_bar(aes(x = days, fill = mut), stat = "count", position = "stack")
p <- p + guides(color=guide_legend(title = "Cells")) + theme_bw()
p <- p + theme(axis.text.x = element_text(size = 12),
               axis.text.y = element_text(size=12), 
               axis.title=element_text(size=15), 
               legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + ylab("Number of cells") + xlab("Days")
p

# distance

d.after <- read.table("AML_res/AML328/AML328_RNN_after_2.csv", sep=",", header=TRUE)
d.ori <- read.table("AML_res/AML328/AML328_RNN_before_2.csv", sep=",", header=TRUE)

id <- read.table("AML_res/AML328/AML328_predict.csv", sep=",", header=TRUE)

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
pre[id=="1"] <- "T"
pre[id=="2"] <- "NK"
pre[id=="3"] <- "Prog & GMP"
pre[id=="4"] <- "Prog & GMP"
pre[id=="5"] <- "Mono"
pre[id=="6"] <- "Prog & GMP"
pre[id=="7"] <- "T"
pre[id=="8"] <- "Ery"
pre[id=="9"] <- "Plasma & B"

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
#d.plot$value <- 100*d.plot$value/max(d.plot$value)
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
#d.plot$value <- 100*d.plot$value/max(d.plot$value)
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
                            legend.text = element_text(size = 15),
                            legend.title = element_text(size = 15),
                            panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + xlab("Distance") + ylab("Density")
p

p <- ggplot(data = d.all[d.all$type=="After",], aes(x=B,y=variable)) + geom_tile(aes(fill = value))
p <- p + scale_fill_continuous(low = "blue", high = "yellow")
p <- p + guides(color=guide_legend(title = "Distance"))
p <- p + theme_bw() + theme(axis.title.x = element_text(size = 18),
                            axis.title.y = element_text(size = 18),
                            axis.text.x = element_text(size = 15),
                            axis.text.y = element_text(size = 15),
                            legend.text = element_text(size = 15),
                            legend.title = element_text(size = 15),
                            panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + xlab("Cells") + ylab("Cells")
p








