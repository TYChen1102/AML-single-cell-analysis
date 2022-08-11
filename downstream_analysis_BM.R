expression <- read.table("AML_data/BM4/BM4_data.csv", sep="\t", header=TRUE)
expression <- expression[,2:ncol(expression)]

latent <- read.table("AML_data/BM4/BM4_cluster.csv", sep=",", header=TRUE)
predict.labels <- read.table("AML_data/BM4/BM4_predict.csv", sep=",", header=TRUE)
truth.labels <- read.table("AML_data/BM4/BM4_celltype_merge.csv", sep="\t", header=TRUE)
truth.labels <- as.data.frame(t(truth.labels))
truth.labels <- truth.labels[2:nrow(truth.labels),]

latent$predict <- as.character(predict.labels$X0)
latent$truth <- as.character(truth.labels)
latent <- cbind(latent, expression)

unique(latent$predict)

latent$truth[latent$truth == "TRUE"] <- "T"

latent.sample <- latent[latent$truth == "B" | latent$truth == "lateEry" | latent$truth == "T" | latent$truth == "GMP" | latent$truth == "Mono" | latent$truth == "NK" | latent$truth == "Plasma",]

library(ggplot2)
library(dplyr)

#p <- latent %>% mutate(predict = factor(predict, levels = c(0:26))) %>% ggplot() + geom_point(mapping = aes(x = X0, y = X1, color = predict))
p <- ggplot(data = latent) + geom_point(mapping = aes(x = X0, y = X1, color = predict), size = 1)
p <- p + guides(color=guide_legend(title = "Clusters")) + theme_bw()
p <- p + theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + labs(x = "t-SNE1", y = "t-SNE2")
p

p <- ggplot(data = latent.sample) + geom_point(mapping = aes(x = X0, y = X1, color = predict), size = 1)
p <- p + guides(color=guide_legend(title = "Clusters")) + theme_bw()
p <- p + theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + labs(x = "t-SNE1", y = "t-SNE2")
p

p <- ggplot(data = latent) + geom_point(mapping = aes(x = X0, y = X1, color = truth), size = 1)
p <- p + guides(color=guide_legend(title = "Cells")) + theme_bw()
p <- p + theme(axis.text=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + labs(x = "t-SNE1", y = "t-SNE2")
p

p <- ggplot(data = latent.sample) + geom_point(mapping = aes(x = X0, y = X1, color = truth), size = 1)
p <- p + guides(color=guide_legend(title = "Cells")) + theme_bw()
p <- p + theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + labs(x = "t-SNE1", y = "t-SNE2")
p

# count
cell.count <- matrix(data = 0, nrow = max(predict.labels$X0) + 1, ncol = 1)
cell.count <- as.data.frame(cell.count)
for (i in 1:max(predict.labels$X0)+1)
{
  count0 <- nrow(expression[predict.labels$X0 == (i-1),])
  cell.count[i,] <- count0
}
count0 <- nrow(expression[predict.labels$X0 == 0,])
cell.count[1,] <- count0
colnames(cell.count) <- c("number")
rownames(cell.count) <- as.character(c(0:(nrow(cell.count)-1)))

# mean
expression.mean <- matrix(data = 0, nrow = max(predict.labels$X0) + 1, ncol = ncol(expression))
expression.mean <- as.data.frame(expression.mean)
colnames(expression.mean) <- colnames(expression)
rownames(expression.mean) <- rownames(cell.count)
for (i in 1:max(predict.labels$X0)+1)
{
  test0.mean <- apply(latent[latent$predict == as.character(i-1),5:ncol(latent)], 2, mean)
  expression.mean[i,] <- test0.mean
}
test0.mean <- apply(latent[latent$predict == "0",5:ncol(latent)], 2, mean)
expression.mean[1,] <- test0.mean

test0.mean.var <- apply(expression.mean, 2, var)
expression.mean <- rbind(expression.mean, test0.mean.var)
row.names(expression.mean)[nrow(expression.mean)] <- "var"
expression.mean.sorted <- expression.mean[,order(-expression.mean[nrow(expression.mean),])]

# gene expression on scatter plot
# see expression of some marker genes

ggplot(data = latent.sample) + geom_point(mapping = aes(x = X0, y = X1, color = HBB))
ggplot(data = latent.sample) + geom_point(mapping = aes(x = X0, y = X1, color = CD14))

# z_score
latent.sample[5:ncol(latent.sample)] <- scale(latent.sample[5:ncol(latent.sample)])

# z-score mean
expression.mean <- matrix(data = 0, nrow = max(predict.labels$X0) + 1, ncol = ncol(expression))
expression.mean <- as.data.frame(expression.mean)
colnames(expression.mean) <- colnames(expression)
rownames(expression.mean) <- rownames(cell.count)
for (i in 1:max(predict.labels$X0)+1)
{
  test0.mean <- apply(latent.sample[latent.sample$predict == as.character(i-1),5:ncol(latent.sample)], 2, mean)
  expression.mean[i,] <- test0.mean
}
test0.mean <- apply(latent.sample[latent.sample$predict == "0",5:ncol(latent.sample)], 2, mean)
expression.mean[1,] <- test0.mean
test0.mean.var <- apply(expression.mean, 2, var)
expression.mean <- rbind(expression.mean, test0.mean.var)
row.names(expression.mean)[nrow(expression.mean)] <- "var"
expression.mean.sorted <- expression.mean[,order(-expression.mean[nrow(expression.mean),])]

latent.sample$new_labels <- latent.sample$predict
latent.sample$new_labels[latent.sample$predict=="0"] <- "GMPs"
latent.sample$new_labels[latent.sample$predict=="1"] <- "T cells"
latent.sample$new_labels[latent.sample$predict=="2"] <- "GMPs"
latent.sample$new_labels[latent.sample$predict=="3"] <- "T cells"
latent.sample$new_labels[latent.sample$predict=="4"] <- "Monocytes"
latent.sample$new_labels[latent.sample$predict=="5"] <- "T cells"
latent.sample$new_labels[latent.sample$predict=="6"] <- "Erythrocytes"
latent.sample$new_labels[latent.sample$predict=="7"] <- "Monocytes"
latent.sample$new_labels[latent.sample$predict=="8"] <- "GMPs"
latent.sample$new_labels[latent.sample$predict=="9"] <- "Monocytes"
latent.sample$new_labels[latent.sample$predict=="10"] <- "Monocytes"
latent.sample$new_labels[latent.sample$predict=="11"] <- "Monocytes"
latent.sample$new_labels[latent.sample$predict=="12"] <- "Erythrocytes"
latent.sample$new_labels[latent.sample$predict=="13"] <- "Monocytes"
latent.sample$new_labels[latent.sample$predict=="14"] <- "T cells"
latent.sample$new_labels[latent.sample$predict=="15"] <- "GMPs"
latent.sample$new_labels[latent.sample$predict=="16"] <- "Plasma and B cells"
latent.sample$new_labels[latent.sample$predict=="17"] <- "GMPs"
latent.sample$new_labels[latent.sample$predict=="18"] <- "NK cells"
latent.sample$new_labels[latent.sample$predict=="19"] <- "T cells"
latent.sample$new_labels[latent.sample$predict=="20"] <- "GMPs"
latent.sample$new_labels[latent.sample$predict=="21"] <- "Monocytes"
latent.sample$new_labels[latent.sample$predict=="22"] <- "Erythrocytes"
latent.sample$new_labels[latent.sample$predict=="23"] <- "GMPs"

p <- ggplot(data = latent.sample) + geom_point(mapping = aes(x = X0, y = X1, color = new_labels), size = 1)
p <- p + guides(color=guide_legend(title = "Cells")) + theme_bw()
p <- p + theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + labs(x = "t-SNE1", y = "t-SNE2")
p

# small sample
latent.small <- latent.sample[,1:4]
latent.small$new_labels <- latent.sample$new_labels
select.genes <- colnames(expression.mean.sorted)[1:2000]
select.expression <- latent.sample[,select.genes]
latent.small <- cbind(latent.small, select.expression)

write.table(latent.small, file='AML_data/BM4/BM4_latent_small_zscore.csv', quote=FALSE, sep='\t', col.names = colnames(latent.small))

# confusion matrix
p <- latent.sample %>% mutate(truth = factor(truth, levels = c("lateEry","GMP","Mono","NK","Plasma","B","T"))) %>% ggplot() + geom_count(mapping = aes(x = new_labels, y = truth))
p <- p + guides(color=guide_legend(title = "Count")) + theme_bw()
p <- p + theme(axis.text=element_text(size=7.5), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + ylab("Labels from paper") + xlab("Our labels")
p

# visualize marker violin
p <- ggplot(data = latent.sample, mapping = aes(x = new_labels, y = CCR7, fill = new_labels)) + geom_violin()
p <- p + guides(color=guide_legend(title = "Cells")) + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text=element_text(size=15), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), panel.grid.major=element_line(colour=NA), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), panel.grid.minor = element_blank())
p <- p + labs(x = NULL)
p
