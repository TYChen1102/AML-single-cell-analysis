expression <- read.table("GSM3587923_AML1012-D0_normalize.csv", sep=",", header=TRUE)

expression.data <- read.table("GSM3587923_AML1012-D0_data.csv", sep="\t", header=TRUE)

row.names(expression.data) <- expression.data$X
expression.data <- subset(expression.data, select = -c(X))

rownames(expression) <- rownames(expression.data)
colnames(expression) <- colnames(expression.data)

latent <- read.table("GSM3587923_AML1012-D0_clusters.csv", sep=",", header=TRUE)
predict.labels <- read.table("GSM3587923_AML1012-D0_predict.csv", sep=",", header=TRUE)
truth.labels <- read.table("GSM3587923_AML1012-D0_info_new.csv", sep="\t", header=TRUE)
truth.labels <- as.data.frame(t(truth.labels))
truth.labels <- truth.labels[2:nrow(truth.labels),]

latent$predict <- as.character(predict.labels$X0)
latent$truth <- as.character(truth.labels)
latent <- cbind(latent, expression)

library(ggplot2)
ggplot(data = latent) + geom_point(mapping = aes(x = X0, y = X1, color = predict)) #+ guides(color=guide_legend(title = "truth"))
ggplot(data = latent) + geom_point(mapping = aes(x = X0, y = X1, color = EEF1A1))

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

# sum
expression.sum <- matrix(data = 0, nrow = max(predict.labels$X0) + 1, ncol = ncol(expression))
expression.sum <- as.data.frame(expression.sum)
colnames(expression.sum) <- colnames(expression)
rownames(expression.sum) <- rownames(cell.count)
for (i in 1:max(predict.labels$X0)+1)
{
  test0.sum <- apply(latent[latent$predict == i-1,4:ncol(latent)], 2, sum)
  expression.sum[i,] <- test0.sum
}
test0.sum <- apply(latent[latent$predict == 0,4:ncol(latent)], 2, sum)
expression.sum[1,] <- test0.sum

test0.sum.var <- apply(expression.sum, 2, var)
expression.sum <- rbind(expression.sum, test0.sum.var)
row.names(expression.sum)[nrow(expression.sum)] <- "var"
expression.sum.sorted <- expression.sum[,order(-expression.sum[nrow(expression.sum),])]

# mean
expression.mean <- matrix(data = 0, nrow = max(predict.labels$X0) + 1, ncol = ncol(expression))
expression.mean <- as.data.frame(expression.mean)
colnames(expression.mean) <- colnames(expression)
rownames(expression.mean) <- rownames(cell.count)
for (i in 1:max(predict.labels$X0)+1)
{
  test0.mean <- apply(latent[latent$predict == i-1,4:ncol(latent)], 2, mean)
  expression.mean[i,] <- test0.mean
}
test0.mean <- apply(latent[latent$predict == 0,4:ncol(latent)], 2, mean)
expression.mean[1,] <- test0.mean

test0.mean.var <- apply(expression.mean, 2, var)
expression.mean <- rbind(expression.mean, test0.mean.var)
row.names(expression.mean)[nrow(expression.mean)] <- "var"
expression.mean.sorted <- expression.mean[,order(-expression.mean[nrow(expression.mean),])]

# gene expression on scatter plot
ggplot(data = latent) + geom_point(mapping = aes(x = X0, y = X1, color = RPS6))

# gene expression density for each cluster

# predict
ggplot(data = latent, mapping = aes(x = RPS6, y=..density.., color = predict)) + geom_freqpoly(binwidth = 5)
ggplot(data = latent,aes(x = predict, y = LYZ, fill = predict)) + geom_violin()

# truth
ggplot(data = latent, mapping = aes(x = RPS6, y=..density.., color = truth)) + geom_freqpoly(binwidth = 5)
ggplot(data = latent,aes(x = truth, y = LYZ, fill = truth)) + geom_violin()
