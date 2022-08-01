drug_test <- read.table("AML-029-08-1E_data.csv", sep="\t", header=TRUE)
drug_normalize <- read.table("AML-029-08-1E_data.csv", sep="\t", header=TRUE)

row.names(drug_test) <- drug_test$X
drug_test <- subset(drug_test, select = -c(X))

drug_ed <- read.table("AML-029-08-1E_ed.csv", sep=",", header=TRUE)
drug_label <- read.table("AML-029-08-1E_label.csv", sep=",", header=TRUE)

genes <- read.table("AML-029-08-1E.DEGs.for.all.celltypes.tsv", sep="\t", header=TRUE)
genes <- rownames(genes)
drug_test_genes <- drug_test[,colnames(drug_test) %in% genes]

#drug_test <- log2(1+drug_test)
#drug_normalize <- log2(1+drug_normalize)
#colnames(drug_normalize) <- colnames(drug_test)

drug_ed$label <- as.character(drug_label$X0)
drug_ed_1 <- cbind(drug_ed, drug_test_genes)

library(ggplot2)
ggplot(data = drug_ed) + geom_point(mapping = aes(x = X0, y = X1, color = label))
ggplot(data = drug_ed_1) + geom_point(mapping = aes(x = X0, y = X1, color = IL7R))

drug_test <- drug_normalize

test0 <- drug_test[drug_label == 0,]
test0_sum <- apply(test0, 2, mean)
test0 <- rbind(test0, test0_sum)
rownames(test0)[nrow(test0)] <- "average"
ggplot(data = drug_ed_1) + geom_point(mapping = aes(x = X0, y = X1, color = EEF1A1))

# count
cell_count <- matrix(data = 0, nrow = max(drug_label) + 1, ncol = 1)
cell_count <- as.data.frame(cell_count)
for (i in 1:max(drug_label)+1)
{
  count0 <- nrow(drug_test[drug_label == (i-1),])
  cell_count[i,] <- count0
}
count0 <- nrow(drug_test[drug_label == 0,])
cell_count[1,] <- count0
colnames(cell_count) <- c("number")
rownames(cell_count) <- as.character(c(0:(nrow(cell_count)-1)))

# sum
test_sum <- matrix(data = 0, nrow = max(drug_label) + 1, ncol = ncol(drug_ed_1) - 3)
test_sum <- as.data.frame(test_sum)
for (i in 1:max(drug_label)+1)
{
  test0_sum <- apply(drug_ed_1[drug_ed_1$label == i-1,4:ncol(drug_ed_1)], 2, sum)
  test_sum[i,] <- test0_sum
}
colnames(test_sum) <- colnames(drug_ed_1)[4:ncol(drug_ed_1)]
rownames(test_sum) <- rownames(cell_count)
test0_sum <- apply(drug_ed_1[drug_ed_1$label == 0,4:ncol(drug_ed_1)], 2, sum)
test_sum[1,] <- test0_sum

test0_sum_var <- apply(test_sum, 2, var)
test_sum <- rbind(test_sum, test0_sum_var)
row.names(test_sum)[nrow(test_sum)] <- "var"
test_sum_sorted <- test_sum[,order(-test_sum[nrow(test_sum),])]

#test_sum_sorted <- test_sum[,order(test_sum[1,])]

# average
test_avg <- matrix(data = 0, nrow = max(drug_label) + 1, ncol = ncol(drug_ed_1) - 3)
test_avg <- as.data.frame(test_avg)
for (i in 1:max(drug_label)+1)
{
  test0_avg <- apply(drug_ed_1[drug_ed_1$label == (i-1),4:ncol(drug_ed_1)], 2, mean)
  test_avg[i,] <- test0_avg
}
colnames(test_avg) <- colnames(drug_ed_1)[4:ncol(drug_ed_1)]
rownames(test_avg) <- rownames(cell_count)
test0_avg <- apply(drug_ed_1[drug_ed_1$label == 0,4:ncol(drug_ed_1)], 2, mean)
test_avg[1,] <- test0_avg

test0_avg_var <- apply(test_avg, 2, var)
test_avg <- rbind(test_avg, test0_avg_var)
row.names(test_avg)[nrow(test_avg)] <- "var"
test_avg_sorted <- test_avg[,order(-test_avg[nrow(test_avg),])]

ggplot(data = drug_ed_1) + geom_point(mapping = aes(x = X0, y = X1, color = IL7R))

label_new <- drug_label
label_new[drug_label == 9] <- "CD4+ effector"
label_new[drug_label == 12] <- "CD4+ naive"
label_new[drug_label == 6] <- "CD4+ naive"
#label_new[drug_label == 4] <- "CD4+ naive"
label_new[drug_label == 8] <- "CD4+ Tcm"
label_new[drug_label == 11] <- "CD8+ effector"
label_new[drug_label == 15] <- "Plasma cells"
label_new[drug_label == 0] <- "Erythrocytes"
label_new[drug_label == 1] <- "Erythrocytes"

drug_ed$new_label <- label_new$X0
ggplot(data = drug_ed) + geom_point(mapping = aes(x = X0, y = X1, color = new_label))
