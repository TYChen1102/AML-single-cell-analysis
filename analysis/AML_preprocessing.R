test.data.0 <- read.table("AML_data/AML328/AML328-D0.dem.txt", sep="\t", header=TRUE)
test.data.29 <- read.table("AML_data/AML328/AML328-D29.dem.txt", sep="\t", header=TRUE)
test.data.113 <- read.table("AML_data/AML328/AML328-D113.dem.txt", sep="\t", header=TRUE)
test.data.171 <- read.table("AML_data/AML328/AML328-D171.dem.txt", sep="\t", header=TRUE)

test.data.0.anno <- read.table("AML_data/AML328/GSM3587932_AML328-D0.anno.txt", sep="\t", header=TRUE)
test.data.29.anno <- read.table("AML_data/AML328/GSM3587938_AML328-D29.anno.txt", sep="\t", header=TRUE)
test.data.113.anno <- read.table("AML_data/AML328/GSM3587934_AML328-D113.anno.txt", sep="\t", header=TRUE)
test.data.171.anno <- read.table("AML_data/AML328/GSM3587936_AML328-D171.anno.txt", sep="\t", header=TRUE)

test.data.0 <- as.data.frame(t(test.data.0))
colnames(test.data.0) <- test.data.0[1,]
test.data.0 <- test.data.0[2:nrow(test.data.0),]
test.data.0 <- apply(test.data.0, 2, as.numeric)
test.data.0 <- as.data.frame(test.data.0)

test.data.29 <- as.data.frame(t(test.data.29))
colnames(test.data.29) <- test.data.29[1,]
test.data.29 <- test.data.29[2:nrow(test.data.29),]
test.data.29 <- apply(test.data.29, 2, as.numeric)
test.data.29 <- as.data.frame(test.data.29)

test.data.113 <- as.data.frame(t(test.data.113))
colnames(test.data.113) <- test.data.113[1,]
test.data.113 <- test.data.113[2:nrow(test.data.113),]
test.data.113 <- apply(test.data.113, 2, as.numeric)
test.data.113 <- as.data.frame(test.data.113)

test.data.171 <- as.data.frame(t(test.data.171))
colnames(test.data.171) <- test.data.171[1,]
test.data.171 <- test.data.171[2:nrow(test.data.171),]
test.data.171 <- apply(test.data.171, 2, as.numeric)
test.data.171 <- as.data.frame(test.data.171)

test.data.all <- rbind(test.data.0, test.data.29, test.data.113, test.data.171)
test.data.all <- test.data.all[,which(colSums(test.data.all)>0)]

####
write.table(test.data.all, file='AML_data/AML328/AML328-all_data.csv', quote=FALSE, sep='\t', col.names = colnames(test.data.all))
####

test.data.0.anno <- test.data.0.anno$CellType
test.data.29.anno <- test.data.29.anno$CellType
test.data.113.anno <- test.data.113.anno$CellType
test.data.171.anno <- test.data.171.anno$CellType

test.data.0.size <- nrow(test.data.0)
test.data.29.size <- nrow(test.data.29)
test.data.113.size <- nrow(test.data.113)
test.data.171.size <- nrow(test.data.171)

test.data.all.anno <- c(test.data.0.anno, test.data.29.anno, test.data.113.anno, test.data.171.anno)

write.table(t(test.data.all.anno), file='AML_data/AML328/AML328-all_celltype.csv', quote=FALSE, sep='\t', col.names = NA)

unique(test.data.all.anno)

test.data.all.anno[test.data.all.anno == "GMP" | test.data.all.anno == "GMP-like"] <- "GMP"
test.data.all.anno[test.data.all.anno == "Prog" | test.data.all.anno == "Prog-like"] <- "Prog"
test.data.all.anno[test.data.all.anno == "ProMono" | test.data.all.anno == "ProMono-like"] <- "ProMono"
test.data.all.anno[test.data.all.anno == "cDC" | test.data.all.anno == "cDC-like"] <- "cDC"
test.data.all.anno[test.data.all.anno == "Mono" | test.data.all.anno == "Mono-like"] <- "Mono"
test.data.all.anno[test.data.all.anno == "HSC" | test.data.all.anno == "HSC-like"] <- "HSC"
test.data.all.anno[test.data.all.anno == "T"] <- "T"
test.data.all.anno[test.data.all.anno == "CTL"] <- "CTL"
test.data.all.anno[test.data.all.anno == "Plasma"] <- "Plasma"
test.data.all.anno[test.data.all.anno == "ProB"] <- "ProB"
test.data.all.anno[test.data.all.anno == "earlyEry"] <- "earlyEry"
test.data.all.anno[test.data.all.anno == "lateEry"] <- "lateEry"
test.data.all.anno[test.data.all.anno == "NK"] <- "NK"
test.data.all.anno[test.data.all.anno == "pDC"] <- "pDC"

####
write.table(t(test.data.all.anno), file='AML_data/AML328/AML328-all_celltype_merge.csv', quote=FALSE, sep='\t', col.names = NA)
####

test.name <- unique(test.data.all.anno)
test.name <- as.data.frame(test.name)
rownames(test.name) <- c(1:nrow(test.name))
test.name$number <- c(1:nrow(test.name))

for (i in 1:nrow(test.name)){test.data.all.anno[test.data.all.anno == test.name[i,1]] <- (test.name[i,2] - 1)}

write.table(t(test.data.all.anno), file='AML_data/AML328/AML328-all_celltype_number.csv', quote=FALSE, sep='\t', col.names = NA)

test.data.all.sample <- test.data.all.anno
test.data.all.sample[1:1094] <- "D0"
test.data.all.sample[1095:2974] <- "D29"
test.data.all.sample[2975:5003] <- "D113"
test.data.all.sample[5004:6405] <- "D171"
test.days <- c(rep("D0", 1094), rep("D29", 1880), rep("D113", 2029), rep("D171", 1402))

#####
write.table(t(test.data.all.sample), file='AML_data/AML328/AML328-all_days_1.csv', quote=FALSE, sep='\t', col.names = NA)
#####

#write.table(t(colnames(test.data.all)), file='AML_data/AML707B/all/AML707B-all_genes.csv', quote=FALSE, sep='\t', col.names = NA)

test.data <- read.table("AML_data/AML707B/all/AML707B-all_normalize.csv", sep=",", header=TRUE)
test.label <- read.table("AML_data/AML707B/all/AML707B-all_celltype.csv", sep="\t", header=TRUE)
test.label <- subset(test.label, select = -c(X))

data.selected <- sample(1:nrow(test.data), floor(nrow(test.data)/3), replace = FALSE)
data.selected <- sort(data.selected, decreasing = FALSE)

test.data <- test.data[data.selected,]
test.label <- test.label[,data.selected]

data.selected <- sample(1:ncol(test.data), floor(ncol(test.data)/3), replace = FALSE)
data.selected <- sort(data.selected, decreasing = FALSE)

test.data <- test.data[,data.selected]

write.table(test.data, file='AML_data/AML707B/all/AML707B-all_normalize_select.csv', quote=FALSE, sep='\t', col.names = NA)
write.table(test.label, file='AML_data/AML707B/all/AML707B-all_celltype_select.csv', quote=FALSE, sep='\t', col.names = NA)



test.data.0.anno <- read.table("AML_data/AML707B/GSM3587970_AML707B-D0.anno.txt", sep="\t", header=TRUE)
test.data.18.anno <- read.table("AML_data/AML707B/GSM3587974_AML707B-D18.anno.txt", sep="\t", header=TRUE)
test.data.41.anno <- read.table("AML_data/AML707B/GSM3587976_AML707B-D41.anno.txt", sep="\t", header=TRUE)
test.data.97.anno <- read.table("AML_data/AML707B/GSM3587978_AML707B-D97.anno.txt", sep="\t", header=TRUE)
test.data.113.anno <- read.table("AML_data/AML707B/GSM3587972_AML707B-D113.anno.txt", sep="\t", header=TRUE)
test.days <- c(rep("D0", 1586), rep("D18", 1673), rep("D41", 387), rep("D97", 84), rep("D113", 708))
write.table(test.days, file='AML_data/AML707B/AML707B-all_days_1.csv', quote=FALSE, sep='\t', col.names = NA)

test.data.0.anno <- test.data.0.anno$CellType
test.data.18.anno <- test.data.18.anno$CellType
test.data.41.anno <- test.data.41.anno$CellType
test.data.97.anno <- test.data.97.anno$CellType
test.data.113.anno <- test.data.113.anno$CellType

test.data.all.anno <- c(test.data.0.anno, test.data.18.anno, test.data.41.anno, test.data.97.anno, test.data.113.anno)
cells[cells=="TRUE"] <- "T"


test.data.0.anno <- read.table("AML_data/BM4/GSM3588001_BM4.anno.txt", sep="\t", header=TRUE)
test.data.0.anno <- read.table("AML_data/BM4/BM4.dem.txt", sep="\t", header=TRUE)
rownames(test.data.0.anno) <- test.data.0.anno$Gene
test.data.0.anno <- test.data.0.anno[,2:ncol(test.data.0.anno)]
test.data.0.anno <- as.data.frame(t(test.data.0.anno))

