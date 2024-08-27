test.data.0 <- read.table("AML_data/AML707B/GSM3587969_AML707B-D0.dem.txt", sep="\t", header=TRUE)
test.data.18 <- read.table("AML_data/AML707B/GSM3587973_AML707B-D18.dem.txt", sep="\t", header=TRUE)
test.data.41 <- read.table("AML_data/AML707B/GSM3587975_AML707B-D41.dem.txt", sep="\t", header=TRUE)
test.data.97 <- read.table("AML_data/AML707B/GSM3587977_AML707B-D97.dem.txt", sep="\t", header=TRUE)
test.data.113 <- read.table("AML_data/AML707B/GSM3587971_AML707B-D113.dem.txt", sep="\t", header=TRUE)

test.data.0.anno <- read.table("AML_data/AML707B/GSM3587970_AML707B-D0.anno.txt", sep="\t", header=TRUE)
test.data.18.anno <- read.table("AML_data/AML707B/GSM3587974_AML707B-D18.anno.txt", sep="\t", header=TRUE)
test.data.41.anno <- read.table("AML_data/AML707B/GSM3587976_AML707B-D41.anno.txt", sep="\t", header=TRUE)
test.data.97.anno <- read.table("AML_data/AML707B/GSM3587978_AML707B-D97.anno.txt", sep="\t", header=TRUE)
test.data.113.anno <- read.table("AML_data/AML707B/GSM3587972_AML707B-D113.anno.txt", sep="\t", header=TRUE)

test.data.0 <- as.data.frame(t(test.data.0))
colnames(test.data.0) <- test.data.0[1,]
test.data.0 <- test.data.0[2:nrow(test.data.0),]
test.data.0 <- apply(test.data.0, 2, as.numeric)
test.data.0 <- as.data.frame(test.data.0)

test.data.18 <- as.data.frame(t(test.data.18))
colnames(test.data.18) <- test.data.18[1,]
test.data.18 <- test.data.18[2:nrow(test.data.18),]
test.data.18 <- apply(test.data.18, 2, as.numeric)
test.data.18 <- as.data.frame(test.data.18)

test.data.41 <- as.data.frame(t(test.data.41))
colnames(test.data.41) <- test.data.41[1,]
test.data.41 <- test.data.41[2:nrow(test.data.41),]
test.data.41 <- apply(test.data.41, 2, as.numeric)
test.data.41 <- as.data.frame(test.data.41)

test.data.97 <- as.data.frame(t(test.data.97))
colnames(test.data.97) <- test.data.97[1,]
test.data.97 <- test.data.97[2:nrow(test.data.97),]
test.data.97 <- apply(test.data.97, 2, as.numeric)
test.data.97 <- as.data.frame(test.data.97)

test.data.113 <- as.data.frame(t(test.data.113))
colnames(test.data.113) <- test.data.113[1,]
test.data.113 <- test.data.113[2:nrow(test.data.113),]
test.data.113 <- apply(test.data.113, 2, as.numeric)
test.data.113 <- as.data.frame(test.data.113)

test.data.all <- rbind(test.data.0, test.data.18, test.data.41, test.data.97, test.data.113)
test.data.all <- test.data.all[,which(colSums(test.data.all)>0)]

write.table(test.data.all, file='AML_data/AML707B/all/AML707B-all_data.csv', quote=FALSE, sep='\t', col.names = NA)

test.data.0.anno <- test.data.0.anno$CellType
test.data.18.anno <- test.data.18.anno$CellType
test.data.41.anno <- test.data.41.anno$CellType
test.data.97.anno <- test.data.97.anno$CellType
test.data.113.anno <- test.data.113.anno$CellType

test.data.0.size <- nrow(test.data.0)
test.data.18.size <- nrow(test.data.18)
test.data.41.size <- nrow(test.data.41)
test.data.97.size <- nrow(test.data.97)
test.data.113.size <- nrow(test.data.113)

test.data.all.anno <- c(test.data.0.anno, test.data.18.anno, test.data.41.anno, test.data.97.anno, test.data.113.anno)

write.table(t(test.data.all.anno), file='AML_data/AML707B/all/AML707B-all_celltype.csv', quote=FALSE, sep='\t', col.names = NA)

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

write.table(t(test.data.all.anno), file='AML_data/AML707B/all/AML707B-all_celltype_merge.csv', quote=FALSE, sep='\t', col.names = NA)

test.name <- unique(test.data.all.anno)
test.name <- as.data.frame(test.name)
rownames(test.name) <- c(1:nrow(test.name))
test.name$number <- c(1:nrow(test.name))

for (i in 1:nrow(test.name)){test.data.all.anno[test.data.all.anno == test.name[i,1]] <- (test.name[i,2] - 1)}

write.table(t(test.data.all.anno), file='AML_data/AML707B/all/AML707B-all_celltype_number.csv', quote=FALSE, sep='\t', col.names = NA)

test.data.all.sample <- test.data.all.anno
test.data.all.sample[1:1586] <- "D0"
test.data.all.sample[1587:3259] <- "D18"
test.data.all.sample[3260:3646] <- "D41"
test.data.all.sample[3647:3730] <- "D97"
test.data.all.sample[3731:4438] <- "D113"

write.table(t(test.data.all.sample), file='AML_data/AML707B/all/AML707B-all_celltype_sample.csv', quote=FALSE, sep='\t', col.names = NA)

write.table(t(colnames(test.data.all)), file='AML_data/AML707B/all/AML707B-all_genes.csv', quote=FALSE, sep='\t', col.names = NA)
