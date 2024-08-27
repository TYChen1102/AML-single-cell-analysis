test.data.0 <- read.table("AML_data/BM4/GSM3588000_BM4.dem.txt", sep="\t", header=TRUE)

test.data.0.anno <- read.table("AML_data/BM4/GSM3588001_BM4.anno.txt", sep="\t", header=TRUE)

test.data.0 <- as.data.frame(t(test.data.0))
colnames(test.data.0) <- test.data.0[1,]
test.data.0 <- test.data.0[2:nrow(test.data.0),]
test.data.0 <- apply(test.data.0, 2, as.numeric)
test.data.0 <- as.data.frame(test.data.0)

test.data.0 <- test.data.0[,which(colSums(test.data.0)>0)]

write.table(test.data.0, file='AML_data/BM4/BM4_data.csv', quote=FALSE, sep='\t', col.names = colnames(test.data.0))

test.data.0.anno <- test.data.0.anno$CellType

write.table(t(test.data.0.anno), file='AML_data/BM4/BM4_celltype.csv', quote=FALSE, sep='\t', col.names = NA)

unique(test.data.0.anno)

test.data.0.anno[test.data.0.anno == "GMP" | test.data.0.anno == "GMP-like"] <- "GMP"
test.data.0.anno[test.data.0.anno == "Prog" | test.data.0.anno == "Prog-like"] <- "Prog"
test.data.0.anno[test.data.0.anno == "ProMono" | test.data.0.anno == "ProMono-like"] <- "ProMono"
test.data.0.anno[test.data.0.anno == "cDC" | test.data.0.anno == "cDC-like"] <- "cDC"
test.data.0.anno[test.data.0.anno == "Mono" | test.data.0.anno == "Mono-like"] <- "Mono"
test.data.0.anno[test.data.0.anno == "HSC" | test.data.0.anno == "HSC-like"] <- "HSC"
test.data.0.anno[test.data.0.anno == "T"] <- "T"
test.data.0.anno[test.data.0.anno == "CTL"] <- "CTL"
test.data.0.anno[test.data.0.anno == "Plasma"] <- "Plasma"
test.data.0.anno[test.data.0.anno == "ProB"] <- "ProB"
test.data.0.anno[test.data.0.anno == "earlyEry"] <- "earlyEry"
test.data.0.anno[test.data.0.anno == "lateEry"] <- "lateEry"
test.data.0.anno[test.data.0.anno == "NK"] <- "NK"
test.data.0.anno[test.data.0.anno == "pDC"] <- "pDC"

write.table(t(test.data.0.anno), file='AML_data/BM4/BM4_celltype_merge.csv', quote=FALSE, sep='\t', col.names = NA)

test.name <- unique(test.data.0.anno)
test.name <- as.data.frame(test.name)
rownames(test.name) <- c(1:nrow(test.name))
test.name$number <- c(1:nrow(test.name))

for (i in 1:nrow(test.name)){test.data.0.anno[test.data.0.anno == test.name[i,1]] <- (test.name[i,2] - 1)}

write.table(t(test.data.0.anno), file='AML_data/BM4/BM4_celltype_number.csv', quote=FALSE, sep='\t', col.names = NA)
