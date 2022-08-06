test.data.1 <- read.table("AML_D0/GSM3587923_AML1012-D0.dem.txt", sep="\t", header=TRUE)
test_data.2 <- read.table("AML_D0/GSM3587925_AML210A-D0.dem.txt", sep="\t", header=TRUE)
test_data.3 <- read.table("AML_D0/GSM3587927_AML314-D0.dem.txt", sep="\t", header=TRUE)

test_data.4 <- read.table("AML-029-08-1E_data.csv", sep="\t", header=TRUE)
test_data.5 <- read.table("AML-029-08-1E_data.csv", sep="\t", header=TRUE)
test.data.6 <- read.table("AML-029-08-1E_data.csv", sep="\t", header=TRUE)
test_data.7 <- read.table("AML-029-08-1E_data.csv", sep="\t", header=TRUE)
test_data.8 <- read.table("AML_D0/GSM3587953_AML420B-D0.dem.txt", sep="\t", header=TRUE)
test_data.9 <- read.table("AML-029-08-1E_data.csv", sep="\t", header=TRUE)
test.data.10 <- read.table("AML_D0/GSM3587969_AML707B-D0.dem.txt", sep="\t", header=TRUE)
test.data.11 <- read.table("AML-029-08-1E_data.csv", sep="\t", header=TRUE)
test.data.12 <- read.table("AML_D0/GSM3587988_AML916-D0.dem.txt", sep="\t", header=TRUE)
test_data.13 <- read.table("AML-029-08-1E_data.csv", sep="\t", header=TRUE)
test.data.14 <- read.table("AML_D0/GSM3587959_AML475-D0.dem.txt", sep="\t", header=TRUE)
test_data.15 <- read.table("AML-029-08-1E_data.csv", sep="\t", header=TRUE)
test.data.16 <- read.table("AML_D0/GSM3587990_AML921A-D0.dem.txt", sep="\t", header=TRUE)

test.data.1.anno <- read.table("AML_D0_anno/GSM3587924_AML1012-D0.anno.txt", sep="\t", header=TRUE)
test.data.8.anno <- read.table("AML_D0_anno/GSM3587954_AML420B-D0.anno.txt", sep="\t", header=TRUE)
test.data.10.anno <- read.table("AML_D0_anno/GSM3587970_AML707B-D0.anno.txt", sep="\t", header=TRUE)
test.data.12.anno <- read.table("AML_D0_anno/GSM3587989_AML916-D0.anno.txt", sep="\t", header=TRUE)
test.data.14.anno <- read.table("AML_D0_anno/GSM3587960_AML475-D0.anno.txt", sep="\t", header=TRUE)
test.data.16.anno <- read.table("AML_D0_anno/GSM3587991_AML921A-D0.anno.txt", sep="\t", header=TRUE)

test.data.1 <- test.data.10
test.data.1.anno <- test.data.10.anno

test.data.1 <- as.data.frame(t(test.data.1))
colnames(test.data.1) <- test.data.1[1,]
test.data.1 <- test.data.1[2:nrow(test.data.1),]
test.data.1.new <- apply(test.data.1, 2, as.numeric)
test.data.1.new <- as.data.frame(test.data.1.new)

write.table(test.data.1.new, file='GSM3587969_AML707B-D0_data.csv', quote=FALSE, sep='\t', col.names = NA)

test.label.1 <- test.data.1.anno$CellType
write.table(t(test.label.1), file='GSM3587969_AML707B-D0_info.csv', quote=FALSE, sep='\t', col.names = NA)

unique(test.data.1.anno$CellType)
test.label.1[test.label.1 == "GMP" | test.label.1 == "GMP-like"] <- "GMP"
test.label.1[test.label.1 == "Prog" | test.label.1 == "Prog-like"] <- "Prog"
test.label.1[test.label.1 == "ProMono" | test.label.1 == "ProMono-like"] <- "ProMono"
test.label.1[test.label.1 == "cDC" | test.label.1 == "cDC-like"] <- "cDC"
test.label.1[test.label.1 == "Mono" | test.label.1 == "Mono-like"] <- "Mono"
test.label.1[test.label.1 == "HSC" | test.label.1 == "HSC-like"] <- "HSC"
test.label.1[test.label.1 == "T"] <- "T"
test.label.1[test.label.1 == "CTL"] <- "CTL"
test.label.1[test.label.1 == "Plasma"] <- "Plasma"
test.label.1[test.label.1 == "ProB"] <- "ProB"
test.label.1[test.label.1 == "earlyEry"] <- "earlyEry"
test.label.1[test.label.1 == "NK"] <- "NK"
test.label.1[test.label.1 == "pDC"] <- "pDC"

write.table(t(test.label.1), file='GSM3587969_AML707B-D0_info_new.csv', quote=FALSE, sep='\t', col.names = NA)

test.name <- unique(test.label.1)
test.name <- as.data.frame(test.name)
rownames(test.name) <- c(1:nrow(test.name))
test.name$number <- c(1:nrow(test.name))

for (i in 1:nrow(test.name)){test.label.1[test.label.1 == test.name[i,1]] <- test.name[i,2]}

write.table(t(test.label.1), file='GSM3587969_AML707B-D0_info_number.csv', quote=FALSE, sep='\t', col.names = NA)
