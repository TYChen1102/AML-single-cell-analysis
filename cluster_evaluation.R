#latent_data <- read.table("AML-029-08-1E_ed.csv", sep=",", header=TRUE)
#predict_label <- read.table("AML-029-08-1E_label.csv", sep=",", header=TRUE)
truth_label <- read.table("colon_cell_info.csv", sep="\t", header=TRUE)

dist <- read.table("colon_dist.csv", sep=",", header=TRUE)

library(ggplot2)
library(reshape2)
library(dplyr)

# heatmap
dist_plot <- data.frame(dist[1:400,1:400])
colnames(dist_plot) <- c(1:400)
rownames(dist_plot) <- c(1:400)
dist_plot <- dist_plot %>% mutate(B=colnames(dist_plot)) %>% melt()
dist_plot$B <- as.integer(dist_plot$B)
dist_plot$variable <- as.integer(dist_plot$variable)
p1 <- ggplot(dist_plot,aes(x=variable,y=B,fill=value)) 
p2 <- p1 + geom_raster() + scale_fill_gradient2(low="yellow", high="blue", mid="white")
p2
p2 + xlab("Cells") + ylab("Cells")

# histogram
ggplot(data = dist_plot, mapping = aes(x = value)) + geom_histogram(binwidth = 0.1)
