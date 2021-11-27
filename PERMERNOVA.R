setwd("D:\\CNHS\\results_paper\\PCOA")

library(permute)
library(lattice)
library(vegan)
library(ggplot2)
library(plyr)

otu <- read.delim('OTU.txt', row.names = 1, sep = '\t', stringsAsFactors = F, check.names = F)
group <- read.delim('group.txt', sep = '\t', stringsAsFactors = F)

distance <- vegdist(otu, method = 'bray')
bray_curtis <- as.matrix(distance)

design <- group

idx = design$V1 %in% colnames(bray_curtis) 
sub_design = design[idx,]
bray_curtis = bray_curtis[rownames(sub_design), rownames(sub_design$V1)] # subset and reorder distance matrix

pcoa = cmdscale(bray_curtis, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points = as.data.frame(pcoa$points) # get coordinate string, format to dataframme
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig
points = cbind(points, sub_design[match(rownames(points), sub_design$eid), ])
fix(design)

design$V2 <- as.factor(design$District)

p = ggplot(points, aes(x=x, y=y, color=design$V2)) +
  geom_point(alpha=.7, size=2) + 
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title="PCoA of Bray-Curtis distance of selected microbial composition")+stat_ellipse(level = 0.85, show.legend = F, size = 1.3)
p
ggsave("beta_diversity_north_south.pdf", p)

group$District <- as.factor(group$District)

model1_perm <- adonis(distance~group$District, data = group, permutations = 1000, method="bray")

model1_perm
