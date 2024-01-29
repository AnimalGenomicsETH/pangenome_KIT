setwd("/sotiria/Documents/")
install.packages("ape")
install.packages("tidyverse")
install.packages("phylotools")
#load library
library(ape)
library(tidyr)
library(dplyr)
library(stringr)


#read the distance matrix from the tsv file
datdis <- read.csv("6.dist.tsv", header=FALSE, stringsAsFactors=FALSE, sep="\t")

# Rename the header
colnames(datdis)  <- c("anim1","anim2","distr","comp4","comp5")

# Extract the correct assembly name from anim1 and anim2
datdis$anim1c <- str_extract(datdis$anim1, ".*")
datdis$anim2c <- str_extract(datdis$anim2, ".*")

# Make the distance matrix into a wide matrix
datsel  <- datdis  %>% select(anim1c, anim2c, distr)
datwide  <- datsel  %>% pivot_wider(names_from = anim2c, values_from = distr)
datmat  <- as.matrix(datwide  %>% select(-anim1c))
rownames(datmat)  <- datwide$anim1c

ref <- "ARS_UCD1.2#0#6"

outgroup  <- which.max(datmat[ref,])  %>% names()
print(outgroup)

# Apply neighbor joining
tr  <- nj(datmat)

# Visualize the tree
pdf("output.pdf", width=12, height=10)
plot.phylo(root(tr, outgroup=outgroup), cex=2, edge.width=2)
axisPhylo(backward=FALSE, cex.axis=2)
dev.off()
