

library(poppr)
library(vcfR)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(cowplot)
library(adegenet)
library(dartR)
library("genepopedit")
library(stringr)
library(hierfstat)
library(pegas)
library(ape)

#############################################################################################
#### define file names, software paths ####
work_dir <- "/Users/macbook2017/Desktop/dapc/ps/" 
setwd(work_dir)
all.gen  = "ps.11pop.neutral.17517.gen.txt"
all.renamed.gen="all.renamed.gen"

#############################################################################################
#### data manipulation, reformat genepop individual names to get population names automatically ####
genepop_ID(genepop=all.gen, path=paste0(work_dir, all.renamed.gen))  ##BJYQ01 > BJYQ_01

##get variable from genepop using genepop_detective
PopNames.all <- genepop_detective(all.renamed.gen, variable="Pops")
#PopCounts.all <- genepop_detective(all.renamed.gen, variable="PopNum")
#SampleIDs.all <- genepop_detective(all.renamed.gen, variable="Inds")
#LociNames.all <- genepop_detective(all.renamed.gen, variable="Loci")
#metadata.all <- genepop_detective(all.renamed.gen,variable="All")
#Alleles.all <- genepop_detective(all.renamed.gen,variable="Allele")

all.genind <- read.genepop(all.renamed.gen, ncode=3)
##set population to geneind file
SampleIDs <- genepop_detective(all.renamed.gen, variable="Inds")
SamplePop <- str_replace_all(SampleIDs, "_", "")
SamplePop <- str_replace_all(SamplePop, "0", "")
SamplePop <- str_replace_all(SamplePop, "1", "")
SamplePop <- str_replace_all(SamplePop, "2", "")
SamplePop <- str_replace_all(SamplePop, "3", "")
SamplePop <- str_replace_all(SamplePop, "4", "")
SamplePop <- str_replace_all(SamplePop, "5", "")
SamplePop <- str_replace_all(SamplePop, "6", "")
SamplePop <- str_replace_all(SamplePop, "7", "")
SamplePop <- str_replace_all(SamplePop, "8", "")
SamplePop <- str_replace_all(SamplePop, "9", "")
pop(all.genind) <- SamplePop
##subset populations for use from genind file##
PopNames.all
pop.sublist <- "ALL" #c("CQCQ", "GDGZ", "HNHK", "SCDY", "SCLS", "SXTY", "YNLJ", "YNYX") #  
pop.sublistname <- paste(pop.sublist, collapse ="_")
subset.genind = paste(pop.sublistname, ".genind", sep = "") 
subset.genind <- popsub(all.genind, sublist = pop.sublist) # sublist=1:10, sublist=1:10, blacklist="Bazadais", sublist=c(1:6, 11:15)

#############################################################################################
#### set the final used data in subsequent analysis ####
used.genind = subset.genind
#############################################################################################

## set color for figures and population names used
PopNames.used <- popNames(used.genind)
cols <- rainbow(nPop(used.genind))
#toRemove <- is.na(glMean(used.genlight, alleleAsUnit = T))
#used.genlight <- used.genlight[, !toRemove]

#### PCA from genelight file ####
pca.hierfstat <- indpca(used.hierfstat)  #Carry out a PCA on the centered matrix of individual's allele frequencies.
pdf(file=paste(pop.sublistname, ".pca.hierfstat.pdf", sep=""))
plot(pca.hierfstat, cex = 0.7)
dev.off()

glpca <- glPca(used.genlight,parallel = TRUE, n.cores = 7L,  nf = 2)
pdf(file=paste(pop.sublistname, ".pca.genlight.barplot.pdf", sep=""))
barplot(100*glpca$eig/sum(glpca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)
dev.off()

pdf(file=paste(pop.sublistname, ".pca.genlight.scatter2.pdf", sep=""))
scatter(glpca,posi="bottomleft")
dev.off()

glpca.scores <- as.data.frame(glpca$scores)
glpca.scores$pop <- pop(used.genlight)
library(ggplot2)
set.seed(9)

pdf(file=paste(pop.sublistname, ".pca.genlight.scatter.ggplot.pdf", sep=""))
p <- ggplot(glpca.scores, aes(x=PC1, y=PC2, colour=pop))  ##add x, y
p <- p + geom_point(size=2) ##add dot
p <- p + stat_ellipse(level = 0.95, size = 0.4) ## add circle
p <- p + scale_color_manual(values = cols) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()
p
dev.off()
