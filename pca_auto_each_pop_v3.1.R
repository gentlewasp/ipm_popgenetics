### do PCA for each pop

library(poppr)
library(RColorBrewer)
library(ggplot2)
library(adegenet)
library(genepopedit)
library(stringr)
library(dartR)

#############################################################################################
#### define file names, software paths ####
work_dir <- "~/Desktop/dapc/" 
setwd(work_dir)
all.gen  = "ec.mac2.thin1000.7109.287.gen"
all.renamed.gen="all.renamed.gen"
species="ec"
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

for (i in 1:9) {
  SamplePop <- str_replace_all(SamplePop, as.character(i), "")
}

pop(all.genind) <- SamplePop
##subset populations for use from genind file##
PopNames.all
 
for (i in PopNames.all) {
  pop.sublist <- "CQCQ" #c("CQCQ", "GDGZ", "HNHK", "SCDY", "SCLS", "SXTY", "YNLJ", "YNYX") #
  pop.sublistname <- paste(species, pop.sublist, collapse ="_", sep=".")
  subset.genind <- popsub(all.genind, sublist = pop.sublist) # sublist=1:10, sublist=1:10, blacklist="Bazadais", sublist=c(1:6, 11:15)
  pop.sublistname

#############################################################################################
#### set the final used data in subsequent analysis ####
  used.genind = subset.genind
  ##convert genind to genlight format
  used.genlight <- gi2gl(used.genind)
#############################################################################################

## set color for figures and population names used
  PopNames.used <- popNames(used.genind)
  cols <- rainbow(nPop(used.genind))
#toRemove <- is.na(glMean(used.genlight, alleleAsUnit = T))
#used.genlight <- used.genlight[, !toRemove]

#### PCA from genelight file ####
  glpca <- glPca(used.genlight,parallel = TRUE, n.cores = 7L,  nf = 2)
  
  #pdf(file=paste(pop.sublistname, ".pca.genlight.barplot.pdf", sep=""))
  #barplot(100*glpca$eig/sum(glpca$eig), col = heat.colors(50), main="PCA Eigenvalues")
  #title(ylab="Percent of variance\nexplained", line = 2)
  #title(xlab="Eigenvalues", line = 1)
  #dev.off()

  glpca.scores <- as.data.frame(glpca$scores)
  glpca.scores$pop <- pop(used.genlight)
  glpca.scores$ind_id <- indNames(used.genind)
  #glpca.eig <- as.data.frame(glpca$eig)
  
  set.seed(9)
  pdf(file=paste(pop.sublistname, ".pca.genlight.scatter.ggplot.pdf", sep=""))
  p <- ggplot(glpca.scores, aes(x=PC1, y=PC2, colour=pop, label=ind_id))  ##add x, y
  p <- p + geom_point(size=1) ##add dot
  p <- p + stat_ellipse(level = 0.95, size = 0.4) ## add circle
  p <- p + scale_color_manual(values = cols) 
  p <- p + geom_hline(yintercept = 0) 
  p <- p + geom_vline(xintercept = 0) 
  p <- p + theme_bw() + geom_text(size=3)
  print(p)
  dev.off()
}
