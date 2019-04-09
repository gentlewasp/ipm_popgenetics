##https://grunwaldlab.github.io/Population_Genetics_in_R/DAPC.html
##https://github.com/grunwaldlab/Population_Genetics_in_R/blob/master/gbs_analysis.Rmd#principal-components-analysis
#source("https://bioconductor.org/biocLite.R")
#biocLite("qvalue")
#install.packages("dartR")
devtools::install_github(c("thibautjombart/adegenet", "grunwaldlab/poppr"))

library(poppr)
library(vcfR)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(cowplot)
library(adegenet)
#library(snowfall)
library(dartR)
library("genepopedit")
library(stringr)
library(hierfstat)
library(pegas)
library(ape)

#############################################################################################
#### define file names, software paths ####
#############################################################################################
work_dir <- "/Users/macbook2017/Desktop/dapc/ps/" 
setwd(work_dir)
all.gen  = "ps.11pop.neutral.17517.gen"
#############################################################################################

#plinkpath="/Users/macbook2017/Desktop/softwares/plink_mac_20190304"
#pgdspiderpath="/Users/macbook2017/Desktop/softwares/PGDSpider_2.1.1.5"
all.renamed.gen="all.renamed.gen"
#all.vcf <- read.vcfR("ec.ddrad.mac2.thin1000.7109snps.287samples.recode.vcf") ##import all populations vcf
#removePOP.gen="removePOP.gen"
#removePOP_IND.gen="removePOP_IND.gen.gen"

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
all.genind <- genepop.to.genind(all.renamed.gen, ncode=3)

all.genind <- read.genepop("ps.11pop.neutral.17517.gen", ncode=3)
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
pop.sublist <- "AHHF" #"AHHF" "FJPT" "GDGZ" "GXNN" "HNCS" "HNYY" "JSWX" "JXGZ" "JXNC" "XJTL" "ZJHZ"
pop.sublistname <- paste(pop.sublist, collapse ="_")
subset.genind = paste(pop.sublistname, ".genind", sep = "") 
subset.genind <- popsub(all.genind, sublist = pop.sublist) # sublist=1:10, sublist=1:10, blacklist="Bazadais", sublist=c(1:6, 11:15)

used.genind = subset.genind
used.genlight <- gi2gl(used.genind)

library(HWxtest)
used.df <- genind2df(used.genind)

hwdf(hwx.test(all.genind),options(mc.cores = 8)) #
used.df
used.genind

hwe <- hw.test(used.genind,res.type="matrix") ## test HWE
hwe$fca90$P10
