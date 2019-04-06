##A recent special issue in Molecular Ecology Resources provides a nice overview of the arsenal of tools available in R (Paradis et al., 2017). New tools have become available in R for analyzing HTS data including adegenet (Jombart, 2008), ape (Paradis, Claude & Strimmer, 2004), vcfR (Knaus & Gr??nwald, 2017), and poppr (Kamvar, Tabima & Gr??nwald, 2014; Kamvar, Brooks & Gr??nwald, 2015) among others. 
##pipe in R: ">", see https://rforcats.net/#pipes

##install.packages(c("poppr", "mmod", "magrittr", "treemap"), repos = "http://cran.rstudio.com", dependencies = TRUE)
library("poppr")
library(vcfR)
library(RColorBrewer)
library(pinfsc50)
library(ape)
library(reshape2)
#setwd("~/Desktop/ML/ddrad/PS/filter")
getwd()

########################Part 1: import and read vcf file into R, check data#################
cs.vcf.raw <- vcfR::read.vcfR("gm.ddRAD.274.5X.0.999.15067snp.vcf", verbose = FALSE) #time-costing
cs.vcf <- cs.vcf.raw
#strwrap(cs.vcf@meta[1:10])
#queryMETA(cs.vcf)
#queryMETA(cs.vcf, element = 'DP')
#queryMETA(cs.vcf, element = 'FORMAT=<ID=DP')
#head(getFIX(cs.vcf))
cs.vcf@gt[1:10, 1:4]

#########################################Part 2: Quality Control############################
##check summary data, including missing data
cs.vcf

##extract depth of genotype#
dp <- extract.gt(cs.vcf, element = "DP", as.numeric=TRUE) #time-costing

##barplot of missingness for every sample#
myMiss <- apply(dp, MARGIN = 2, function(x){ sum(is.na(x)) }) ###data from bwa and bowtie2, miss data will represent with NA, while data from bowtie, miss data will be indecated by 0.
#myMiss <- apply(dp, MARGIN = 2, function(x){ sum(x==0) })   ####data from bowtie?

myMiss <- myMiss/nrow(cs.vcf)
pdf(file="gm_ddRAD5X_missing_barplot.pdf", width=30, height=8)
palette(brewer.pal(n=12, name = 'Set3'))
par(mar = c(16,4,4,2))
barplot(myMiss, las = 2, col = 1:16)
title(ylab = "Missingness (%)")
dev.off()

##histogram of missingness for every sample#
myMiss <- apply(dp, MARGIN = 1, function(x){ sum(is.na(x)) }) ###data from bwa, miss data will represent with NA, while data from bowtie, miss data will be indecated by 0.
#myMiss <- apply(dp, MARGIN = 1, function(x){ sum(x==0) })   ####data from bowtie
myMiss <- myMiss/ncol(cs.vcf@gt[,-1])
pdf(file="gm_ddRAD5X_missing_histogram.pdf", width=10, height=8)
hist(myMiss, col = "#8DD3C7", xlab = "Missingness (%)", main = "")
dev.off()

##draw barplot of depth for every sample#
pdf(file="gm_ddRAD5X_depth_boxplot.pdf", width=30, height=8)
par(mar=c(12,4,4,2))
boxplot(dp, col=2:8, las=3)
title(ylab = "Depth (DP)")
dev.off()

##draw violin plot of depth#
library(reshape2)
library(ggplot2) 
library(cowplot)

##Melt our matrix into a long form data.frame.
dpf <- melt(dp, varnames=c('Index', 'Sample'), value.name = 'Depth', na.rm=TRUE)
dpf <- dpf[ dpf$Depth > 0,]

##Create a row designator, You may want to adjust this
samps_per_row <- 25                               ####it is better to be divided evenly by numbel of individuals
myRows <- ceiling(length(levels(dpf$Sample))/samps_per_row)
myList <- vector(mode = "list", length = myRows)

for(i in 1:myRows){
  myIndex <- c(i*samps_per_row - samps_per_row + 1):c(i*samps_per_row)
  myIndex <- myIndex[myIndex <= length(levels(dpf$Sample))]
  myLevels <- levels(dpf$Sample)[myIndex]
  myRegex <- paste(myLevels, collapse = "$|^")
  myRegex <- paste("^", myRegex, "$", sep = "")
  myList[[i]] <- dpf[grep(myRegex, dpf$Sample),]
  myList[[i]]$Sample <- factor(myList[[i]]$Sample)
}

# Create the plot.
	pdf("gm_ddRAD5X_depth_violin.pdf", width = 10, height=8)
	myPlots <- vector(mode = "list", length = myRows)
	for(i in 1:myRows){
	  myPlots[[i]] <- ggplot(myList[[i]], aes(x=Sample, y=Depth)) + 
	    geom_violin(fill="#8dd3c7", adjust=1.0, scale = "count", trim=TRUE)
	  
	  myPlots[[i]] <- myPlots[[i]] + theme_bw()
	  myPlots[[i]] <- myPlots[[i]] + theme(axis.title.x = element_blank(), 
	                                       axis.text.x = element_text(angle = 60, hjust = 1))
	  myPlots[[i]] <- myPlots[[i]] + scale_y_continuous(trans=scales::log2_trans(), 
	                                                    breaks=c(1, 10, 100, 800),
	                                                    minor_breaks=c(1:10, 2:10*10, 2:8*100))
	  myPlots[[i]] <- myPlots[[i]] + theme( panel.grid.major.y=element_line(color = "#A9A9A9", size=0.6) )
	  myPlots[[i]] <- myPlots[[i]] + theme( panel.grid.minor.y=element_line(color = "#C0C0C0", size=0.2) )
	}
	myPlots
	dev.off()

##heatmap of depth#
pdf("gm_ddRAD5X_depth_raw_variant_heatmap.pdf", width = 40, height=40)
heatmap.bp((dp[1:1000,])^(1/4),rlabels = F)           #########some sequenced samples are so deep that it was hard to read by heatmap, so the depth was divided by the square root---Lijun
dev.off()

##ommit data, replace low and high coverage SNPs (DP value) with NA#
quants <- apply(dp, MARGIN=2, quantile, probs=c(0.1, 0.999), na.rm=TRUE) #probs=c(1, 2) is key to define a DP value range, can be adjusted.
quants[,1:11] #adjust probs in last setp according to the quantile coverages 
#dp2 <- sweep(dp, MARGIN=2, FUN = "-", quants[1,])
#dp[dp2 < 0] <- NA #remove the SNPs coverage lower than 10%, do not change???. if we us dp<4, this commond should be deactivated.
dp[dp < 8] <- NA #select SNPsat least 4X coverage, canbe adjusted.
dp2 <- sweep(dp, MARGIN=2, FUN = "-", quants[2,])
dp[dp2 > 0] <- NA #remove the SNPs coverage higher than 99.99%, do not change???

cs.vcf@gt[1:10, 1:4]
cs.vcf@gt[,-1][ is.na(dp) == TRUE ] <- NA    #use -1 to omit the first column of the gt matrix,which is the ???FORMAT??? column
cs.vcf@gt[1:10, 1:4]

pdf("gm_ddRAD_depth_ommit_data_heatmap.pdf", width = 40, height=40)
heatmap.bp((dp[1:1000,])^(1/4),rlabels = F) 
dev.off()

##ommiting samples#
myMiss <- apply(dp, MARGIN = 2, function(x){ sum( is.na(x) ) } ) 
myMiss <- myMiss / nrow(dp)


pdf("gm_ddRAD_depth_ommit_sample_barplot.pdf", width = 40, height=40)
barplot(myMiss, las = 3)
dev.off()

cs.vcf@gt <- cs.vcf@gt[, c(TRUE, myMiss < 0.2)] #myMiss is key to remove samples, determined according ps_ddRAD_depth_ommit_sample_barplot.tiff

cs.vcf

dp <- extract.gt(cs.vcf, element = "DP", as.numeric=TRUE)

pdf("gm_ddRAD_depth_ommit_sample_heatmap.pdf", width = 40, height=40)
heatmap.bp((dp[1:1000,])^(1/4), rlabels = FALSE)
dev.off()

##ommiting variants (SNP)#
myMiss <- apply(dp, MARGIN = 1, function(x){ sum( is.na(x) ) } ) 
myMiss <- myMiss / ncol(dp)
vcf_final <- cs.vcf[myMiss < 0.05, ] #change value to keep final number of SNPS.
vcf_final

dp <- extract.gt(vcf_final, element = "DP", as.numeric=TRUE)


pdf("gm_ddRAD_depth_ommit_variant_heatmap.pdf", width = 40, height=40)
heatmap.bp(dp^(1/4), rlabels = FALSE)
dev.off()

write.vcf(vcf_final, file = "gm_ddRAD_258_8X_0.999_8135snp.vcf.gz")

#######################################################
#####We can also ommit variants and then samples, if too many samples were deleted
####################################
cs.vcf <- cs.vcf.raw
##extract depth of genotype#
dp <- extract.gt(cs.vcf, element = "DP", as.numeric=TRUE) #time-costing

quants <- apply(dp, MARGIN=2, quantile, probs=c(0.1, 0.999), na.rm=TRUE) #probs=c(1, 2) is key to define a DP value range, can be adjusted.
dp[dp < 5] <- NA #select SNPsat least 4X coverage, canbe adjusted.
dp2 <- sweep(dp, MARGIN=2, FUN = "-", quants[2,])
dp[dp2 > 0] <- NA #remove the SNPs coverage higher than 99.99%, do not change???

cs.vcf@gt[1:10, 1:4]
cs.vcf@gt[,-1][ is.na(dp) == TRUE ] <- NA    #use -1 to omit the first column of the gt matrix,which is the ???FORMAT??? column
cs.vcf@gt[1:10, 1:4]

myMiss <- apply(dp, MARGIN = 1, function(x){ sum( is.na(x) ) } ) 
myMiss <- myMiss / ncol(dp)
cs.vcf  <- cs.vcf[myMiss < 0.1, ] #change value to keep final number of SNPS.
cs.vcf

dp <- extract.gt(cs.vcf, element = "DP", as.numeric=TRUE)
pdf("gm_ddRAD_depth_ommit_variant_heatmap1.pdf", width = 40, height=40)
heatmap.bp(dp^(1/4), rlabels = FALSE)
dev.off()

####ommited individuals
myMiss <- apply(dp, MARGIN = 2, function(x){ sum( is.na(x) ) } ) 
myMiss <- myMiss / nrow(dp)
cs.vcf@gt<- cs.vcf@gt[, c(TRUE, myMiss < 0.2)] #myMiss is key to remove samples, determined according ps_ddRAD_depth_ommit_sample_barplot.tiff
cs.vcf

dp <- extract.gt(cs.vcf, element = "DP", as.numeric=TRUE)

pdf("gm_ddRAD_depth_ommit_sample_heatmap1.pdf", width = 40, height=40)
heatmap.bp((dp[1:1000,])^(1/4), rlabels = FALSE)
dev.off()

cs.vcf

write.vcf(cs.vcf, file = "gm.ddRAD.274.5X.0.999.15067snp.vcf.gz")

####*****terminal: vcftools --gzvcf cs_ddRAD_203_10X_0.999_12668snp.vcf.gz  --maf 0.05 --recode --out cs_ddRAD_203_10X_0.999_12663_maf0.05
####*****terminal: vcftools --gzvcf cs_ddRAD_203_10X_0.999_12668snp.vcf.gz  --maf 0.1 --recode --out cs_ddRAD_203_10X_0.999_12663_maf0.1
