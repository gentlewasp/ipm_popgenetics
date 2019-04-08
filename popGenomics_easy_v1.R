##https://grunwaldlab.github.io/Population_Genetics_in_R/DAPC.html
##https://github.com/grunwaldlab/Population_Genetics_in_R/blob/master/gbs_analysis.Rmd#principal-components-analysis
#source("https://bioconductor.org/biocLite.R")
#biocLite("qvalue")
#install.packages("dartR")

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
work_dir <- "/Users/macbook2017/Desktop/dapc/" 
all.gen  = "px_ddRAD_2011_184_728_all_4X_0.999_thin_genepop.txt"
#############################################################################################

setwd(work_dir)
plinkpath="/Users/macbook2017/Desktop/softwares/plink_mac_20190304"
pgdspiderpath="/Users/macbook2017/Desktop/softwares/PGDSpider_2.1.1.5"
all.renamed.gen="all.renamed.gen"
#all.vcf <- read.vcfR("ec.ddrad.mac2.thin1000.7109snps.287samples.recode.vcf") ##import all populations vcf
#removePOP.gen="removePOP.gen"
#removePOP_IND.gen="removePOP_IND.gen.gen"

#############################################################################################
#### data manipulation, reformat genepop individual names to get population names automatically ####
genepop_ID(genepop=all.gen, path=paste0(work_dir, all.renamed.gen)) 

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
used.genind = all.genind
#############################################################################################

## set color for figures and population names used
#cols <- brewer.pal(n = nPop(used.genind), name = "Dark2")  ## maximum is 8, if more than 8 populatins, changed it manually
PopNames.used <- popNames(used.genind)
cols <- rainbow(nPop(used.genind))
##convert genind to genlight format
used.genlight <- gi2gl(used.genind)
toRemove <- is.na(glMean(used.genlight, alleleAsUnit = T))
used.genlight <- used.genlight[, !toRemove]

#### remove population, individual, loci using radiator pacakge ####
#library("radiator")
## set population to be removed
#PopNames.all <- genepop_detective(all.renamed.gen, variable="Pops") # check original population
#PopNames.all
#popnumber.all <- length(PopNames.all)+1
#PopKeep <- c("CQCQ", "GDGZ", "HNHK", "SCDY", "SCLS", "SXTY", "YNLJ", "YNYX") # to be kept populations
#PopKeep <- setdiff(PopNames.all, c("CCC","GGG")) # to be remvoed populations
## set loci to be removed
#LociNames <- genepop_detective(removePOP.gen, variable="Loci")
#subloci <- setdiff(LociNames,c("691:6:+"	,	"1114:19:+"	,	"2456:7:-" , "2774:5:+"	,	"3216:7:-" , "3342:6:-"	, "4138:7:-" , "4327:9:-"	,	"4674:7:-" , "5184:8:+"	,	"5875:5:-" , "6079:5:-"	, "6210:6:-" , "6834:5:+"	,	"6996:6:-" , "7281:9:-"	,	"7453:7:-" , "7846:6:+"	, "8844:20:-"	,	"9061:29:+"	,	"9202:118:+" , "9385:39:-" , "9988:5:-"	,	"11167:7:-"	, "12628:266:+", "12742:73:+"	,	"13157:5:-"	,	"13341:6:+"	,	"13579:7:-"	,	"16112:14:+" , "16279:5:+"	,	"16591:13:-" , "17762:5:-" , "18171:10:+"	,	"18316:16:-" , "18966:6:-" , "19158:101:+", "20683:8:-" , "21647:11:-"	,	"21895:28:+" , "22478:11:+"	,	"22685:5:-"	, "23396:6:+"	,	"25047:6:-"	,	"25429:15:+" , "25506:280:+",	"26940:78:+" , "27253:5:-" , "28157:7:+"	,	"28626:18:+" , "29342:6:-" , "29874:5:-" , "30035:11:+"	,	"30776:7:+"	, "31891:10:+" , "32432:10:-"	,	"32561:11:+" , "32763:5:-" , "32946:7:-" , "33097:11:-"	, "33445:30:-" , "34028:9:-" , "34325:5:+" , "34835:41:-"	,	"35111:15:+" , "35281:10:-"	, "36231:5:+"	,	"36487:5:-"	,	"37945:14:+" , "38432:5:+" , "38793:6:+" , "39139:9:+" , "40232:5:-"	,	"42708:132:+", "42876:7:-"))
#remove populations
#subset_genepop(genepop= all.renamed.gen, keep = TRUE,spop = PopKeep, path = paste0(work_dir,removePOP.gen))  #subs = subloci,
## set individuals to be removed, and remove individuals
#SampleIDs <- genepop_detective(removePOP.gen, variable="Inds")
#SampleIDs
#subid <- c("GXNN_02","GXNN_09","GXNN_05","GDSZ_13","GDGZ_13","HNHK_13","HNHK_13","GDGC_13","GDGB_14", "YNYX_14", "GXNY_14", "YNDH_14")
#subset_genepop_individual(genepop= removePOP.gen,indiv = subid, keep = FALSE,path = paste0(work_dir,removePOP_IND.gen))

#############################################################################################
#### basic statistics ####
poppr(used.genind, plot = TRUE, legend = TRUE)  #sublist=c("specific population to be analyzed"),
div <- summary(used.genind)

pdf(file = paste(pop.sublistname, ".Ho.perLocus.pdf", sep=""))
plot(div$Hobs, xlab="Loci number", ylab="Observed Heterozygosity", 
     main="Observed heterozygosity per locus")
dev.off()

pdf(file=paste(pop.sublistname, ".HeHo.pdf", sep=""))
plot(div$Hobs,div$Hexp, xlab="Hobs", ylab="Hexp", 
     main="Expected heterozygosity as a function of observed heterozygosity per locus")
dev.off()

bartlett.test(list(div$Hexp, div$Hobs))  # result p-value, difference between expected and observed heterozygosity.

used.hierfstat <- genind2hierfstat(used.genind) # convert to hierfstat format

basicstat <- basic.stats(used.hierfstat, diploid = TRUE, digits = 2) 
names(basicstat) 
basicstat

boot.ppfis(used.hierfstat) 

#hw.test(used.genind, B = 1000) ## test HWE

basic.stats(used.genind)  # Fst following Nei (1987) on genind object
wc(used.genind) # Weir and Cockerham's estimate

##Hierarchical Fst tests (=AMOVA for SNP dataset)
loci <- used.genind[, -1] # Remove the population column
varcomp.glob(levels = data.frame(population, county), loci, diploid = TRUE) 
test.g(loci, level = population)
test.between(loci, test.lev = population, rand.unit = county, nperm = 100) 

##Pairwise Fst
#genet.dist(used.genind, method = "WC84")  #slow

## individual-based genetic distance
distgenEUCL <- dist(used.genind, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
pdf(file=paste(pop.sublistname, ".alle.freqBased.indiv.genetic.distance.euclidean.pdf", sep=""))
hist(distgenEUCL)
dev.off()

#The option pairwise.deletion = FALSE in the command dist.gene() removes all loci with one missing values : you an see on the histogram that we get a maximum distance of 3 loci out of 100.
#We can see that we get 98 loci with at least one sample missing. Then using the option pairwise.deletion = TRUE in the command dist.gene() allows you to keep loci with one missing value.
used.loci <-genind2loci(used.genind)
distgenDIFF <- dist.gene(used.loci, method="pairwise", pairwise.deletion = FALSE, variance = FALSE)
pdf(file=paste(pop.sublistname, ".number.of.loci.for.which.individuals.differpairwise.deletionFALSE.pdf", sep=""))
hist(distgenDIFF)
dev.off()

distgenDIFF <- dist.gene(used.loci, method="pairwise", pairwise.deletion = TRUE, variance = FALSE)
pdf(file=paste(pop.sublistname, ".number.of.loci.for.which.individuals.differpairwise.deletionTRUE.pdf", sep=""))
hist(distgenDIFF)
dev.off()

distgenDISS <- diss.dist(used.genind, percent = FALSE, mat = FALSE)
pdf(file=paste(pop.sublistname, ".number.of.allelic.differences.between.two.individuals.pdf", sep=""))
hist(distgenDISS)
dev.off()

# Get percent missing data per population
missing_data <- info_table(used.genind, type = "missing")
sum(missing_data["Total", 1:1000] > 0)
pdf(file=paste(pop.sublistname, ".missing.data.barplot.pdf", sep=""))
barplot(missing_data["Total", 1:100], xlab = "Locus", ylab = "Percent Missing")
dev.off()

##summary Conclusions drawn from the analysis
pdf(file=paste(pop.sublistname, ".distgenEUCL.distgenDIFF.distgenDISS.pdf", sep=""))
boxplot(distgenEUCL, distgenDIFF, distgenDISS)
dev.off()

#############################################################################################
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

pdf(file=paste(pop.sublistname, ".pca.genlight.scatter.pdf", sep=""))
scatter(glpca,posi="bottomright")
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

##plot genelight file
pdf(file=paste(pop.sublistname, ".genlight.plot.pdf", sep=""))
glPlot(used.genlight,posi="bottomleft")
dev.off()

##phylogenetic tree
tree <- aboot(used.genlight, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T) # nj
pdf(file=paste(pop.sublistname, ".indiv.phyloTree.upgma.pdf", sep=""))
plot.phylo(tree,  cex = 0.8, font = 2, adj = 0, tip.color =  cols[pop(used.genlight)]) # type ="unrooted",
nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.6, font = 3, xpd = TRUE)
#legend(35,10,PopNames.used,cols, border = FALSE, bty = "n")
legend('topleft', legend =PopNames.used, fill = cols, border = FALSE, bty = "n", cex = 1)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")
dev.off()

##minimum spanning network
library(igraph)
dist <- bitwise.dist(used.genlight)
msn <- poppr.msn(used.genlight, dist, showplot = FALSE, include.ties = T)
node.size <- rep(2, times = nInd(used.genlight))
names(node.size) <- indNames(used.genlight)
vertex.attributes(msn$graph)$size <- node.size

set.seed(9)
pdf(file=paste(pop.sublistname, ".minimum.spanning.network.pdf", sep=""))
plot_poppr_msn(used.genlight, msn , palette = rainbow(nPop(used.genlight)), gadj = 70)
dev.off()

#####  find optimal number of cluster #####
#pdf(file=paste(pop.sublistname, ".dapc.findClusters.pdf", sep=""))
cluster <- find.clusters(used.genlight, n.clust=NULL,                                                                        
                         max.n.clust=30,                                                                                      
                         stat=c("AIC"),                                                                                       
                         n.iter=1e9, n.start=1e3, # 1e9, 1e3                                                                  
                         #truenames=TRUE,
                         scale=FALSE,
                         parallel=TRUE)
#dev.off()

PC=2                #number of principle componetskept, for dapc analysis
numberofcluster = 4   #number of clusters kept, for dapc analysis

#####DAPC#####
dapc <- dapc(used.genlight, n.da=numberofcluster, n.pca=PC)

pdf(file=paste(pop.sublistname, ".dapc.assignscatter.pdf", sep=""))
scatter(dapc,cstar=0,                                                                                                   
        #mstree=TRUE,                                                                                                          
        #posi.da= "bottomleft", 
        posi.pca = "topleft", scree.pca=TRUE,                                                         
        scree.da=FALSE,                                                                                                      
        pch=20,                                                                                                               
        leg=TRUE,                                                                                                             
        col=seasun(14),                                                                                                       
        clab=0.8, # population names in the circle                                                                              
        ratio.pca=0.3, solid=.6, cex=3)                                                                                       
dev.off() 

set.seed(4)
contrib <- loadingplot(dapc$var.contr, axis = 2, thres = 0.07, lab.jitter = 1)

##corss-valiation
set.seed(999)
pdf(file=paste(pop.sublistname, ".dapc.corss.valid.pdf", sep=""))
pramx <- xvalDapc(tab(used.genlight, NA.method = "mean"), pop(used.genlight))
dev.off()

set.seed(9999)
pdf(file=paste(pop.sublistname, ".dapc.corss.valid2.pdf", sep=""))
pramx <- xvalDapc(tab(used.genlight, NA.method = "mean"), pop(used.genlight),
                              n.pca = 90:110, n.rep = 100,
                              parallel = "multicore", ncpus = 7L)
dev.off()

names(pramx)
pramx[-1]

pdf(file=paste(pop.sublistname, ".dapc.assignscatter2.pdf", sep=""))
scatter(dapc, col=seasun(14), cex = 2, legend = TRUE,
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)
dev.off()


pdf(file=paste(pop.sublistname, ".dapc.compoplot.pdf", sep=""))
compoplot(dapc,col = cols, posi = 'top')
dev.off()

mycol= rainbow(numberofcluster)
pdf(file=paste(pop.sublistname, ".dapc.compoplot2.pdf", sep=""))
compoplot(dapc,posi="bottomright",
          txt.leg=paste("Cluster", 1:numberofcluster),
          col=mycol, xlab="individuals", 
          lab=FALSE)
dev.off()

pdf(file=paste(pop.sublistname, ".dapc.barplot.pdf", sep=""))
barplot(t(dapc$posterior), col = cols, las = 3, space = 0, border = NA, cex.names = 0.6)
dev.off()

pdf(file=paste(pop.sublistname, ".dapc.cluster.table.pdf", sep=""))
table.value(table(pop(used.genlight), cluster$grp), 
            col.lab=paste("Cluster", 1:numberofcluster),
            row.lab=PopNames.used) 
dev.off()

## draw dapc assignment using ggplot,This bar plot shows us a more organized perspective of our data set by contrasting the population membership probability assignments against their original populations.
dapc.results <- as.data.frame(dapc$posterior)
dapc.results$pop <- pop(used.genind)
dapc.results$indNames <- rownames(dapc.results)
library(reshape2)
dapc.results <- melt(dapc.results)

pdf(file=paste(pop.sublistname, ".dapc.compoplot3.pdf", sep=""))
coul = brewer.pal(4, "Dark2") 
colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")
p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = colorRampPalette(coul)(15))  ##values = cols
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p
dev.off()

