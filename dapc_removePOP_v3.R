##https://grunwaldlab.github.io/Population_Genetics_in_R/DAPC.html
#source("https://bioconductor.org/biocLite.R")
#biocLite("qvalue")
#install.packages("dartR")

library("poppr")
library(vcfR)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(cowplot)
library(adegenet)
#library(snowfall)
library(dartR)
library("radiator")
library("genepopedit")
library(stringr)

##define file names, software paths
output_dir <- "/Users/macbook2017/Desktop/dapc/" 
oldgenepop=   "px_ddRAD_2018_432_5074_all_4X_0.999_thin_genepop.txt"

setwd(output_dir)
removepop.genepop="removePOP.gen"
removepop_ind.genepop="removepop_ind.genepop.gen"
plinkpath="/Users/macbook2017/Desktop/softwares/plink_mac_20190304"
pgdspiderpath="/Users/macbook2017/Desktop/softwares/PGDSpider_2.1.1.5"
genepop_ID(genepop=oldgenepop, path=paste0(output_dir, oldgenepop)) ##change genepop individual name to detect pop nane

##get variable from genepop using genepop_detective
#PopNames <- genepop_detective(oldgenepop, variable="Pops")
#PopCounts <- genepop_detective(oldgenepop, variable="PopNum")
#SampleIDs <- genepop_detective(oldgenepop, variable="Inds")
#LociNames <- genepop_detective(oldgenepop, variable="Loci")
#metadata <- genepop_detective(oldgenepop,variable="All")
#Alleles <- genepop_detective(oldgenepop,variable="Allele")

## set population to be removed
PopNames <- genepop_detective(oldgenepop, variable="Pops") # check original population
PopNames
PopKeep <- c("HNHK", "GDGZ", "YNKM", "SCLS", "CQCQ", "SCYA", "SCNC", "SCCD", "SCGY") # to be kept populations
#PopKeep <- setdiff(PopNames, c("CCC","GGG")) # to be remvoed populations
## set loci to be removed
#LociNames <- genepop_detective(removepop.genepop, variable="Loci")
#subloci <- setdiff(LociNames,c("691:6:+"	,	"1114:19:+"	,	"2456:7:-" , "2774:5:+"	,	"3216:7:-" , "3342:6:-"	, "4138:7:-" , "4327:9:-"	,	"4674:7:-" , "5184:8:+"	,	"5875:5:-" , "6079:5:-"	, "6210:6:-" , "6834:5:+"	,	"6996:6:-" , "7281:9:-"	,	"7453:7:-" , "7846:6:+"	, "8844:20:-"	,	"9061:29:+"	,	"9202:118:+" , "9385:39:-" , "9988:5:-"	,	"11167:7:-"	, "12628:266:+", "12742:73:+"	,	"13157:5:-"	,	"13341:6:+"	,	"13579:7:-"	,	"16112:14:+" , "16279:5:+"	,	"16591:13:-" , "17762:5:-" , "18171:10:+"	,	"18316:16:-" , "18966:6:-" , "19158:101:+", "20683:8:-" , "21647:11:-"	,	"21895:28:+" , "22478:11:+"	,	"22685:5:-"	, "23396:6:+"	,	"25047:6:-"	,	"25429:15:+" , "25506:280:+",	"26940:78:+" , "27253:5:-" , "28157:7:+"	,	"28626:18:+" , "29342:6:-" , "29874:5:-" , "30035:11:+"	,	"30776:7:+"	, "31891:10:+" , "32432:10:-"	,	"32561:11:+" , "32763:5:-" , "32946:7:-" , "33097:11:-"	, "33445:30:-" , "34028:9:-" , "34325:5:+" , "34835:41:-"	,	"35111:15:+" , "35281:10:-"	, "36231:5:+"	,	"36487:5:-"	,	"37945:14:+" , "38432:5:+" , "38793:6:+" , "39139:9:+" , "40232:5:-"	,	"42708:132:+", "42876:7:-"))
#remove populations
subset_genepop(genepop= oldgenepop, keep = TRUE, 
               spop = PopKeep, 
               #subs = subloci,
               path = paste0(output_dir,removepop.genepop))  

## set individuals to be removed, and remove individuals
#SampleIDs <- genepop_detective(removepop.genepop, variable="Inds")
#SampleIDs
subid <- c("GXNN_02","GXNN_09","GXNN_05","GDSZ_13","GDGZ_13","HNHK_13","HNHK_13","GDGC_13","GDGB_14", "YNYX_14", "GXNY_14", "YNDH_14")
subset_genepop_individual(genepop= removepop.genepop, 
                          indiv = subid, 
                          keep = FALSE,  
                          path = paste0(output_dir,removepop_ind.genepop))

##read in removed populations genepop file
genindData <- read.genepop(removepop_ind.genepop, ncode=3)
##set population to geneind file
SampleIDs <- genepop_detective(removepop_ind.genepop, variable="Inds")
SampleIDs
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
SamplePop
pop(genindData) <- SamplePop
##convert genind to genlight format
genlight <- gi2gl(genindData)
#toRemove <- is.na(glMean(genlight, alleleAsUnit = T))
#genlight <- genlight[, !toRemove]

#####  find optimal number of cluster #####
cluster <- find.clusters(genlight, n.clust=NULL,                                                                        
                         max.n.clust=100,                                                                                      
                         stat=c("BIC"),                                                                                       
                         n.iter=1e9, n.start=1e3, # 1e9, 1e3                                                                  
                         #truenames=TRUE,
                         scale=FALSE,
                         parallel=TRUE)
PC=100                #number of principle componetskept, for dapc analysis
numberofcluster = 2   #number of clusters kept, for dapc analysis

#####DAPC#####
dapc <- dapc(genlight, n.da=numberofcluster, n.pca=PC)

pdf(file="dapc.assignscatter_pca50.pdf")                                                                                      
scatter(dapc,cstar=0,                                                                                                   
        mstree=TRUE,                                                                                                          
        posi.da="bottomright", posi.pca="bottomleft", scree.pca=TRUE,                                                         
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
pramx <- xvalDapc(tab(genlight, NA.method = "mean"), pop(genlight))

set.seed(9999)
system.time(pramx <- xvalDapc(tab(genlight, NA.method = "mean"), pop(genlight),
                              n.pca = 50:70, n.rep = 30,
                              parallel = "multicore", ncpus = 7L))

names(pramx)
pramx[-1]

scatter(dapc, col=seasun(14), cex = 2, legend = TRUE,
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1,leg=TRUE,)


mycol= rainbow(numberofcluster)
pdf(file="compoplot.cluters.pdf")
compoplot(dapc,posi="bottomright",
          txt.leg=paste("Cluster", 1:numberofcluster),
          ncol=2,  col=mycol, xlab="individuals", 
          lab=FALSE)
dev.off()


pdf(file="dapc.cluster.table.pdf")                                                                                            
table.value(table(pop(genlight), cluster$grp), 
            col.lab=paste("Cluster", 1:numberofcluster),
            row.lab=c("HNHK", "GDGZ", "HNCS", "JXNC", "SCLS", "CQCQ", "SCYA", "SCNC", "SCCD", "SCGY", "HBXG", "HBWH", "AHHF", "BJYQ", "HBZJ")) 
dev.off()
