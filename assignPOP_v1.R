################  establish refpops

library(assignPOP)
library(klaR)
library(mime)
library(genepopedit)

work_dir <- "/Users/macbook2017/Desktop/dapc/" 
setwd(work_dir)
all.gen  = "ec.140.2.500.6437.neutral_3pop.gen"
all.renamed.gen="all.renamed.gen"

genepop_ID(genepop=all.gen, path=paste0(work_dir, all.renamed.gen))  ##BJYQ01 > BJYQ_01

PopNames.all <- genepop_detective(all.renamed.gen, variable="Pops")
PopCounts.all <- genepop_detective(all.renamed.gen, variable="PopNum")
PopCounts.all
all_ref <- read.Genepop(all.renamed.gen, pop.names=PopNames.all)

#all_inc <- read.Genepop( "nokin_allq9_inc.txt")
#all_ref_rd <- reduce.allele(all_ref, p = 0.95)

assign_trial = assign.MC(all_ref, train.inds=c(0.9), train.loci=c(0.8), 
                         loci.sample="fst", iterations=100, dir="mc_randomForest1/", model="randomForest" )   

accuMC200 <- accuracy.MC(dir = "mc_randomForest1/")
accuracy.plot(accuMC200, pop = PopNames.all) 

