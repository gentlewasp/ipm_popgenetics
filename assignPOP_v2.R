################  establish refpops

library(assignPOP)
library(klaR)
library(mime)
library(genepopedit)

work_dir <- "/Users/macbook2017/Desktop/dapc/" 
setwd(work_dir)
all.gen  = "px_ddRAD_2018_432_6230_mac3_thin500-twogroups2.gen"
all.renamed.gen="all.renamed.gen"
species = "px_2018"
genepop_ID(genepop=all.gen, path=paste0(work_dir, all.renamed.gen)) 
PopNames.all <- genepop_detective(all.renamed.gen, variable="Pops")
PopCounts.all <- genepop_detective(all.renamed.gen, variable="PopNum")
PopCounts.all
pop.sublist <- c("GD", "YN")
pop.sublistname <- paste(pop.sublist, collapse ="_", sep=".")
subset_genepop(genepop= all.renamed.gen, keep = TRUE, 
                 spop = pop.sublist, 
                 path = paste0(work_dir, paste(pop.sublistname, ".pop.gen", sep="")))
SampleIDs <- genepop_detective(paste(pop.sublistname, ".pop.gen", sep=""), variable="Inds")
SampleIDs
subid <- c("GD_GA03", "GD_GA19", "GD_GA11", "GD_GA09", "GD_GA01", "GD_GA10", "GD_GA08", "GD_GA12", "GD_GA20", "GD_GA14", "GD_GA02", "GD_GA07","GD_GA13",
           "GD_SZ13", "GD_SZ18",  "GD_SZ20", "GD_SZ10", "GD_SZ06", "GD_SZ19", "GD_SZ14", 
           "GD_GB02", "GD_GC13", "GD_GC12", "GD_GC07", "GD_GC02",
           "YN_KM04", "YN_KM13", "YN_KM08", "YN_KM14", "YN_KM17", "YN_KM07", "YN_KM19", "YN_KM15",
           "YN_YX08", "YN_YX19", "YN_YX02", "YN_YX06", "YN_YX01", "YN_YX14") #"GXNN_01", "GXNN_06", "GDGC_01", "GDGC_13","YNYX_13", "YNYX_14", "YNKM_13", "YNKM_17"
subset_genepop_individual(genepop= paste(pop.sublistname, ".pop.gen", sep=""), 
                          indiv = subid, 
                          keep = FALSE,  
                          path = paste0(work_dir,paste(pop.sublistname, ".gen", sep="")))
used.SampleIDs <- genepop_detective(paste(pop.sublistname, ".gen", sep=""), variable="Inds")
used.SampleIDs
used.PopNames <- genepop_detective(paste(pop.sublistname, ".gen", sep=""), variable="Pops")
used.PopCounts <- genepop_detective(paste(pop.sublistname, ".gen", sep=""), variable="PopNum")
used.PopCounts
all_ref <- read.Genepop(paste(pop.sublistname, ".gen", sep=""), pop.names=used.PopNames)
#all_inc <- read.Genepop( "nokin_allq9_inc.txt")
#all_ref_rd <- reduce.allele(all_ref, p = 0.95)

## MCMC ##
model.set = "svm"  #"lda", "svm", "naiveBayes", "tree", and "randomForest". Default is "svm"(support vector machine).
assign_trial = assign.MC(all_ref, train.inds=c(0.92), train.loci=c(0.8), 
                         loci.sample="fst", iterations=50, 
                         dir = paste(species, pop.sublistname, model.set, "result/", sep="."), 
                         model = model.set )   
accuMC <- accuracy.MC(dir = paste(species, pop.sublistname, model.set, "result/", sep="."))
pdf(file=paste(species, pop.sublistname, model.set, "assignpop.plot.pdf", sep="."))
accuracy.plot(accuMC, pop = c("all", used.PopNames))
dev.off()

## Kfold validation ##
assign.kfold( all_ref, k.fold=c(2,3), train.loci=c(1), 
              loci.sample="random", model=model.set, 
              dir = paste(species, pop.sublistname, "kfold_result/", sep="."))
accuKF <- accuracy.kfold(dir = paste(species, pop.sublistname, "kfold_result/", sep=".")) #Use this function for K-fold cross-validation results
pdf(file=paste(species, pop.sublistname, "kfold.membership.plot.pdf", sep="."), width=30, height=8)
membership.plot(dir = paste(species, pop.sublistname, "kfold_result/", sep=".")) 
dev.off()

check.loci(dir = paste(species, pop.sublistname, model.set, "result/", sep="."), top.loci = 1000)

#used.PopNames
file.remove(paste(pop.sublistname, ".gen", sep=""))
#file.remove(paste(pop.sublistname, ".pop.gen", sep=""))
#file.remove("all.renamed.gen")

##assigntest##
PopNames.all
pop.sublist.test <- c("GDHK", "GDG",  "GDGC", "GDNN" , "FJFZ", "YNDH", "YNPZ", "SCPZ", "HNCS", "JXNA", "JXNB", "JXNC", "ZJNB", "SCCD", "SCGY", "NMAL", "BJYA", "HBZA")
pop.sublistname.test <- paste(pop.sublist.test, collapse ="_", sep=".")
subset_genepop(genepop= all.renamed.gen, keep = TRUE, 
               spop = pop.sublist.test, 
               path = paste0(work_dir, paste(pop.sublistname.test, ".pop.gen", sep="")))
used.PopNames.test <- genepop_detective(paste(pop.sublistname.test, ".pop.gen", sep=""), variable="Pops")
all_inc <- read.Genepop(paste(pop.sublistname.test, ".pop.gen", sep=""), pop.names=used.PopNames.test)
pdf(file=paste(pop.sublistname, "(ref)", pop.sublistname.test, model.set, "assignpop.plot.pdf",sep="."))
assign.X(x1=all_ref, x2=all_inc, 
         dir=paste(pop.sublistname, "(ref)", pop.sublistname.test, model.set, "result/", sep="."), 
         model = model.set, mplot=TRUE)
dev.off()
