################  establish refpops

library(assignPOP)
library(klaR)
library(mime)
library(genepopedit)

work_dir <- "/Users/macbook2017/Desktop/dapc/" 
setwd(work_dir)
all.gen  = "px_ddRAD_2018_432_6230_mac3_thin500-twogroups.gen"
all.renamed.gen="all.renamed.gen"
species = "px_2018."
genepop_ID(genepop=all.gen, path=paste0(work_dir, all.renamed.gen)) 
PopNames.all <- genepop_detective(all.renamed.gen, variable="Pops")
PopCounts.all <- genepop_detective(all.renamed.gen, variable="PopNum")
PopNames.all
PopCounts.all
pop.sublist <- c("GD", "YN")
pop.sublistname <- paste(pop.sublist, collapse ="_", sep=".")
subset_genepop(genepop= all.renamed.gen, keep = TRUE, 
                 spop = pop.sublist, 
                 path = paste0(work_dir, paste(pop.sublistname, ".pop.gen", sep="")))
SampleIDs <- genepop_detective(paste(pop.sublistname, ".pop.gen", sep=""), variable="Inds")
SampleIDs
subid <- c() #"GXNN_01", "GXNN_06", "GDGC_01", "GDGC_13","YNYX_13", "YNYX_14", "YNKM_13", "YNKM_17"
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
assign_trial = assign.MC(all_ref, train.inds=c(0.92), train.loci=c(0.8), 
                         loci.sample="fst", iterations=50, 
                         dir = paste(species, pop.sublistname, ".mc_svm/", sep=""), 
                         model="svm" )   
accuMC <- accuracy.MC(dir = paste(species, pop.sublistname, ".mc_svm/", sep=""))
pdf(file=paste(species, pop.sublistname, ".assignpop.mc.plot.pdf", sep=""))
accuracy.plot(accuMC, pop = c("all", used.PopNames))
dev.off()
#used.PopNames
file.remove(paste(pop.sublistname, ".gen", sep=""))
#file.remove(paste(pop.sublistname, ".pop.gen", sep=""))
#file.remove("all.renamed.gen")
PopNames.all
pop.sublist.test <- c("FJFZ", "SCPZ", "HNCS", "JXNA", "NMAL", "BJYA", "HBZA")
pop.sublistname.test <- paste(pop.sublist.test, collapse ="_", sep=".")
subset_genepop(genepop= all.renamed.gen, keep = TRUE, 
               spop = pop.sublist.test, 
               path = paste0(work_dir, paste(pop.sublistname.test, ".pop.gen", sep="")))
used.PopNames.test <- genepop_detective(paste(pop.sublistname.test, ".pop.gen", sep=""), variable="Pops")
all_inc <- read.Genepop(paste(pop.sublistname.test, ".pop.gen", sep=""), pop.names=used.PopNames.test)
pdf(file=paste(pop.sublistname, "(ref)", pop.sublistname.test, ".assignpop.plot.pdf",sep="."))
assign.X(x1=all_ref, x2=all_inc, dir=paste(pop.sublistname, "(ref)", pop.sublistname.test, "svm/",sep="."), model="svm",mplot=TRUE)
dev.off()
