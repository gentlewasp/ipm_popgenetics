################  establish refpops

library(assignPOP)
library(klaR)
library(mime)
library(genepopedit)

work_dir <- "/Users/macbook2017/Desktop/dapc/" 
setwd(work_dir)
all.gen  = "ec.287.2.500.6347.remove.outliers.gen.txt"
all.renamed.gen="all.renamed.gen"
species = "ec"
genepop_ID(genepop=all.gen, path=paste0(work_dir, all.renamed.gen)) 
PopNames.all <- genepop_detective(all.renamed.gen, variable="Pops")
PopCounts.all <- genepop_detective(all.renamed.gen, variable="PopNum")
PopCounts.all
for (a in PopNames.all) {
  for (b in PopNames.all) {
    if(a<b){
      pop.sublist <- c(a,b)
    }
  }
  pop.sublistname <- paste(pop.sublist, collapse ="_", sep=".")
  subset_genepop(genepop= all.renamed.gen, keep = TRUE, 
                 spop = pop.sublist, 
                 path = paste0(work_dir, paste(pop.sublistname, ".gen", sep="")))
  used.PopNames <- genepop_detective(paste(pop.sublistname, ".gen", sep=""), variable="Pops")
  all_ref <- read.Genepop(paste(pop.sublistname, ".gen", sep=""), pop.names=used.PopNames)
  #all_inc <- read.Genepop( "nokin_allq9_inc.txt")
  #all_ref_rd <- reduce.allele(all_ref, p = 0.95)
  assign_trial = assign.MC(all_ref, train.inds=c(0.92), train.loci=c(0.8), 
                         loci.sample="fst", iterations=200, dir = paste(species, pop.sublistname, "mc_svm/", sep="."), model="svm" )   
  accuMC <- accuracy.MC(dir = paste(species, pop.sublistname, "mc_svm/", sep="."))
  pdf(file=paste(species, pop.sublistname, "assignpop.mc.plot.pdf", sep="."))
  plot.result <- accuracy.plot(accuMC, pop = c("all", used.PopNames))
  print(plot.result)
  dev.off()
  }
file.remove("all.renamed.gen")

