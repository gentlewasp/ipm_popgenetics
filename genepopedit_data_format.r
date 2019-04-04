
##use genepopedit to convert snp data format
install.packages("devtools")
library(devtools)
install_github("rystanley/genepopedit")

#Use genepopedit to subset a SNP dataset by:
#    removing specified loci.
#    removing specified populations.
#    removing specified individuals.
#    renaming populations.
#    grouping populations.
#    reordering populations.
#    extracting meta-data:
#       population names.
#        population counts.
#        sample IDs.
#        loci names.
#        allele calls.
#        allele frequencies by population grouping.
#        loci linkage and global Weir and Cockerham's Fst.
#       lists of unlinked loci maximizing pairwise global Weir and Cockerham's Fst.
#    create datasets for training, assignment, and outlier detection, according to a population stratified random sample.
#    convert Genepop to STRUCTURE, FSTAT, NEWHYBRIDS, ASSIGNER, BGC, TREEMIX, COLONY, Genetic Stock Identification (gsi_sim), HZAR, and flattened/unflattened format.
##simulate individual genotypes using pooled DNA allele frequencies.

setwd("~/Dropbox/dbm_migration_ddRAD/data/genpop")

### import vcf file###
vcf <- vcfR::read.vcfR("px_ddRAD_2018_432_5074_all_4X_0.999_thin.vcf", verbose = FALSE) 
####converge to genlight####
genlight <- vcfR2genlight(vcf)
####define population of each sample####
pop(genlight) <- as.factor(c("HNHK","HNHK","HNHK","HNHK","HNHK","HNHK","HNHK","HNHK","HNHK","HNHK","HNHK","GDSZ","GDSZ","GDSZ","GDSZ","GDSZ","GDSZ","GDSZ","GDSZ","GDSZ","GDSZ","GDSZ","GDSZ","GDSZ","GDSZ","GDSZ","GDSZ","GDSZ","GDSZ","GXNN","GXNN","GXNN","GXNN","GXNN","GXNN","GXNN","GXNN","GXNN","GXNN","GXNN","GXNN","GXNN","GXNN","GXNN","GXNN","GXNN","GXNN","GXNN","GXNN","GDGA","GDGA","GDGA","GDGA","GDGA","GDGA","GDGA","GDGA","GDGA","GDGA","GDGA","GDGA","GDGA","GDGA","GDGA","GDGA","GDGA","GDGA","GDGA","GDGA","GDGB","GDGB","GDGB","GDGB","GDGB","GDGB","GDGB","GDGB","GDGB","GDGB","GDGB","GDGB","GDGB","GDGB","GDGB","GDGB","GDGB","GDGB","GDGC","GDGC","GDGC","GDGC","GDGC","GDGC","GDGC","GDGC","GDGC","GDGC","GDGC","GDGC","GDGC","GDGC","GDGC","GDGC","GDGC","YNYX","YNYX","YNYX","YNYX","YNYX","YNYX","YNYX","YNYX","YNYX","YNYX","YNYX","YNYX","YNYX","YNYX","YNYX","YNDH","YNDH","YNDH","YNDH","YNDH","YNDH","YNDH","YNDH","YNDH","YNDH","YNDH","YNDH","YNDH","YNDH","YNKM","YNKM","YNKM","YNKM","YNKM","YNKM","YNKM","YNKM","YNKM","YNKM","YNKM","YNKM","YNKM","YNKM","YNKM","FJFZ","FJFZ","FJFZ","FJFZ","FJFZ","FJFZ","FJFZ","FJFZ","FJFZ","FJFZ","FJFZ","FJFZ","FJFZ","FJFZ","FJFZ","FJFZ","FJFZ","FJFZ","FJFZ","FJFZ","SCPZ","SCPZ","SCPZ","SCPZ","SCPZ","SCPZ","SCPZ","SCPZ","SCPZ","SCPZ","SCPZ","SCPZ","SCPZ","SCPZ","SCPZ","SCPZ","SCPZ","SCPZ","SCPZ","SCPZ","HNCS","HNCS","HNCS","HNCS","HNCS","HNCS","HNCS","HNCS","HNCS","HNCS","HNCS","HNCS","HNCS","HNCS","HNCS","HNCS","HNCS","HNCS","HNCS","HNCS","JXNA","JXNA","JXNA","JXNA","JXNA","JXNA","JXNA","JXNA","JXNA","JXNA","JXNA","JXNA","JXNA","JXNA","JXNA","JXNA","JXNA","JXNA","JXNA","JXNA","JXNB","JXNB","JXNB","JXNB","JXNB","JXNB","JXNB","JXNB","JXNB","JXNB","JXNB","JXNB","JXNB","JXNB","JXNB","JXNB","JXNB","JXNB","JXNB","JXNB","JXNC","JXNC","JXNC","JXNC","JXNC","JXNC","JXNC","JXNC","JXNC","JXNC","JXNC","JXNC","JXNC","JXNC","JXNC","JXNC","JXNC","JXNC","JXNC","ZJNB","ZJNB","ZJNB","ZJNB","ZJNB","ZJNB","ZJNB","ZJNB","ZJNB","ZJNB","ZJNB","ZJNB","ZJNB","ZJNB","SCCD","SCCD","SCCD","SCCD","SCCD","SCCD","SCCD","SCCD","SCCD","SCCD","SCCD","SCCD","SCCD","SCCD","SCCD","SCCD","SCCD","SCGY","SCGY","SCGY","SCGY","SCGY","SCGY","SCGY","SCGY","SCGY","SCGY","SCGY","SCGY","SCGY","SCGY","SCGY","SCGY","SCGY","SCGY","SCGY","SCGY","NMAL","NMAL","NMAL","NMAL","NMAL","NMAL","NMAL","NMAL","NMAL","NMAL","NMAL","NMAL","NMAL","NMAL","NMAL","NMAL","NMAL","NMAL","BJYA","BJYA","BJYA","BJYA","BJYA","BJYA","BJYA","BJYA","BJYA","BJYA","BJYA","BJYA","BJYA","BJYA","BJYA","BJYA","BJYA","BJYA","BJYB","BJYB","BJYB","BJYB","BJYB","BJYB","BJYB","BJYB","BJYB","BJYB","BJYB","BJYB","BJYB","BJYB","BJYB","BJYB","BJYB","BJYB","BJYB","BJYC","BJYC","BJYC","BJYC","BJYC","BJYC","BJYC","BJYC","BJYC","BJYC","BJYC","BJYC","BJYC","BJYC","BJYC","BJYC","BJYC","BJYC","BJYC","BJYC","HBZB","HBZB","HBZB","HBZB","HBZB","HBZB","HBZB","HBZB","HBZB","HBZB","HBZB","HBZB","HBZB","HBZB","HBZB","HBZB","HBZB","HBZB","HBZB","HBZA","HBZA","HBZA","HBZA","HBZA","HBZA","HBZA","HBZA","HBZA","HBZA","HBZA","HBZA","HBZA","HBZA","HBZA","HBZA","HBZA","HBZA","HBZA","HBZA"))
#####remove NA#####
#toRemove <- is.na(glMean(genlight, alleleAsUnit = T))
#genlight <- genlight[, !toRemove]

library(radiator)
library(genepopedit)

radiator::detect_genomic_format("px_ddRAD_2018_432_5074_all_4X_0.999_thin_south1.gen")
pathtopgdspider="/usr/local/bin/"

pgdSpideR(input = paste0(output_dir,"px_ddRAD_2018_432_5074_all_4X_0.999_thin_south1.gen"),
          input_format="GENEPOP",
          output = paste0(output_dir,"px_ddRAD_2018_432_5074_all_4X_0.999_thin_south1.vcf",
                          output_format="vcf",
                          spid="c:/Users/YOURNAME/Documents/spids/GENEPOP_vcf.spid",
                          where.pgdspider=pathtopgdspider))

str <- genepop_structure(genepop="px_ddRAD_2018_432_5074_all_4X_0.999_thin_south1.gen",
                  path = paste0(output_dir,
                  "px_ddRAD_2018_432_5074_all_4X_0.999_thin_south1.str"))

genlight <- genomic_converter(str = "px_ddRAD_2018_432_5074_all_4X_0.999_thin_south1.str", 
                              strata = "2018_south1.strata.txt",
                              output = c("genlight"))

##automatically get population names from genepop file
filepath="px_ddRAD_2018_432_5074_all_4X_0.999_thin_genepop.txt"
output_dir <- "~/Dropbox/dbm_migration_ddRAD/data/genpop/" 
genepopfile="px_ddRAD_2018_432_5074_all_4X_0.999_thin_SampleID_fixed.genepop"
genepop_ID(genepop=filepath, path=paste0(output_dir, genepopfile))


##get variable from genepop using genepop_detective
PopNames <- genepop_detective(genepopfile, variable="Pops")
PopCounts <- genepop_detective(genepopfile, variable="PopNum")
SampleIDs <- genepop_detective("px_ddRAD_2018_432_5074_all_4X_0.999_thin_south1.genepop", variable="Inds")
LociNames <- genepop_detective(genepopfile, variable="Loci")
metadata <- genepop_detective(genepopfile,variable="All")
Alleles <- genepop_detective(genepopfile,variable="Allele")

##genepop_allelefreq 
AlleleFreq_default <- genepop_allelefreq(genepopfile)

pathtoplink=""
pathtopgdspider="/usr/local/bin/"
##genepop_filter_maf 
genepop_filter_maf(genepop=genepopfile, 
                   where.plink=pathtoplink,
                   where.pgdspider=pathtopgdspider,
                   maf=0.07,
                   path=paste0(output_dir,
                  "Genepop_maf_0-07.txt")) 

##genepop_toploci 
TopLoci <- genepop_toploci(genepop=genepopfile,
                           where.plink=pathtoplink,
                           where.pgdspider=pathtopgdspider,
                           fst.threshold = 0.04)

##subset_genepop
#inspect the names of the loci within the Genepop data
genepop_detective(genepopfile,"Loci")

#subset the Genepop file and 'keep' the specified loci names.
subloci <-c("Loci03","Loci15","Loci23","Loci49","Loci62","Loci81","Loci88","Loci94")
subset_genepop(genepop= genepopfile, keep = TRUE, subs = subloci, path = paste0(output_dir,"Genepop_Loci_selection.txt"))
#subset the Genepop file and do not 'keep' the specified 'subs'
subset_genepop(genepop= genepopfile, keep = FALSE, subs = subloci, path = paste0(output_dir,"Genepop_Loci_neutral.txt"))

##remove population
#use genepop_detective to find the population names
PopNames <- genepop_detective(genepopfile,"Pops")
PopNames
#create a list of populations you want to keep manual
PopKeep <- c("HNHK", "GDSZ", "GXNN", "GDGA", "GDGB", "GDGC", "YNYX", "YNDH", "YNKM", "FJFZ", "SCPZ", "HNCS","JXNA", "JXNB", "JXNC", "ZJNB", "SCCD")
#or if you have fewer pops to remove
#PopKeep <- setdiff(PopNames, c("CCC","GGG"))

#vector of loci to keep
#subloci <- c("Loci03","Loci15","Loci23","Loci49","Loci62","Loci81","Loci88","Loci94")
subset_genepop(genepop= genepopfile, keep = TRUE, 
               spop = PopKeep, 
               path = paste0(output_dir,"px_ddRAD_2018_432_5074_all_4X_0.999_thin_south1.genepop"))  ##subs = subloci,

##subset_genepop_rename
#view the population names
genepop_detective(GenePopData,"Pops")

#create a dataframe for renaming. 
#column 1 = the original pop names.
#column 2 = the new names. 
PopRename <- data.frame(oPop = c("AAA","BBB","CCC","DDD","EEE","FFF","GGG","HHH","III","JJJ"),
                        newname = c("AAA","BBB","CCC","YYY","EEE","FFF","GGG","ZZZ","III","JJJ"))

#rename populations
subset_genepop_rename(genepop= GenePopData, path = paste0(output_dir,"Genepop_renamed.txt"),nameframe = PopRename)

##Populations can also be grouped by a common name. The result is similar to subset-genepop-aggregate except that population names will be replaced. Here it is useful to enable renumbering (renumber = TRUE) or there individuals with the same number will be assigned the same ID, if populations are combined.
#view the population names
genepop_detective(GenePopData,"Pops")

#create a dataframe for renaming. 
#column 1 = the original pop names.
#column 2 = the new names. 
PopRename_group <- data.frame(oPop = c("AAA","BBB","CCC","DDD","EEE","FFF","GGG","HHH","III","JJJ"),newname = c("Pop1","Pop1","Pop1","Pop1","Pop1","Pop2","Pop2","Pop2","Pop2","Pop2"))

#rename populations
subset_genepop_rename(genepop= GenePopData, nameframe = PopRename_group, renumber = TRUE, meta="Pop",path = paste0(output_dir,"Genepop_renamed_renumbered.txt"))

##Alternatively the sampleIDs can be renamed directly.
#create a renaming vector (column 1 = old column 2 = new)
IndRename_group <- data.frame(old = c("AAA_10","BBB_10","CCC_10","DDD_10"),new = c("AAA_99","BBB_99","CCC_99","DDD_99"))

#rename sampleIDs
subset_genepop_rename(genepop= GenePopData, nameframe = IndRename_group, meta="Ind",path = paste0(output_dir,"Genepop_ID_renamed.txt"))

subset_genepop_aggregate

Now lets change group some populations together. In this example we will combine populations DDD/EEE & FFF/HHH and remove populations JJJ/GGG (not listed in column 1 of PopRename). No loci will be removed from this dataframe (subs = NULL). This is useful when you don't want a clustering program (e.g. STRUCTURE or Bayescan) to assume differences among population groups. Function description

#Use genepop_detective to find the meta-data for populations
genepop_detective(GenePopData,"Pops")

#create a dataframe for renaming. 
#column 1 = the original pop names and the list of pops required.
#column 2 = Grouping variables. 
PopAggregate <-data.frame(oP = c("AAA","BBB","CCC","DDD","EEE","FFF","HHH","III"),
agname = c("AAA","Pop1","CCC","Pop1","DDD","Pop2","HHH","Pop2"))

#re-cast Genepop format to group populations based on 'PopAggregate'
subset_genepop_aggregate(genepop= GenePopData, subs = NULL, path = paste0(output_dir,"Genepop_grouped.txt"),agpopframe = PopAggregate)

Now that we have two grouped populations, we can use the population rename function to give them a common name.

#read in the grouped Genepop file
GenePopData2 <-read.table(paste0(output_dir,"Genepop_grouped.txt"),header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)

#investigate the population parameters
genepop_detective(GenePopData2,"Pops")

#we can see that "BBB" & "DDD" and "FFF" & "III" are at the end but still have different names. 

#create renaming frame
PopRename <- data.frame(oPop = c("AAA","CCC","EEE","HHH","BBB","DDD","FFF","III"),
newname = c("AAA","CCC","EEE","HHH","Group1","Group1","Group2","Group2"))

#Rename the grouped populations
subset_genepop_rename(genepop= GenePopData2, nameframe = PopRename, path = paste0(output_dir,"Genepop_grouped_renamed.txt"))


subset_genepop_individual

Remove some individual samples from the dataframe. Here we consider the 'individual' ID to be an alpha-numeric code which at the start of each row of data. The keep parameter defines whether the specified list is removed or retained (default: keep = FALSE). This function is best run at the beginning of the analysis prior to removal of specific loci or populations. Function description

#Use genepop_detective to see the naming structure.
genepop_detective(GenePopData,"Inds")

#vector of sample IDs to remove
subid <- c("AAA_1","AAA_3","BBB_20","CCC_21","EEE_3","EEE_26","HHH_25","JJJ_4")

subset_genepop_individual(genepop= GenePopData, indiv = subid, keep = FALSE, path = paste0(output_dir,"Genepop_IDsubset.txt"))

genepop_reorder

Reorder the genepop file sequentially according to a specified ordering (see below for examples). Function description

#Identify populaton names
pops <- genepop_detective(genepop="Genepop_IDsubset.txt")

#reorder populatons. In this example the first and last example will be switched 
genepop_reorder(genepop="Genepop_IDsubset.txt",reorder=pops[c(length(pops),2:(length(pops)-1),1)]),path=paste0(output_dir,"GenePop_reordered.txt")

Sampling
genepop_sample

Function which will allow you to create subsets of the genepop. This might be useful for creating training and leave-out datasets for assignment. Function description

#view how many individuals you have in each population
genepop_detective(GenePopData,"PopNum")

There are several ways we can sub-sample each population

Use a fixed number (e.g. 5 from each population)

SubSamp <- genepop_sample(GenePopData, nsample = 5)

or

Use a fixed fraction (e.g. 25% from each population)

SubSamp <- genepop_sample(GenePopData, nsample = 0.25)

or

Use a fraction which varies for each population

#vector of populations in the Genepop file
Pops <- genepop_detective(GenePopData,"Pops")

#create a dataframe which defines how many you want to sample from each population. Note here we create a column which has an entry for each population.
subdf <- data.frame(Pops = Pops, nsample = c(5, 2, 3, 5, 6, 5, 2, 3, 5, 6))

SubSamp <- genepop_sample(GenePopData, nsample = subdf)

or

Use a fraction which varies for each population

#vector of populations in the Genepop file
Pops <- genepop_detective(GenePopData,"Pops")

#create a dataframe which defines how many you want to sample from each population. Note here we create a column which has an entry for each population.
subdf <- data.frame(Pops = Pops, nsample = c(5, 2, 3, 5, 6, 5, 2, 3, 5, 6)/10)

SubSamp <- genepop_sample(GenePopData, nsample = subdf)

Once you select a method to sub-sample each population, you can use the function subset-genepop-indiv to create the sampled genepop file.

For example if you wanted to create a training and assignment dataset in genepop format using the random stratified sampled we selected using subset-genepop-indiv:

#Create training dataset (keep = TRUE)
subset_genepop_individual(GenePopData, indiv = SubSamp, keep = TRUE, path = paste0(output_dir,"Genepop_training.txt"))

#Create an assignment dataset using the remaining individuals (keep = FALSE)
subset_genepop_individual(GenePopData, indiv = SubSamp, keep = FALSE, path = paste0(output_dir,"Genepop_assignment.txt"))

Conversion
genepop_structure

If you are interested in investigating population structure you can convert your modified Genepop object or path to saved file directly into a STRUCTURE formatted text (.str) file <>. Function description

#convert Genepop format to STRUCTURE to default groupings
genepop_structure(genepop="Genepop_IDsubset.txt",path = paste0(output_dir,"Sturcture_IDsubset_groups.txt"))

# Specify population groupings (group CCC/DDD/EEE and III/JJJ)
pGroups <- data.frame(pops = c("AAA","BBB","CCC","DDD","EEE","FFF","GGG","HHH","III","JJJ"),groups = c("1","2","3","3","3","4","5","6","7","7"))

#convert Genepop format to STRUCTURE with new groupings
genepop_structure(genepop="Genepop_IDsubset.txt",popgroup = pGroups, path = paste0(output_dir,"Sturcture_IDsubset_groups.txt"))

genepop_fstat

If you are interested in doing calculating gene diversities and-or F-statistics with the cleaned loci datasets, you can convert from Genepop to FSTAT format This data can be used by the R package hierfstat or the program FSTAT Function description

#convert Genepop format to FSTAT (.dat)
genepop_fstat(genepop="Genepop_IDsubset.txt",path = paste0(output_dir,"Fstat_IDsubset.dat"))

#convert Genepop format to FSTAT but keep the data in the workspace
genepop_fstat(genepop="Genepop_IDsubset.txt",addworkspace = TRUE)

genepop_newhybrids

If you are interested in testing for hybridization you can convert from Genepop to the format necessary for the program New Hybrids Function description

#convert Genepop format to New Hybrids format (.txt)
genepop_newhybrids(genepop="Genepop_IDsubset.txt",path = paste0(output_dir,"NewHybrids_IDsubset.txt"))

genepop_assigner

If you are interested conducting assignment analysis, you can convert to the format required for the R package assigner Function description

#convert Genepop format to assigner format (.txt) using the populations each as their own assessment level
genepop_assigner(genepop="Genepop_IDsubset.txt",path = paste0(output_dir,"assigner_IDsubset.txt"))

# Set new population assignment groupings (group CCC/DDD/EEE and III/JJJ)
pGroups <- data.frame(pops = c("AAA","BBB","CCC","DDD","EEE","FFF","GGG","HHH","III","JJJ"),groups = c("1","2","3","3","3","4","5","6","7","7"))

genepop_assigner(genepop="Genepop_IDsubset.txt",popgroup = pGroups, path = paste0(output_dir,"assigner_IDsubset_NewGroups.txt"))


genepop_colony

If you are interested in investigating sibship or parentage, you can convert to the format required for the software Colony. Returned are functions for translating the Colony output and those required by the maximum likelihood analysis. 'Individual' and 'Loci' conversion files can be used to match the loci and individual names formatted for input to Colony to those provided as input by the user.'MarkerTypeErrorRT' file defines the allele dropout rates as calculated using the 'missing' function in PLINK. Rates of other error types are assumed to be equal to the estimated allele dropout rates. 'GENOTYPES' are the geneotypes in two columns per loci format. Output from Colony is based on the input for of the 'GENOTYPES' file. Original loci and individual names can be assigned using the conversion output files. All files are returned to path and are named according tot he input filename. Function description

#convert Genepop format to the files necessary for Colony 
genepop_colony(genepop="Genepop_IDsubset.txt",where.plink="c:/Users/YOURNAME/Documents/Programs/plink/",where.pgdspider="c:/Users/YOURNAME/Documents/Programs/PGDSpider_2.0.9.0/",denote.missing = "000",path = output_dir)


genepop_bgc

genepopedit is one of the few tools available to convert data to the format required for input into to the Bayesian estimation of Genomic Clines (BGC) format. BGC can evalute genomic clinal patterns and introgression among loci. To convert to BGC format you must specify which populations you consider to be ancestral or parental (P1, P2) and which you assume could be hybridized (admixed). Function description

#specify which populations are going to be included in the analysis and to which class they belong. Note Pops identified in P1 and-or P2 can also be specified as "Admixed" to test BGC output. 
BGC_groups = data.frame(pops = c("AAA","BBB","CCC","DDD","EEE","FFF","GGG","HHH","III","JJJ"),groups = c("P1","P1","P1","Admixed","Admixed","Admixed","Admixed","P2","P2","P2"))

#convert Genepop to BGC input files (3). Note in this case the variable path is a path to the directory where the input files will be saved.
genepop_bgc(genepop="Genepop_IDsubset.txt",popdef = BGC_groups, fname="BGC_IDsubset",path = output_dir)

genepop_treemix

genepopedit is one of the few tools available to convert data to the format required for input into the Treemix format. It creates a gzipped file that clusters your individuals based on the population ID of the individual's code. Population IDs are extracted using the "_" to differentiate population from sample ID (e.g. BON_01 is population BON sample 01). If samples are not separated refer to genepop_ID. Grouping levels for populations can also be changed using subset_genepop_rename. Note that this is just the first of two steps in the process; the gzipped output from this function needs to be run through the Python script (Python 2.7+ must be installed first) that accompanies Treemix. This output from "plink2treemix.py" is then ready to run in Treemix. This Python script can be downloaded from link. Treemix is used to infer migration weight and directionality among your populations, as well as detecting population splits and mixtures. Function description

#convert Genepop to Treemix input and keep intermediary conversion files (default: FALSE)
genepop_treemix(genepop="Genepop_IDsubset.txt",where.plink="c:/Users/YOURNAME/Documents/Programs/plink/",where.pgdspider="c:/Users/YOURNAME/Documents/Programs/PGDSpider_2.0.9.0/",keep_inter = TRUE, path = paste0(output_dir,"Treemix_IDsubset.txt"))

genepop_GSIsim

If you are interested in assessing the accuracy of a genetic stock identification analysis, given a genetic baseline, you can convert directly from GENEPOP to a GSI_sim formatted text file. genepopedit is among the few tools currenlty available to format data for GSI_sim. Function description

#convert Genepop to GSIsim format
genepop_GSIsim(genepop="Genepop_IDsubset.txt",path = paste0(output_dir,"GSIsim_IDsubset.txt"))

genepop_hzar

If you are interested in exploring clines in allele frequency over distance, you can convert from Genepop to HZAR format.This data can be used by the R package hzar Function description

#convert Genepop format to HZAR (.csz)
genepop_hzar(genepop="Genepop_IDsubset.txt",distances=data.frame(Pop=genepop_detective("Genepop_IDsubset.txt"), Distance=c(5,8,22,3,84,39,45.5,67,101),stringsAsFactors=FALSE), path = paste0(output_dir,"HZAR_Input.csv"))

genepop_flatten

#Flatten your geneotype data for other calculations. This function will convert and return the genepop format into a dataframe. Because this function returns the dataframe, it needs to be assigned a variable ID. Note that this function does not have any subsetting capabilities. If you would like to remove loci and-or populations, this must be before hand using subset_genepop or its sister functions (see below for examples). Function description

#Flatten the dataframe
GenePop_df<- genepop_flatten(GenePopData)

#inspect the output for the first 10 columns 
head(GenePop_df[,1:10])

genepop_unflatten

If you are working with a flattened dataframe in your workspace you can convert it back to Genepop format. **Note that the first column of this data.frame should correspond to the sampleID (e.g. "BON_01") and the remaining columns should be Loci. Function description

# create a Genepop format which can be used by other genepopedit functions
genepop_unflatten(GenePop_df, path = paste0(output_dir,"GenePop_UNFLATTENED.txt"))

Conversion using PGDspider
PGDspideR

There are a broad range of data formats which can be converted to and from using PGDspider. If you have multiple files which you are going to be converting using the same conversion parameters, PGDspideR.R can be useful, providing a code based interface to exploit the conversion functions of PGDspider. Function description

#convert between GENEPOP and FSTAT format using PGDspider
pgdSpideR(input = paste0(output_dir,"Genepop_IDsubset.txt"),
          input_format="GENEPOP",
          output = paste0(output_dir,"Genepop_IDsubset_FSTAT.dat",
                          output_format="FSTAT",
                          spid="c:/Users/YOURNAME/Documents/spids/GENEPOP_FSTAT.spid",
                          where.pgdspider="c:/Users/YOURNAME/Documents/Programs/PGDSpider_2.0.9.0/")
          
          #Simulate genotypes using pooled DNA allele frequencies
          alleleotype_genepop
          
          Individual geneotypes can be simulated using the alleleotype_genepop() function and formatted to GENEPOP for use in conventional genomic analyses (e.g. assigner, DAPC, STRUCTURE). This function will simulate individuals for each population provided by the input file. The output can then be diagnosed, sampled, and formatted using other genepopedit functions. Function description
          
          #simulate individual geneotypes allele frequencies stratified by population.
          #Create an example pooled DNA input file
          pooledDNA <- data.frame(Pop = c("Pop1","Pop2","Pop3","Pop4","Pop5"),Loci1 = round(runif(5, 0, 1),2),Loci2 = round(runif(5, 0, 1),2),Loci3 = round(runif(5, 0, 1),2),Loci4 = round(runif(5, 0, 1),2),Loci5 = round(runif(5, 0, 1),2),Loci6 = round(runif(5, 0, 1),2))
          
          #inspect the input
          pooledDNA
          
          #simulate 200 individual geneotypes for each population
          alleleotype_genepop(pooledDNA, numsim = 200, path = paste0(output_dir,"SimulatedGeneotypes.txt"))
          
          #To validate output simulate 100 (default) geneotypes per population and compare allele frequencies to those of the input file. Estimated frequencies per population and SNP should match those of the input file within +/- 1%.
          
          alleleotype_genepop(pooledDNA, numsim = 100, path = paste0(output_dir,"SimulatedGeneotypesValidate.txt"))
          
          genepop_allelefreq(paste0(output_dir,"SimulatedGeneotypesValidate.txt"))
          
          #Note you can simulate new geneotypes using the summary from genepop_allelefreq(... , wide = TRUE).
          GeneFreq <- genepop_allelefreq("Genepop_IDsubset.txt",wide = TRUE)
          
          alleleotype_genepop(GeneFreq, numsim = 100, path = paste0(output_dir,"SimulatedGeneotypes_Genepop_IDsubset.txt"))




