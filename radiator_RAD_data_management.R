##home page:https://github.com/thierrygosselin/radiator
citation("radiator")
##installation
install.packages("future")
if (!require("devtools")) install.packages("devtools") # to install
devtools::install_github("thierrygosselin/radiator")
##if failed install using above method, download radiator_pkg_install.R from following link, and run ...
#download.file(url="https://www.dropbox.com/s/7ekjvqx2qahg8mg/radiator_pkg_install.R", destfile="radiator_pkg_install.R")
#source("radiator_pkg_install.R")
#rad <- radiator_pkg_install()

library(radiator)
##It's a tab separated file, e.g. radiator.strata.tsv.
##A minimum of 2 columns: INDIVIDUALS and STRATA is required.
##The STRATA column identifies the individuals stratification, the hierarchical groupings: populations, sampling sites or any grouping you want.
##It's like stacks population map file with header...
strata <- radiator::read_strata("px_ddRAD_2018_432_5074_all_4X_0.999_thin_.strata.tsv")
names(strata)

##input and output files
##test input data format
radiator::detect_genomic_format("px_ddRAD_2018_432_5074_all_4X_0.999_thin.vcf")
##filter data
data <- radiator::filter_rad(data = "px_ddRAD_2018_432_5074_all_4X_0.999_thin.vcf", 
                             strata = "px_ddRAD_2018_432_5074_all_4X_0.999_thin_.strata.tsv", 
                             output = c("genepop"))
##other formats: "genind", "hierfstat", "genepop","plink", "structure", "faststructure", "arlequin", "hierfstat", "gtypes", "bayescan", "betadiv", "pcadapt", "hzar", "fineradstructure", "related", "seqarray", "snprelate","maverick"

radiator::data.structure <- write_structure(data)

radiator::genomic_converter(data="px_ddRAD_2018_432_5074_all_4X_0.999_thin.vcf",
                            strata = "px_ddRAD_2018_432_5074_all_4X_0.999_thin_.strata.tsv",
                            filename="px_ddRAD_2018_432_5074_all_4X_0.999_thin.genpop",
                            output = c("structure"))
