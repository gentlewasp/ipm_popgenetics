
#http://membres-timc.imag.fr/Olivier.Francois/LEA/tutorial.htm
#http://membres-timc.imag.fr/Olivier.Francois/snmf/tutorial.htm
install.packages(c("fields","RColorBrewer","mapplots"))
source("https://bioconductor.org/biocLite.R")
biocLite("LEA")
biocLite("vignette")
biocLite("tess3r")
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")

##import data
### import vcf file###
radiator::genomic_converter(data="px_ddRAD_2018_432_5074_all_4X_0.999_thin_south1.gen",
                            strata = "px_ddRAD_2018_432_5074_all_4X_0.999_thin_.strata.tsv",
                            filename="px_ddRAD_2018_432_5074_all_4X_0.999_thin",
                            output = c("structure"))

LEA::vcf2geno("px_ddRAD_2018_432_5074_all_4X_0.999_thin.vcf", "px_ddRAD_2018_432_5074_all_4X_0.999_thin.geo")

library(LEA)
kValue=2
obj.snmf = snmf("px_ddRAD_2018_432_5074_all_4X_0.999_thin_south1.gen", K = kValue, alpha = 100, project = "new")
qmatrix = Q(obj.snmf, K = kValue)

barplot(t(qmatrix), col = c("black","orange","violet","lightgreen","red"), border = NA, space = 0,
        xlab = "Individuals", ylab = "Admixture coefficients")

coord = read.table("coordinates.coord")
pop = rep(1:60, each = 10)

K = 3
Npop = length(unique(pop))
qpop = matrix(NA, ncol = K, nrow = Npop)
coord.pop = matrix(NA, ncol = 2, nrow = Npop)
for (i in unique(pop)){
  qpop[i,] = apply(qmatrix[pop == i,], 2, mean)
  coord.pop[i,] = apply(coord[pop == i,], 2, mean)}

library(mapplots)
plot(coord, xlab = "Longitude", ylab = "Latitude", type = "n")
map(add = T, col = "grey90", fill = TRUE)
for (i in 1:Npop){
  add.pie(z = qpop[i,], x = coord.pop[i,1], y = coord.pop[i,2], labels = "",
          col = c("orange","violet","lightgreen"))}

obj.snmf = snmf("secondary_contact.geno", K = 1:8, ploidy = 2, entropy = T,
                alpha = 100, project = "new")
plot(obj.snmf, col = "blue4", cex = 1.4, pch = 19)

