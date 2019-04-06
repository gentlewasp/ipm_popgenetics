library("poppr")
library("vcfR")
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(cowplot)
library(adegenet)
library("snowfall")
install.packages("snowfall")
### import vcf file###
vcf <- vcfR::read.vcfR("populations.snps.vcf", verbose = FALSE) 
####converge to genlight####
genlight <- vcfR2genlight(vcf)
####define population of each sample####
pop(genlight) <- as.factor(c("BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX",	"BJYQX","BJYQX","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP","BJYQP"))
#####################
#####simple PCA####
######################
x <- tab(genlight,NA.method="mean")
pca2 <- prcomp(x, scale=T)
#x11(height=6, width=12, pointsize=12); par(mfrow=c(1,2)) 
mycolors <- rainbow(19) 
plot(pca2$x, pch=20, col=mycolors)
plot(pca2$x, type="n"); text(pca2$x, rownames(pca2$x), cex=0.8, col=mycolors) 
#install.packages("geneplotter")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("geneplotter", version = "3.8")

library(geneplotter); smoothScatter(pca2$x) # Same as above, but generates a smooth scatter plot that shows the density of the data points.
pairs(pca$x[,1:4], pch=20, col=mycolors[sort(rep(1:5, 500))]) 
# Plots scatter plots for all combinations between the first four principal components.
biplot(pca2) 
# Plots a scatter plot for the first two principal components plus the corresponding eigen vectors that are stored in pca$rotation.
library(scatterplot3d) # Loads library scatterplot3d.
scatterplot3d(pca$x[,1:3], pch=20, color=mycolors[sort(rep(1:5, 500))]) 
# Same as above, but plots the first three principal components in 3D scatter plot.
library(rgl); rgl.open(); offset <- 50; par3d(windowRect=c(offset, offset, 640+offset, 640+offset)); rm(offset); rgl.clear(); rgl.viewpoint(theta=45, phi=30, fov=60, zoom=1); spheres3d(pca$x[,1], pca$x[,2], pca$x[,3], radius=0.3, color=mycolors, alpha=1, shininess=20); aspect3d(1, 1, 1); axes3d(col='black'); title3d("", "", "PC1", "PC2", "PC3", col='black'); bg3d("white") 


#####remove NA#####
toRemove <- is.na(glMean(genlight, alleleAsUnit = T))
head(toRemove)
genlight <- genlight[, !toRemove]
View(genlight)



#####find optimal number of cluster####
cluster <- find.clusters(genlight, n.clust=NULL,                                                                        
                         max.n.clust=12,                                                                                      
                         stat=c("BIC"),                                                                                       
                         n.iter=1e9, n.start=1e3, # 1e9, 1e3                                                                  
                         #truenames=TRUE,
                         scale=FALSE,
                         parallel=T)

numberofcluster = 5

#####DAPC#####
dapc <- dapc(genlight, n.da=numberofcluster, n.pca=150)

pdf(file="dapc.assignscatter.pdf")                                                                                      
scatter(dapc,cstar=0,                                                                                                   
        #mstree=TRUE,                                                                                                          
        # posi.da="bottomright", posi.pca="bottomleft", scree.pca=TRUE,                                                         
        scree.da=FALSE,                                                                                                       
        pch=20,                                                                                                               
        leg=TRUE,                                                                                                             
        # col=seasun(14),                                                                                                       
        clab=0.8, # population names in the circle                                                                              
        ratio.pca=0.3, solid=.6, cex=3)                                                                                       
dev.off() 


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
            row.lab=c("BJYQ","HBYC","HLHE","HLMD", "LNXC", "NXWZ", "SCAB","SDLK","SDTA","SXJZ","SXXY","YNKM")) 
dev.off()
