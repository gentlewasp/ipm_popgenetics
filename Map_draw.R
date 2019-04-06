################################################################################
######### This is a script for drawing the sample site distribution map ########
################################################################################

### install packages below if required
library(ggplot2)
library(scatterpie)
library(readxl)
setwd("~/Desktop")
### read the Excel file as input, which should in your working dictionary
### type getwd() for looking up the current working dictionary
my_data <- read_excel("test_data_map.xlsx")

### extract the clusters names as a vector for later use
a <- colnames(my_data)[1:4]
b <- colnames(my_data)
clusters <- b[! b %in% a]

### plot the sample site region map as the first layer
world <- map_data('world2')
p2 <- ggplot(world, aes(long, lat)) +
  geom_map(map=world, aes(map_id=region), fill="darkgrey", color="#7f7f7f") +
  coord_equal(xlim=c(min(my_data$long)-2,max(my_data$long)+2), ylim=c(min(my_data$lat)-2,max(my_data$lat)+2)) + 
  theme_bw() + geom_point(aes(x=long, y=lat), data=my_data, color="red", size=0.6) +
  geom_text(aes(x=long-1, y=lat, label=pop), data=my_data, size = 2, hjust = "right") + 
  theme(panel.grid.major=element_line(linetype ="dotted"), panel.grid.minor=element_line(linetype="dotted"), 
        axis.title.x=element_blank(), axis.title.y=element_blank())

### plot the pie chart representing populations and clusters structure
p2 + geom_scatterpie(aes(x=long+1, y=lat+1.5, group=pop, r=number*0.05), data=my_data, cols=clusters, color=NA, alpha=0.8) + 
  theme(legend.title = element_blank(), legend.key.size = unit(0.2, "cm"), 
        legend.key.width = unit(0.3,"cm"), legend.text = element_text(size=7), legend.position = c(0.9,0.1))

