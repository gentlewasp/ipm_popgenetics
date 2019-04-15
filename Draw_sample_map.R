######################################################################################
# This is a script for drawing the sample site distribution map, writtern by Yi Chen #
######################################################################################

### install packages below if required
library(ggplot2)
library(maps)
library(mapdata)
library(dismo)
library(scatterpie)
library(readxl)
library(grid)
library(rasterVis)

setwd("~/Desktop")
# read the Excel file as input, which should in your working dictionary
# type getwd() for looking up the current working dictionary
my_data <- read_excel("test_data_map.xlsx")

# extract the clusters names as a vector for later use
a <- colnames(my_data)[1:4]
b <- colnames(my_data)
clusters <- b[! b %in% a]

# load the world map data
world <- map_data('world2')  
# download global data on altitude at 10' resolution
altM <- getData("alt", country = "China", mask= FALSE, res=5)
alt <- altM[[1]]
# download China geographical data at level 1
china <- getData('GADM' , country="CHN", level=1)
china <- fortify(china)


# plot the elevation, map, path, point and text layers as p1, which will be the buttom plot
p1 <- gplot(alt) +
  geom_raster(aes(x = x, y = y, fill = value), show.legend = T, interpolate = T) +
  scale_fill_gradientn(colours = rev(terrain.colors(10)), name="Altitude") +
  geom_map(data=world, map=world, aes(x=long, y=lat, map_id = region), fill=NA, color="#7f7f7f") +
  geom_path(data=china, aes(long,lat, group=group), color="#7f7f7f", size=0.3) +
  geom_point(aes(x=long, y=lat), data=my_data, color="red", size=0.6) +
  geom_text(aes(x=long-1, y=lat, label=pop), data=my_data, size = 2, hjust = "right") +
  coord_equal(xlim = c(73.5, 134.9), ylim = c(18.1, 53.7), expand = F) +
  #coord_equal(xlim=c(min(my_data$long)-2,max(my_data$long)+2), ylim=c(min(my_data$lat)+0.5,max(my_data$lat)+2)) + 
  theme_classic() + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), 
        legend.key.size = unit(0.4, "cm"), legend.key.width = unit(0.3,"cm"), 
        legend.text = element_text(size=7), legend.title = element_text(size=7), legend.background = element_blank()) 


# create an empty theme to make elements on the top plot be transparent, which will be used in p2
new_theme_empty <- theme_bw()
new_theme_empty$line <- element_blank()
new_theme_empty$rect <- element_blank()
new_theme_empty$strip.text <- element_blank()
new_theme_empty$axis.text <- element_blank()
new_theme_empty$plot.title <- element_blank()
new_theme_empty$axis.title <- element_blank()
# t=0, r=0, b=0.5, l=0.5
new_theme_empty$plot.margin <- structure(c(0, 0, 0.5, 0.5), unit = "lines", valid.unit = 3L, class = "unit")


# plot the scatter pie layer as p2, which will be the top plot
p2 <- ggplot() + 
  geom_scatterpie(data=my_data, aes(x=long+1, y=lat+1.5, group=pop, r=number*0.04), cols=clusters, color=NA, alpha=0.8) +
  #geom_scatterpie_legend(my_data$number*0.04, x=80, y=20, n=3, labeller=function(x) x*25) +
  coord_equal(xlim = c(73.5, 134.9), ylim = c(18.1, 53.7), expand = F) + new_theme_empty + 
  theme(legend.title = element_blank(), legend.key.size = unit(0.2, "cm"), legend.key.width = unit(0.3,"cm"), 
        legend.text = element_text(size=7), legend.background = element_blank())


# extract the legend data of p1
tmp1 <- ggplot_gtable(ggplot_build(p1))
leg1 <- which(sapply(tmp1$grobs, function(x) x$name) == "guide-box")
legend1 <- tmp1$grobs[[leg1]]
# extract the legend data of p2
tmp2 <- ggplot_gtable(ggplot_build(p2))
leg2 <- which(sapply(tmp2$grobs, function(x) x$name) == "guide-box")
legend2 <- tmp2$grobs[[leg2]]


# overlay p2 on the top of p1
grid.newpage()
vp1 <- viewport(layout = grid.layout(1, 1, widths = unit(1 , "npc")))
pushViewport(vp1) 
# a bit slow, please wait!
print(p1 + theme(legend.position="none"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
# a bit slow, please wait!
print(p2 + theme(legend.position="none"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))

# overlay p2's legend
upViewport(0)
vp2 <- viewport(width = unit(1, "npc"), x = 0.9, y = 0.185)
pushViewport(vp2)
grid.draw(legend2)
# a bit slow, please wait!
popViewport()

# overlay p1's legend
upViewport(0)
vp3 <- viewport(width = unit(1, "npc"), x = 0.95, y = 0.2)
pushViewport(vp3)
grid.draw(legend1)
# a bit slow, please wait!
popViewport()

