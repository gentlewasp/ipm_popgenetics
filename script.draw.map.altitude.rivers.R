par(mar=rep(0,4))
dat = read.csv(text = "city,jd,wd
SL,103.19,24.51
ND,119.32,26.39
NB,121.32,29.52
WH,114.18,30.35
NJ,118.47,32.03
CD,104.03,30.39
TS,105.43,34.34
YL,108.05,34.16
ZZ,113.37,34.44
YA,109.29,36.35
TG,112.33,37.25
QD,120.22,36.04
SJ,114.30,38.02
BD,115.27,38.52
DX,116.20,39.43
MT,116.06,39.56
PG,117.07,40.08
YQ,115.58,40.27
DL,121.36,38.54
SY,123.25,41.48
KD,128.36,35.52")
library(maps)
library(mapdata)
library("dismo")
# china, world, worldHires, world2,
filename = "ChinaKorea_altitude"
pdffile=paste(filename,".pdf", sep = "")
pdf(file=pdffile, width=10, height=10)
altM <- getData("alt", country = "China", mask= FALSE)  # this will download global data on altitude at 10' resolution
plot(altM)
map("china", 
    col = "darkgray", 
    xlim = c(90,130),
    ylim = c(20,54), 
    panel.first = grid(),
    add = TRUE)    
map("worldHires", 
    col = "darkgray", 
    xlim = c(90,130),
    ylim = c(20,54), 
    panel.first = grid(),
    add = TRUE)
map('rivers', add=TRUE, col="red")
points(dat$jd, dat$wd, 
       pch = 19, 
       col = rgb(0.9, 0, 0, 0.8))
text(dat$jd, dat$wd, 
     dat[, 1], 
     cex = 0.9, 
     col = rgb(0, 0, 0, 1.0), 
     pos = c(2, 4, 4, 4, 3, 4, 2, 3, 4, 2, 4, 2, 2,
    4, 3, 2, 1, 3, 1, 1, 2, 3, 2, 2, 1, 2, 4, 3, 1, 2, 2, 4, 4, 2))
axis(1, lwd = 0); axis(2, lwd = 0); axis(3, lwd = 0); axis(4, lwd = 0)
dev.off()

"
filename = "worldarea"
pdffile=paste(filename,".pdf", sep = "")
pdf(file=pdffile, width=10, height=10)
map("worldHires", 
    col = "darkgray", 
    xlim = c(70,136),
    ylim = c(15,55), 
    panel.first = grid())
dev.off()
"