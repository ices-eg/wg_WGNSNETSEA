library(pacman)
p_load(data.table, ggplot2, mapplots, shapefiles, plyr, openxlsx)

setwd("Q:/20-forskning/20-dfad/users/jostou/home/26_10_27_seaNet")
coast <- read.shapefile("data/europe")

dat <- rbind.fill(read.xlsx("data/LEMON SOLE DATA MIK Q3_FOR MAPS.xlsx", sheet = 1),
                  read.xlsx("data/LEMON SOLE DATA MIK Q3_FOR MAPS.xlsx", sheet = 2),
                  read.xlsx("data/LEMON SOLE DATA MIK Q3_FOR MAPS.xlsx", sheet = 3),
                  read.xlsx("data/LEMON SOLE DATA MIK Q3_FOR MAPS.xlsx", sheet = 4))

#limit of the plot
xlim <- c(-4, 12)
ylim <- c(49, 62)

for (yr in unique(dat$YEAR)) {
  
  xx <- dat[dat$YEAR == yr, ]
  
  dat2 <- xx[xx$LemonSoleN > 0, ]
  zero <- xx[xx$LemonSoleN == 0, ]
  
  nam <- "nr caught"

  # add a point to be used as a sort of label....
  dat2 <- rbind.fill(dat2, data.frame(posLonStartDec = rep(9, 4), 
                                      posLatStartDec = c(49.70, 50, 50.4, 51), 
                                      LemonSoleN = quantile(dat2$LemonSoleN)[-1]))
  
  #plot the plot
  basemap(xlim, ylim, main = paste0(nam, " - ", yr, "- Q3"),
          xaxt = "n", yaxt = "n", xlab="", ylab="")
  draw.rect()
  draw.shape(coast, col="cornsilk")
  rect(7.5, 49.4, 11.5, 52, density = NULL, angle = 45,
       col = "lightblue", border = NULL, lty = par("lty"), lwd = par("lwd"))
  draw.bubble(dat2$posLonStartDec, dat2$posLatStartDec, 
              dat2$LemonSoleN, maxradius = 0.8, pch=21, bg="#ABABAB95")
  axis(1, at = seq(-4, 12, 1), 
       labels = c(paste0(seq(4,1, by = -1), "°W"),
                  "0°", paste0(seq(1, 12, 1), "°E")))
  axis(2, at = seq(49, 62, 1), 
       labels = paste0(seq(49, 62, 1), "°N"))
  text(x = c(7.7, rep(10, 4)), y = c(51.8, 49.8, 50.15, 50.52, 51.1), 
       labels = c(nam, round(quantile(dat2$LemonSoleN)[-1], 2)),
       adj = c(0,1))
  draw.bubble(zero$posLonStartDec, zero$posLatStartDec,
              1, maxradius = 0.05, pch=4, bg="red")
  
  #export the plot
  dev.copy(png, paste0("lemon/mikQ3/lemon_mik_", yr, "_Q3.png"), 
           width = 180, height = 240, units = 'mm', res = 300)
  dev.off()
  
  
}
