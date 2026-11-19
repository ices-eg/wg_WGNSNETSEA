library(pacman)
p_load(data.table, ggplot2, mapplots, shapefiles, plyr, RANN, icesDatras, sf)

setwd("Q:/20-forskning/20-dfad/users/jostou/home/26_10_27_seaNet")

coast <- read.shapefile("data/europe")

cpue1 <- rbind.fill(read.csv("data/NS_IBTS_CPUE.csv"),
                   read.csv("data/BTS_CPUE.csv"),
                   read.csv("data/DYFS_CPUE.csv"))

ca <- rbind(getDATRAS(record = "CA", survey = "NS-IBTS", 2020:2025, c(1, 3)),
            getDATRAS(record = "CA", survey = "BTS", 2020:2025, c(1, 3)))


#ALK
ca$LngtClass[ca$LngtClass == -9] <- NA
ca <- ca[!is.na(ca$LngtClass), ]
ca[ca$LngtCode == ".", "LngtClass"] <- ca[ca$LngtCode == ".", "LngtClass"]/10
ca[ca$LngtCode == "0", "LngtClass"] <- ca[ca$LngtCode == "0", "LngtClass"]/10

ca <- ca[ca$SpecCode == "127140" & 
           !is.na(ca$Age) &
           !is.na(ca$LngtClass), ]


lm_Q1 <- nls(LngtClass ~ SSlogis(Age, Asym, xmid, scal), data = ca[ca$Quarter == 1, ])
lm_Q3 <- nls(LngtClass ~ SSlogis(Age, Asym, xmid, scal), data = ca[ca$Quarter == 3, ])

plot(ca[ca$Quarter == 1, ]$Age, ca[ca$Quarter == 1, ]$LngtClass, col = "white",
     main = "ALK Q1",
     ylab = "LngtClass",
     xlab = "Age",
     cex.main = 1.5,
     cex.lab = 1.5,
     xaxt="n")+
  axis(1,at=unique(ca[ca$Quarter == 1, ]$Age),labels=unique(ca[ca$Quarter == 1, ]$Age))+
  lines(sort(ca[ca$Quarter == 1, ]$Age), fitted(lm_Q1)[order(ca[ca$Quarter == 1, ]$Age)], col='red', lwd=1.5)+
  points(ca[ca$Quarter == 1, ]$Age, ca[ca$Quarter == 1, ]$LngtClass, pch = 16, cex = 1, col = "blue")


plot(ca[ca$Quarter == 3, ]$Age, ca[ca$Quarter == 3, ]$LngtClass, col = "white",
     main = "ALK Q3",
     ylab = "LngtClass",
     xlab = "Age",
     cex.main = 1.5,
     cex.lab = 1.5,
     xaxt="n")+
  axis(1,at=unique(ca[ca$Quarter == 3, ]$Age),labels=unique(ca[ca$Quarter == 3, ]$Age))+
  lines(sort(ca[ca$Quarter == 3, ]$Age), fitted(lm_Q3)[order(ca[ca$Quarter == 3, ]$Age)], col='red', lwd=1.5)+
  points(ca[ca$Quarter == 3, ]$Age, ca[ca$Quarter == 3, ]$LngtClass, pch = 16, cex = 1, col = "blue")



#Apply

yrs <- 2025
cruises <- c("NS-IBTS", "BTS", "DYFS")
qtr <- c(1, 3)
ages <- c(0, 1, 2)

for(yr in yrs){
  
  for (cruise in cruises) {
    
    for (q in qtr){
      
     
      cpue <- cpue1[cpue1$Survey == cruise &
                      cpue1$Year == yr &
                      cpue1$Quarter == q, ]
      
      if (nrow(cpue) == 0)
        next
      
      cpue$val <- cpue$CPUE_number_per_hour
      nam <- "CPUE_number_per_hour"
      if (cruise %in% c("BTS", "DYFS")){
        cpue$val <- cpue$CPUE_number_per_km2
        nam <- "CPUE_number_per_km2"
      }
        
      
      #apply
      cpue$LngtClas <- cpue$LngtClas/10
      
      lm <- lm_Q1
      if (q == 3)
        lm <- lm_Q3
      
      pred <- function(x) predict(lm, list(Age = x))
      
      age <- numeric()
      for (y in cpue$LngtClas) {
        x <- optimize(function(x) (pred(x) - y)^2, range(ca$Age, na.rm = T))
        age <- c(age, round(x[["minimum"]]))
      }
      
      cpue$age <- age
      cpue$station <- paste0(cpue$Year, cpue$Quarter, cpue$Ship, cpue$HaulNo)
      
      for (a in c(0, 1, 2)) {
        
        dat <- data.table(cpue[cpue$age == a & 
                                 cpue$val > 0, ])
        
        if (nrow(dat) == 0)
          next
        
        zero <- cpue[! cpue$station %in% dat$station, ]
        
        
        dat <- dat[ ,. (cpue = sum(val)),
                    by = .(Year, Quarter, Ship, HaulNo, ShootLat, ShootLong)]
        
        #limit of the plot
        xlim <- c(-4, 12)
        ylim <- c(49, 62)
        
        # add a point to be used as a sort of label....
        dat <- rbind.fill(dat, data.frame(ShootLong = rep(9, 4), 
                                          ShootLat = c(49.70, 50, 50.4, 51), 
                                          cpue = quantile(dat$cpue)[-1]))
        
        #plot the plot
        basemap(xlim, ylim, main = paste0(yr, " ", cruise, " Q", q, " age ", a),
                xaxt = "n", yaxt = "n", xlab="", ylab="")
        draw.rect()
        draw.shape(coast, col="cornsilk")
        rect(7.5, 49.4, 11.5, 52, density = NULL, angle = 45,
             col = "lightblue", border = NULL, lty = par("lty"), lwd = par("lwd"))
        draw.bubble(dat$ShootLong, dat$ShootLat, 
                    dat$cpue, maxradius = 0.8, pch=21, bg="#ABABAB95")
        axis(1, at = seq(-4, 12, 1), 
             labels = c(paste0(seq(4,1, by = -1), "°W"),
                        "0°", paste0(seq(1, 12, 1), "°E")))
        axis(2, at = seq(49, 62, 1), 
             labels = paste0(seq(49, 62, 1), "°N"))
        text(x = c(7.7, rep(10, 4)), y = c(51.8, 49.8, 50.15, 50.52, 51.1), 
             labels = c(nam, round(quantile(dat$cpue)[-1], 2)),
             adj = c(0,1))
        draw.bubble(zero$ShootLong, zero$ShootLat,
                    1, maxradius = 0.05, pch=4, bg="red")
        
        #export the plot
        dev.copy(png, paste0("lemon/cpue/",yr, " ", cruise, " Q", q, " age ", a, ".png"), 
                 width = 180, height = 240, units = 'mm', res = 300)
        dev.off()
        
      }
    }
  }
}









