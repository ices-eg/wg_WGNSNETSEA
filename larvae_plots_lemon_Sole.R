library(pacman)
p_load(data.table, ggplot2, mapplots, shapefiles, plyr)

setwd("Q:/20-forskning/20-dfad/users/jostou/home/26_10_27_seaNet")

coast <- read.shapefile("data/europe")

lengths <- read.csv("data/EggsAndLarvae_EM_2026_10_27_12_37_50.csv")
station <- read.csv("data/EggsAndLarvae_EH_2026_10_27_12_37_50.csv")

lengths$Length <- round(as.numeric(lengths$Length))
lengths$Number <- as.numeric(lengths$Number)

lengths <- lengths[lengths$Species == "Microstomus kitt" &
                     ! is.na(lengths$Length) &
                     lengths$Length > 0 &
                     lengths$Length < 50, ]

##
setDT(lengths)
len <- lengths[ ,. (Number = sum(Number)),
                by = .(Length)]

ggplot(len, aes(x=Length, y=Number)) +
  geom_bar(stat="identity")+
  theme_bw()


lengths$cat <- cut(lengths$Length, c(0, 6, 10, 14, 18, 24, 29, 51),
                   labels = c("0-7", "7-10", "11-14", "15-18", 
                              "19-24", "25-29", "29-50"))

setDT(lengths)
len <- lengths[ ,. (Number = sum(Number)),
                by = .(cat)]

ggplot(len, aes(x=cat, y=Number)) +
  geom_bar(stat="identity")+
  theme_bw()


#############################################

dat <- merge(lengths, station, by = "ICES_HaulID", all.y = T)

dat$Length[is.na(dat$Length)] <- 0

dat$cat <- as.character(dat$cat)
dat$cat[is.na(dat$cat)] <- 0

##
dat$ELVolFlag[dat$ELVolFlag == ""] <- "F"

dat$NumberLarvae <- dat$Number*dat$RaisingFactor
dat$L.sqm<- ifelse(dat$ELVolFlag %in% c("f", "F"), 
                   dat$NumberLarvae*dat$DepthLower*dat$FlowIntCalibr/(dat$FlowIntRevs*dat$NetopeningArea), 
                   ifelse(dat$ELVolFlag %in% c("d", "D"), dat$NumberLarvae*dat$DepthLower/(dat$Distance*dat$NetopeningArea), 
                          dat$NumberLarvae*dat$DepthLower/(3002*dat$NetopeningArea)))

dat$L.sqm[is.na(dat$L.sqm)] <- 0

#############################################

dat$survey <- ifelse(dat$Survey == "I4537", "MIK-IBTS", "MIK-Downs")

dat$date <- as.Date(paste0(dat$Yea, "-", 
                           sprintf("%02d", dat$Month), "-",
                           sprintf("%02d", dat$Day)))


dat$quarter <- quarter(dat$date)

setDT(dat)
len <- dat[ ,. (L.sqm = sum(L.sqm)),
            by = .(survey, quarter, ICES_HaulID, cat, StartLatitude, StartLongitude,
                   Year, Month, Day)]

len2 <- dat[ ,. (L.sqm = sum(L.sqm)),
             by = .(survey, quarter, ICES_HaulID, StartLatitude, StartLongitude,
                    Year, Month, Day)]


len2$cat <- "All"
len <- rbind(len, len2)

#limit of the plot
xlim <- c(-4, 12)
ylim <- c(49, 62)


for (yr in 2024:2025) {
  
  xx <- len[len$Year == yr, ]
  
  for (cruise in unique(xx$survey)) { #unique(xx$cat)[unique(xx$cat) != 0]
    
    dat2 <- xx[xx$cat == "All" &
                 xx$survey == cruise &
                 xx$cat != 0, ]
    
    zero <- xx[xx$cat == 0 &
                 xx$survey == cruise, ]
    
    nam <- "nr pr. m2"
    
    # if (category == "All"){
    #   dat2$L.sqm <- log10(1+ dat2$L.sqm)
    #   nam <- "log(1 + nr pr. m2)"
    # }
    
    
    # add a point to be used as a sort of label....
    dat2 <- rbind.fill(dat2, data.frame(StartLongitude = rep(9, 4), 
                                        StartLatitude = c(49.70, 50, 50.4, 51), 
                                        L.sqm = quantile(dat2$L.sqm)[-1]))
    
    #plot the plot
    basemap(xlim, ylim, main = paste0("lengths ", cruise, " - ", yr, "- Q", dat2$quarter[1]),
            xaxt = "n", yaxt = "n", xlab="", ylab="")
    draw.rect()
    draw.shape(coast, col="cornsilk")
    rect(7.5, 49.4, 11.5, 52, density = NULL, angle = 45,
         col = "lightblue", border = NULL, lty = par("lty"), lwd = par("lwd"))
    draw.bubble(dat2$StartLongitude, dat2$StartLatitude, 
                dat2$L.sqm, maxradius = 0.8, pch=21, bg="#ABABAB95")
    axis(1, at = seq(-4, 12, 1), 
         labels = c(paste0(seq(4,1, by = -1), "°W"),
                    "0°", paste0(seq(1, 12, 1), "°E")))
    axis(2, at = seq(49, 62, 1), 
         labels = paste0(seq(49, 62, 1), "°N"))
    text(x = c(7.7, rep(10, 4)), y = c(51.8, 49.8, 50.15, 50.52, 51.1), 
         labels = c(nam, round(quantile(dat2$L.sqm)[-1], 2)),
         adj = c(0,1))
    draw.bubble(zero$StartLongitude, zero$StartLatitude,
                1, maxradius = 0.05, pch=4, bg="red")
    
    #export the plot
    dev.copy(png, paste0("lemon/by_lengths/", cruise, " - ", yr, "- Q", dat2$quarter[1], ".png"), 
             width = 180, height = 240, units = 'mm', res = 300)
    dev.off()
    
  }
}


################################################################################


#############################################

len <- data.table(dat[dat$Survey == "I4537" &
                        dat$Length >= 19, ])
len <- len[ ,. (L.sqm = sum(L.sqm)),
            by = .(ICES_HaulID, StartLatitude, StartLongitude,
                   Year, Month, Day)]

len$date <- as.Date(paste0(len$Year, "-", 
                           sprintf("%02d", len$Month), "-",
                           sprintf("%02d", len$Day)))

len$week <- week(len$date)

len$cat <- "All"

#limit of the plot
xlim <- c(-4, 12)
ylim <- c(49, 62)


for (yr in 2024:2024) {
  
  xx <- len[len$Year == yr, ]
  
  for (wk in unique(xx$week)) {
    
    dat2 <- xx[xx$week == wk &
                 xx$cat != 0, ]
    
    zero <- xx[xx$cat == 0, ]
    
    # add a point to be used as a sort of label....
    dat2 <- rbind.fill(dat2, data.frame(StartLongitude = rep(9, 4), 
                                        StartLatitude = c(49.70, 50, 50.4, 51), 
                                        L.sqm = quantile(dat2$L.sqm)[-1]))
    
    #plot the plot
    basemap(xlim, ylim, main = paste0("week ", wk, " - ", yr ),
            xaxt = "n", yaxt = "n", xlab="", ylab="")
    draw.rect()
    draw.shape(coast, col="cornsilk")
    rect(7.5, 49.4, 11.5, 52, density = NULL, angle = 45,
         col = "lightblue", border = NULL, lty = par("lty"), lwd = par("lwd"))
    draw.bubble(dat2$StartLongitude, dat2$StartLatitude, 
                dat2$L.sqm, maxradius = 0.8, pch=21, bg="#ABABAB95")
    axis(1, at = seq(-4, 12, 1), 
         labels = c(paste0(seq(4,1, by = -1), "°W"),
                    "0°", paste0(seq(1, 12, 1), "°E")))
    axis(2, at = seq(49, 62, 1), 
         labels = paste0(seq(49, 62, 1), "°N"))
    text(x = c(7.7, rep(10, 4)), y = c(51.8, 49.8, 50.15, 50.52, 51.1), 
         labels = c("pr m2", round(quantile(dat2$L.sqm)[-1], 2)),
         adj = c(0,1))
    # draw.bubble(zero$StartLongitude, zero$StartLatitude, 
    #             1, maxradius = 0.1, pch=21, bg="red")
    
    #export the plot
    dev.copy(png, paste0("lemon/by_weeks/", yr, "_", wk, ".png"), 
             width = 180, height = 240, units = 'mm', res = 300)
    dev.off()
    
  }
}




