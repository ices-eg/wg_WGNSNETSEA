

# Script for producing age 1 herring distribution maps from NS-IBTS surveys
# compiled by Patrik Börjesson during WGNSNETSEA 2026-10-29. 


library(DATRASextra) # for downloading and data cleaning
library(DATRAS) # for ALK and plotting
library(maps)
library(mapdata)


# download survey data from DATRAS (done once)
download_datras(path = "data/datras", surveys = "NS-IBTS", years = 2023:2025)

# load survey data
x <- read_datras("data/datras", surveys = "NS-IBTS", years = 2023:2025)

# clean and subset herring data
herQ1Q3 <- clean_datras(x, years = 2023:2025, quarters = c(1,3),
                        aphias = 126417, verbose = FALSE)

# add swept area. could have used haul duration as used for index calculation
herQ1Q3 <- add_swept_area(herQ1Q3, width = "WingSpread",
                          impute_missing = TRUE, full_report = TRUE)

# moving to DATRAS package & add spectrum
d <- herQ1Q3
d <- addSpectrum(d)

# separate quarters
dQ1 <- subset(d, Quarter == 1)
dQ3 <- subset(d, Quarter == 3)

# add ages by quarters
dQ1 <- addNage(dQ1,1:6)
dQ3 <- addNage(dQ3,0:6)

# drop length bins & age 0
dQ1$N <- NULL
dQ3$N <- NULL
dQ3$Nage <-  dQ3$Nage[, colnames(dQ3$Nage) != "0", drop = FALSE]
colnames(dQ3$Nage)

# combine Q1 & Q3 for plotting
dd <- c(dQ1, dQ3)

# calculate cpue for one-year old, again - use haul duration
dd$cpue <- dd$Nage[,1] / dd$SweptArea 


# plot cpue over years and quarters 
for (yy in 2023:2025) {
  for (qq in c(1, 3)) {
    
    jpeg(
      filename = paste0('output_folder/', yy, " Q", qq, " aqe 1.jpg"),
      width = 10, height = 10, units = "in", res = 300, quality = 95
    )
    
    tmp <- subset(dd, Year == yy & Quarter == qq)
    
    bubblePlot(
      tmp,
      response = "cpue",
      scale = 8/max(sqrt(tmp[[2]]$cpue), na.rm = TRUE), # default 4/max(sqrt(d[[2]]$resp.var)
      col.zero = "black",
      pch.zero = "+",
      rim = TRUE,
      col = "green"
    )
    dev.off()
  }
}





