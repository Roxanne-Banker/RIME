


rm(list = ls())

setwd("/home/rbanker/Desktop/RIIME/")

#load("~/Desktop/RIIME/Bayesian-analyses/best-guild-2.RData")
#load("~/Desktop/RIIME/best-guild-backup.RData")

save.image("~/Desktop/RIIME/EPN-final.RData")
#save.image("~/Desktop/RIIME/best-guild-backup.RData")


# load packages
# library(tidyr)
library(dplyr)
library(rjags)
library(BEST)


set.seed(3645)


# I. Data import and manipulation ----

all.data <- read.csv("Datasets/RIME_all_collapsed_data-R-final.csv")
#names(all.data)



# [1] "guild"            "functional_group" "guild_no"         "ntp"              "mean_size"        "An.C"            
# [7] "C.B"              "B.Ap" 

anisian   <- all.data[,c(1,2,3,4,5,6,9,28:30)]
carnian   <- all.data[,c(1,2,3,10,11,12,15,28:30)]
bathonian <- all.data[,c(1,2,3,16,17,18,21,28:30)]
aptian    <- all.data[,c(1,2,3,22,23,24,27,28:30)]
names(anisian)

names(anisian)   <- c("guild", "functional_group", "guild_no", "g_no_sp","ntp", "mean_size", "lc.mean","An.C", "C.B", "B.Ap")
names(carnian)   <- c("guild", "functional_group", "guild_no", "g_no_sp","ntp", "mean_size", "lc.mean","An.C", "C.B", "B.Ap")
names(bathonian) <- c("guild", "functional_group", "guild_no", "g_no_sp","ntp", "mean_size", "lc.mean","An.C", "C.B", "B.Ap")
names(aptian)    <- c("guild", "functional_group", "guild_no", "g_no_sp","ntp", "mean_size", "lc.mean","An.C", "C.B", "B.Ap")

anisian   <- anisian %>% mutate(no_sp_pa = g_no_sp /sum(g_no_sp) * 100 )
carnian   <- carnian %>% mutate(no_sp_pa = g_no_sp /sum(g_no_sp) * 100 )
bathonian <- bathonian %>% mutate(no_sp_pa = g_no_sp /sum(g_no_sp) * 100 )
aptian    <- aptian %>% mutate(no_sp_pa = g_no_sp /sum(g_no_sp) * 100 )

anisian   <- anisian %>% mutate(log_mean_size   = log10(mean_size) )
carnian   <- carnian %>% mutate(log_mean_size   = log10(mean_size) )
bathonian <- bathonian %>% mutate(log_mean_size = log10(mean_size) )
aptian    <- aptian %>% mutate(log_mean_size    = log10(mean_size) )


# II. Guild NTP ----

## A. Anisian - Carnian ----

## Extinct vs Persistent ===

# Extinct: Anisian ntps
evp.extinct.An <- filter(anisian, An.C == "extinct")$ntp %>% .[.>2]
# Persistent: Anisian ntps
evp.persist.An <- filter(anisian, An.C == "persist")$ntp %>% .[.>2]
# Bayesian comparison of groups
AnC.evp.BEST.out <- BESTmcmc(evp.extinct.An, evp.persist.An)
# save(AnC.evp.BEST.out, file = "AnC-evp-BEST-out.RData")

#plot(AnC.evp.BEST.out)

## Persistent vs Persistent ===

# Persistent: Anisian ntps
pvp.persist.An <- filter(anisian, An.C == "persist")$ntp %>% .[.>2]
# Persistent: Carnian ntps
pvp.persist.C <- filter(carnian, An.C == "persist")$ntp %>% .[.>2]
# Bayesian comparison of groups
AnC.pvp.BEST.out <- BESTmcmc(pvp.persist.An, pvp.persist.C)
# save(AnC.pvp.BEST.out, file = "AnC-pvp-BEST-out.RData")
summary(AnC.pvp.BEST.out)


## New vs Persistent ===

# New: Carnian ntps
nvp.new.C <- filter(carnian, An.C == "new")$ntp %>% .[.>2]
# Persistent: Carnian ntps
nvp.persist.C <- filter(carnian, An.C == "persist")$ntp %>% .[.>2]
# Bayesian comparison of groups
AnC.nvp.BEST.out <- BESTmcmc(nvp.new.C, nvp.persist.C)
# save(AnC.nvp.BEST.out, file = "AnC-nvp-BEST-out.RData")



## Extinct vs New ===

# Extinct: Anisian ntps
evn.extinct.An <- filter(anisian, An.C == "extinct")$ntp %>% .[.>2]
# New: Carnian ntps
evn.new.C <- filter(carnian, An.C == "new")$ntp %>% .[.>2]
# Bayesian comparison of groups
AnC.evn.BEST.out <- BESTmcmc(evn.extinct.An, evn.new.C)
# save(AnC.evn.BEST.out, file = "AnC-evn-BEST-out.RData")





## B. Carnian - Bathonian  ----

## Extinct vs Persistent ===

# Extinct: Carnian ntps
evp.extinct.C <- filter(carnian, C.B == "extinct")$ntp %>% .[.>2]
# Persistent: Carnian ntps
evp.persist.C <- filter(carnian, C.B == "persist")$ntp %>% .[.>2]
# Bayesian comparison of groups
CB.evp.BEST.out <- BESTmcmc(evp.extinct.C, evp.persist.C)
# save(CB.evp.BEST.out, file = "CB-evp-BEST-out.RData")




## Persistent vs Persistent ===

# Persistent: Carnian ntps
pvp.persist.C2 <- filter(carnian, C.B == "persist")$ntp %>% .[.>2]
# Persistent: Bathonian ntps
pvp.persist.B <- filter(bathonian, C.B == "persist")$ntp %>% .[.>2]
# Bayesian comparison of groups
CB.pvp.BEST.out <- BESTmcmc(pvp.persist.C2, pvp.persist.B)
# save(CB.pvp.BEST.out, file = "CB-pvp-BEST-out.RData")



## New vs Persistent ===

# New: Bathonian ntps
nvp.new.B <- filter(bathonian, C.B == "new")$ntp %>% .[.>2]
# Persistent: Bathonian ntps
nvp.persist.B <- filter(bathonian, C.B == "persist")$ntp %>% .[.>2]
# Bayesian comparison of groups
CB.nvp.BEST.out <- BESTmcmc(nvp.new.B, nvp.persist.B)
# save(CB.nvp.BEST.out, file = "CB-nvp-BEST-out.RData")


## Extinct vs New ===

# Extinct: Carnian ntps
evn.extinct.C <- filter(carnian, C.B == "extinct")$ntp %>% .[.>2]
# New: Bathonian ntps
evn.new.B <- filter(bathonian, C.B == "new")$ntp %>% .[.>2]
# Bayesian comparison of groups
CB.evn.BEST.out <- BESTmcmc(evn.extinct.C, evn.new.B)
# save(CB.evn.BEST.out, file = "CB-evn-BEST-out.RData")


## C. Bathonian - Aptian  ----

## Extinct vs Persistent ===

# Extinct: Bathonian ntps
evp.extinct.B <- filter(bathonian, B.Ap == "extinct")$ntp %>% .[.>2]
# Persistent: Bathonian ntps
evp.persist.B <- filter(bathonian, B.Ap == "persist")$ntp %>% .[.>2]
# Bayesian comparison of groups
BAp.evp.BEST.out <- BESTmcmc(evp.extinct.B, evp.persist.B)
# save(BAp.evp.BEST.out, file = "BAp-evp-BEST-out.RData")



## Persistent vs Persistent ===

# Persistent: Bathonian ntps
pvp.persist.B2 <- filter(bathonian, B.Ap == "persist")$ntp %>% .[.>2]
# Persistent: Aptian ntps
pvp.persist.Ap <- filter(aptian, B.Ap == "persist")$ntp %>% .[.>2]
# Bayesian comparison of groups
BAp.pvp.BEST.out <- BESTmcmc(pvp.persist.B2, pvp.persist.Ap)
# save(BAp.pvp.BEST.out, file = "BAp-pvp-BEST-out.RData")



## New vs Persistent ===

# New: Aptian ntps
nvp.new.Ap <- filter(aptian, B.Ap == "new")$ntp %>% .[.>2]
# Persistent: Aptian ntps
nvp.persist.Ap <- filter(aptian, B.Ap == "persist")$ntp %>% .[.>2]
# Bayesian comparison of groups
BAp.nvp.BEST.out <- BESTmcmc(nvp.new.Ap, nvp.persist.Ap)
# save(BAp.nvp.BEST.out, file = "BAp-nvp-BEST-out.RData")


## Extinct vs New ===

# Extinct: Bathonian ntps
evn.extinct.B <- filter(bathonian, B.Ap == "extinct")$ntp %>% .[.>2]
# New: Aptian ntps
evn.new.Ap <- filter(aptian, B.Ap == "new")$ntp %>% .[.>2]
# Bayesian comparison of groups
BAp.evn.BEST.out <- BESTmcmc(evn.extinct.B, evn.new.Ap)
# save(BAp.evn.BEST.out, file = "BAp-evn-BEST-out.RData")





## D. NTP Figure ----


# "#9FDA3AFF"
# "#35B779FF"



{par(mfrow = c(3, 4))
  
plot(BAp.evp.BEST.out, 
     showCurve = TRUE,
     credMass = 0.95,
     compVal = 0,
     comparisonColor = "black",
     xlab = expression(mu["1 [Extinct Bathonian]"]-mu["2 [Persistent Bathonian]"]),
     col= "#35B779FF")
plot(BAp.pvp.BEST.out, 
     showCurve = TRUE,
     credMass = 0.95,
     compVal = 0,
     comparisonColor = "black",
     xlab = expression(mu["1 [Persistent Bathonian]"]-mu["2 [Persistent Aptian]"]),
     col= "gray")
plot(BAp.nvp.BEST.out, 
     showCurve = TRUE,
     credMass = 0.95,
     compVal = 0,
     comparisonColor = "black",
     xlab = expression(mu["1 [New Aptian]"]-mu["2 [Persistent Aptian]"]),
     col= "#35B779FF")
plot(BAp.evn.BEST.out, 
     showCurve = TRUE,
     credMass = 0.95,
     compVal = 0,
     comparisonColor = "black",
     xlab = expression(mu["1 [Extinct Bathonian]"]-mu["2 [New Aptian]"]),
     col= "gray")

plot(CB.evp.BEST.out, 
     showCurve = TRUE,
     credMass = 0.95,
     compVal = 0,
     comparisonColor = "black",
     xlab = expression(mu["1 [Extinct Carnian]"]-mu["2 [Persistent Carnian]"]),
     col= "#35B779FF")
plot(CB.pvp.BEST.out, 
     showCurve = TRUE,
     credMass = 0.95,
     compVal = 0,
     comparisonColor = "black",
     xlab = expression(mu["1 [Persistent Carnian]"]-mu["2 [Persistent Bathonian]"]),
     col= "gray")
plot(CB.nvp.BEST.out, 
     showCurve = TRUE,
     credMass = 0.95,
     compVal = 0,
     comparisonColor = "black",
     xlab = expression(mu["1 [New Bathonian]"]-mu["2 [Persistent Bathonian]"]),
     col= "#35B779FF")
plot(CB.evn.BEST.out, 
     showCurve = TRUE,
     credMass = 0.95,
     compVal = 0,
     comparisonColor = "black",
     xlab = expression(mu["1 [Extinct Carnian]"]-mu["2 [New Bathonian]"]),
     col= "gray")

plot(AnC.evp.BEST.out, 
     showCurve = TRUE,
     credMass = 0.95,
     compVal = 0,
     comparisonColor = "black",
     xlab = expression(mu["1 [Extinct Anisian]"]-mu["2 [Persistent Anisian]"]),
     col= "#35B779FF")
plot(AnC.pvp.BEST.out, 
     showCurve = TRUE,
     credMass = 0.95,
     compVal = 0,
     comparisonColor = "black",
     xlab = expression(mu["1 [Persistent Anisian]"]-mu["2 [Persistent Carnian]"]),
     col= "gray")
plot(AnC.nvp.BEST.out, 
     showCurve = TRUE,
     credMass = 0.95,
     compVal = 0,
     comparisonColor = "black",
     xlab = expression(mu["1 [New Carnian]"]-mu["2 [Persistent Carnian]"]),
     col= "gray")
plot(AnC.evn.BEST.out, 
     showCurve = TRUE,
     credMass = 0.95,
     compVal = 0,
     comparisonColor = "black",
     xlab = expression(mu["1 [Extinct Anisian]"]-mu["2 [New Carnian]"]),
     col= "#35B779FF")

  }


# plot(AnC.evn.BEST.out,
#      showCurve = TRUE,
#      credMass = 0.95,
#      compVal = 0,
#      comparisonColor = "black",
#      xlab = expression(paste(mu["1 [Extinct Anisian] "],"\n-",mu["2 [New Carnian]"] ) ),
#      col= "#35B779FF")
# 
# plot(AnC.evn.BEST.out,
#      showCurve = TRUE,
#      credMass = 0.95,
#      compVal = 0,
#      comparisonColor = "black",
#      xlab = expression(atop(mu["1t"], Bootstrap~samples*','~Allianz) ) ,
#      col= "#35B779FF")
# 
# expression(atop("Histogram of "*hat(mu), Bootstrap~samples*','~Allianz))
# 
# expression(mu["1 [Extinct Anisian]"]-mu["2 [New Carnian]"]),





# III. Guild Body Size ----

## A. Anisian - Carnian ----

## Extinct vs Persistent ===

# Extinct: Anisian mean_sizes
evp.extinct.An.bs <- filter(anisian, An.C == "extinct")$mean_size 
evp.extinct.An.bs <- evp.extinct.An.bs[!is.na(evp.extinct.An.bs)]
# Persistent: Anisian mean_sizes
evp.persist.An.bs <- filter(anisian, An.C == "persist")$mean_size 
evp.persist.An.bs <- evp.persist.An.bs[!is.na(evp.persist.An.bs)]
# Bayesian comparison of groups
AnC.evp.BEST.out.bs <- BESTmcmc(evp.extinct.An.bs, evp.persist.An.bs)
# save(AnC.evp.BEST.out.bs, file = "AnC-evp-BEST-out.RData")



## Persistent vs Persistent ===

# Persistent: Anisian mean_sizes
pvp.persist.An.bs <- filter(anisian, An.C == "persist")$mean_size 
pvp.persist.An.bs <- pvp.persist.An.bs[!is.na(pvp.persist.An.bs)]
# Persistent: Carnian mean_sizes
pvp.persist.C.bs <- filter(carnian, An.C == "persist")$mean_size 
pvp.persist.C.bs <- pvp.persist.C.bs[!is.na(pvp.persist.C.bs)]
# Bayesian comparison of groups
AnC.pvp.BEST.out.bs <- BESTmcmc(pvp.persist.An.bs, pvp.persist.C.bs)
# save(AnC.pvp.BEST.out.bs, file = "AnC-pvp-BEST-out.RData")



## New vs Persistent ===

# New: Carnian mean_sizes
nvp.new.C.bs <- filter(carnian, An.C == "new")$mean_size 
nvp.new.C.bs <- nvp.new.C.bs[!is.na(nvp.new.C.bs)]
# Persistent: Carnian mean_sizes
nvp.persist.C.bs <- filter(carnian, An.C == "persist")$mean_size 
nvp.persist.C.bs <- nvp.persist.C.bs[!is.na(nvp.persist.C.bs)]
# Bayesian comparison of groups
AnC.nvp.BEST.out.bs <- BESTmcmc(nvp.new.C.bs, nvp.persist.C.bs)
# save(AnC.nvp.BEST.out.bs, file = "AnC-nvp-BEST-out.RData")



## Extinct vs New ===

# Extinct: Anisian mean_sizes
evn.extinct.An.bs <- filter(anisian, An.C == "extinct")$mean_size 
evn.extinct.An.bs <- evn.extinct.An.bs[!is.na(evn.extinct.An.bs)]
# New: Carnian mean_sizes
evn.new.C.bs <- filter(carnian, An.C == "new")$mean_size 
evn.new.C.bs <- evn.new.C.bs[!is.na(evn.new.C.bs)]
# Bayesian comparison of groups
AnC.evn.BEST.out.bs <- BESTmcmc(evn.extinct.An.bs, evn.new.C.bs)
# save(AnC.evn.BEST.out.bs, file = "AnC-evn-BEST-out.RData")





## B. Carnian - Bathonian  ----

## Extinct vs Persistent ===

# Extinct: Carnian mean_sizes
evp.extinct.C.bs <- filter(carnian, C.B == "extinct")$mean_size 
evp.extinct.C.bs <- evp.extinct.C.bs[!is.na(evp.extinct.C.bs)]
# Persistent: Carnian mean_sizes
evp.persist.C.bs <- filter(carnian, C.B == "persist")$mean_size 
evp.persist.C.bs <- evp.persist.C.bs[!is.na(evp.persist.C.bs)]
# Bayesian comparison of groups
CB.evp.BEST.out.bs <- BESTmcmc(evp.extinct.C.bs, evp.persist.C.bs)
# save(CB.evp.BEST.out.bs, file = "CB-evp-BEST-out.RData")



## Persistent vs Persistent ===

# Persistent: Carnian mean_sizes
pvp.persist.C2.bs <- filter(carnian, C.B == "persist")$mean_size 
pvp.persist.C2.bs <- pvp.persist.C2.bs[!is.na(pvp.persist.C2.bs)]
# Persistent: Bathonian mean_sizes
pvp.persist.B.bs <- filter(bathonian, C.B == "persist")$mean_size 
pvp.persist.B.bs <- pvp.persist.B.bs[!is.na(pvp.persist.B.bs)]
# Bayesian comparison of groups
CB.pvp.BEST.out.bs <- BESTmcmc(pvp.persist.C2.bs, pvp.persist.B.bs)
# save(CB.pvp.BEST.out.bs, file = "CB-pvp-BEST-out.RData")



## New vs Persistent ===

# New: Bathonian mean_sizes
nvp.new.B.bs <- filter(bathonian, C.B == "new")$mean_size 
nvp.new.B.bs <- nvp.new.B.bs[!is.na(nvp.new.B.bs)]
# Persistent: Bathonian mean_sizes
nvp.persist.B.bs <- filter(bathonian, C.B == "persist")$mean_size 
nvp.persist.B.bs <- nvp.persist.B.bs[!is.na(nvp.persist.B.bs)]
# Bayesian comparison of groups
CB.nvp.BEST.out.bs <- BESTmcmc(nvp.new.B.bs, nvp.persist.B.bs)
# save(CB.nvp.BEST.out.bs, file = "CB-nvp-BEST-out.RData")



## Extinct vs New ===

# Extinct: Carnian mean_sizes
evn.extinct.C.bs <- filter(carnian, C.B == "extinct")$mean_size 
evn.extinct.C.bs <- evn.extinct.C.bs[!is.na(evn.extinct.C.bs)]
# New: Bathonian mean_sizes
evn.new.B.bs <- filter(bathonian, C.B == "new")$mean_size 
evn.new.B.bs <- evn.new.B.bs[!is.na(evn.new.B.bs)]
# Bayesian comparison of groups
CB.evn.BEST.out.bs <- BESTmcmc(evn.extinct.C.bs, evn.new.B.bs)
# save(CB.evn.BEST.out.bs, file = "CB-evn-BEST-out.RData")




## C. Bathonian - Aptian  ----

## Extinct vs Persistent ===

# Extinct: Bathonian mean_sizes
evp.extinct.B.bs <- filter(bathonian, B.Ap == "extinct")$mean_size 
evp.extinct.B.bs <- evp.extinct.B.bs[!is.na(evp.extinct.B.bs)]
# Persistent: Bathonian mean_sizes
evp.persist.B.bs <- filter(bathonian, B.Ap == "persist")$mean_size 
evp.persist.B.bs <- evp.persist.B.bs[!is.na(evp.persist.B.bs)]
# Bayesian comparison of groups
BAp.evp.BEST.out.bs <- BESTmcmc(evp.extinct.B.bs, evp.persist.B.bs)
# save(BAp.evp.BEST.out.bs, file = "BAp-evp-BEST-out.RData")



## Persistent vs Persistent ===

# Persistent: Bathonian mean_sizes
pvp.persist.B2.bs <- filter(bathonian, B.Ap == "persist")$mean_size 
pvp.persist.B2.bs <- pvp.persist.B2.bs[!is.na(pvp.persist.B2.bs)]
# Persistent: Aptian mean_sizes
pvp.persist.Ap.bs <- filter(aptian, B.Ap == "persist")$mean_size 
pvp.persist.Ap.bs <- pvp.persist.Ap.bs[!is.na(pvp.persist.Ap.bs)]
# Bayesian comparison of groups
BAp.pvp.BEST.out.bs <- BESTmcmc(pvp.persist.B2.bs, pvp.persist.Ap.bs)
# save(BAp.pvp.BEST.out.bs, file = "BAp-pvp-BEST-out.RData")



## New vs Persistent ===

# New: Aptian mean_sizes
nvp.new.Ap.bs <- filter(aptian, B.Ap == "new")$mean_size 
nvp.new.Ap.bs <- nvp.new.Ap.bs[!is.na(nvp.new.Ap.bs)]
# Persistent: Aptian mean_sizes
nvp.persist.Ap.bs <- filter(aptian, B.Ap == "persist")$mean_size 
nvp.persist.Ap.bs <- nvp.persist.Ap.bs[!is.na(nvp.persist.Ap.bs)]
# Bayesian comparison of groups
BAp.nvp.BEST.out.bs <- BESTmcmc(nvp.new.Ap.bs, nvp.persist.Ap.bs)
# save(BAp.nvp.BEST.out.bs, file = "BAp-nvp-BEST-out.RData")



## Extinct vs New ===

# Extinct: Bathonian mean_sizes
evn.extinct.B.bs <- filter(bathonian, B.Ap == "extinct")$mean_size 
evn.extinct.B.bs <- evn.extinct.B.bs[!is.na(evn.extinct.B.bs)]
# New: Aptian mean_sizes
evn.new.Ap.bs <- filter(aptian, B.Ap == "new")$mean_size 
evn.new.Ap.bs <- evn.new.Ap.bs[!is.na(evn.new.Ap.bs)]
# Bayesian comparison of groups
BAp.evn.BEST.out.bs <- BESTmcmc(evn.extinct.B.bs, evn.new.Ap.bs)
# save(BAp.evn.BEST.out.bs, file = "BAp-evn-BEST-out.RData")




## D. Guild BS Fig ----


{par(mfrow = c(3, 4))
  
  plot(BAp.evp.BEST.out.bs, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [Extinct Bathonian]"]-mu["2 [Persistent Bathonian]"]),
       col= "gray")
  plot(BAp.pvp.BEST.out.bs, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [Persistent Bathonian]"]-mu["2 [Persistent Aptian]"]),
       col= "gray")
  plot(BAp.nvp.BEST.out.bs, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [New Aptian]"]-mu["2 [Persistent Aptian]"]),
       col= "gray")
  plot(BAp.evn.BEST.out.bs, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [Extinct Bathonian]"]-mu["2 [New Aptian]"]),
       col= "gray")
  
  plot(CB.evp.BEST.out.bs, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [Extinct Carnian]"]-mu["2 [Persistent Carnian]"]),
       col= "gray")
  plot(CB.pvp.BEST.out.bs, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [Persistent Carnian]"]-mu["2 [Persistent Bathonian]"]),
       col= "gray")
  plot(CB.nvp.BEST.out.bs, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [New Bathonian]"]-mu["2 [Persistent Bathonian]"]),
       col= "gray")
  plot(CB.evn.BEST.out.bs, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [Extinct Carnian]"]-mu["2 [New Bathonian]"]),
       col= "#35B779FF")
  
  plot(AnC.evp.BEST.out.bs, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [Extinct Anisian]"]-mu["2 [Persistent Anisian]"]),
       col= "gray")
  plot(AnC.pvp.BEST.out.bs, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [Persistent Anisian]"]-mu["2 [Persistent Carnian]"]),
       col= "gray")
  plot(AnC.nvp.BEST.out.bs, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [New Carnian]"]-mu["2 [Persistent Carnian]"]),
       col= "gray")
  plot(AnC.evn.BEST.out.bs, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [Extinct Anisian]"]-mu["2 [New Carnian]"]),
       col= "#35B779FF")
  
}





# IV. Guild long chain ----

## A. Anisian - Carnian ----

## Extinct vs Persistent ===

# Extinct: Anisian 
evp.extinct.An.lc <- filter(anisian, An.C == "extinct")$lc.mean %>% .[.>0] %>% .[!is.na(.)]
# Persistent: Anisian 
evp.persist.An.lc <- filter(anisian, An.C == "persist")$lc.mean %>% .[.>0] %>% .[!is.na(.)]
# Bayesian comparison of groups
AnC.evp.BEST.out.lc <- BESTmcmc(evp.extinct.An.lc, evp.persist.An.lc)
# save(AnC.evp.BEST.out, file = "AnC-evp-BEST-out.RData")



## B. Persistent vs Persistent ===

# Persistent: Anisian 
pvp.persist.An.lc <- filter(anisian, An.C == "persist")$lc.mean %>% .[.>0] %>% .[!is.na(.)]
# Persistent: Carnian 
pvp.persist.C.lc <- filter(carnian, An.C == "persist")$lc.mean %>% .[.>0] %>% .[!is.na(.)]
# Bayesian comparison of groups
AnC.pvp.BEST.out.lc <- BESTmcmc(pvp.persist.An.lc, pvp.persist.C.lc)
# save(AnC.pvp.BEST.out, file = "AnC-pvp-BEST-out.RData")



## C. New vs Persistent ===

# New: Carnian 
nvp.new.C.lc <- filter(carnian, An.C == "new")$lc.mean %>% .[.>0] %>% .[!is.na(.)]
# Persistent: Carnian  
nvp.persist.C.lc <- filter(carnian, An.C == "persist")$lc.mean %>% .[.>0] %>% .[!is.na(.)]
# Bayesian comparison of groups
AnC.nvp.BEST.out.lc <- BESTmcmc(nvp.new.C.lc, nvp.persist.C.lc)
# save(AnC.nvp.BEST.out, file = "AnC-nvp-BEST-out.RData")




## D. Extinct vs New ===

# Extinct: Anisian 
evn.extinct.An.lc <- filter(anisian, An.C == "extinct")$lc.mean %>% .[.>0] %>% .[!is.na(.)]
# New: Carnian 
evn.new.C.lc <- filter(carnian, An.C == "new")$lc.mean %>% .[.>0] %>% .[!is.na(.)]
# Bayesian comparison of groups
AnC.evn.BEST.out.lc <- BESTmcmc(evn.extinct.An.lc, evn.new.C.lc)
# save(AnC.evn.BEST.out, file = "AnC-evn-BEST-out.RData")





## B. Carnian - Bathonian  ----

## Extinct vs Persistent ===

# Extinct: Carnian 
evp.extinct.C.lc <- filter(carnian, C.B == "extinct")$lc.mean %>% .[.>0] %>% .[!is.na(.)]
# Persistent: Carnian 
evp.persist.C.lc <- filter(carnian, C.B == "persist")$lc.mean %>% .[.>0] %>% .[!is.na(.)]
# Bayesian comparison of groups
CB.evp.BEST.out.lc <- BESTmcmc(evp.extinct.C.lc, evp.persist.C.lc)
# save(CB.evp.BEST.out, file = "CB-evp-BEST-out.RData")




## Persistent vs Persistent ===

# Persistent: Carnian 
pvp.persist.C2.lc <- filter(carnian, C.B == "persist")$lc.mean %>% .[.>0] %>% .[!is.na(.)]
# Persistent: Bathonian 
pvp.persist.B.lc <- filter(bathonian, C.B == "persist")$lc.mean %>% .[.>0] %>% .[!is.na(.)]
# Bayesian comparison of groups
CB.pvp.BEST.out.lc <- BESTmcmc(pvp.persist.C2.lc, pvp.persist.B.lc)
# save(CB.pvp.BEST.out, file = "CB-pvp-BEST-out.RData")



## New vs Persistent ===

# New: Bathonian 
nvp.new.B.lc <- filter(bathonian, C.B == "new")$lc.mean %>% .[.>0] %>% .[!is.na(.)]
# Persistent: Bathonian 
nvp.persist.B.lc <- filter(bathonian, C.B == "persist")$lc.mean %>% .[.>0] %>% .[!is.na(.)]
# Bayesian comparison of groups
CB.nvp.BEST.out.lc <- BESTmcmc(nvp.new.B.lc, nvp.persist.B.lc)
# save(CB.nvp.BEST.out, file = "CB-nvp-BEST-out.RData")


## Extinct vs New ===

# Extinct: Carnian 
evn.extinct.C.lc <- filter(carnian, C.B == "extinct")$lc.mean %>% .[.>0] %>% .[!is.na(.)]
# New: Bathonian 
evn.new.B.lc <- filter(bathonian, C.B == "new")$lc.mean %>% .[.>0] %>% .[!is.na(.)]
# Bayesian comparison of groups
CB.evn.BEST.out.lc <- BESTmcmc(evn.extinct.C.lc, evn.new.B.lc)
# save(CB.evn.BEST.out, file = "CB-evn-BEST-out.RData")





## C. Bathonian - Aptian  ----

## Extinct vs Persistent ===

# Extinct: Bathonian 
evp.extinct.B.lc <- filter(bathonian, B.Ap == "extinct")$lc.mean %>% .[.>0] %>% .[!is.na(.)]
# Persistent: Bathonian 
evp.persist.B.lc <- filter(bathonian, B.Ap == "persist")$lc.mean %>% .[.>0] %>% .[!is.na(.)]
# Bayesian comparison of groups
BAp.evp.BEST.out.lc <- BESTmcmc(evp.extinct.B.lc, evp.persist.B.lc)
# save(BAp.evp.BEST.out, file = "BAp-evp-BEST-out.RData")



## Persistent vs Persistent ===

# Persistent: Bathonian 
pvp.persist.B2.lc <- filter(bathonian, B.Ap == "persist")$lc.mean %>% .[.>0] %>% .[!is.na(.)]
# Persistent: Aptian 
pvp.persist.Ap.lc <- filter(aptian, B.Ap == "persist")$lc.mean %>% .[.>0] %>% .[!is.na(.)]
# Bayesian comparison of groups
BAp.pvp.BEST.out.lc <- BESTmcmc(pvp.persist.B2.lc, pvp.persist.Ap.lc)
# save(BAp.pvp.BEST.out, file = "BAp-pvp-BEST-out.RData")



## New vs Persistent ===

# New: Aptian 
nvp.new.Ap.lc <- filter(aptian, B.Ap == "new")$lc.mean %>% .[.>0] %>% .[!is.na(.)]
# Persistent: Aptian 
nvp.persist.Ap.lc <- filter(aptian, B.Ap == "persist")$lc.mean %>% .[.>0] %>% .[!is.na(.)]
# Bayesian comparison of groups
BAp.nvp.BEST.out.lc <- BESTmcmc(nvp.new.Ap.lc, nvp.persist.Ap.lc)
# save(BAp.nvp.BEST.out, file = "BAp-nvp-BEST-out.RData")




## Extinct vs New ===

# Extinct: Bathonian 
evn.extinct.B.lc <- filter(bathonian, B.Ap == "extinct")$lc.mean %>% .[.>0] %>% .[!is.na(.)]
# New: Aptian 
evn.new.Ap.lc <- filter(aptian, B.Ap == "new")$lc.mean %>% .[.>0] %>% .[!is.na(.)]
# Bayesian comparison of groups
BAp.evn.BEST.out.lc <- BESTmcmc(evn.extinct.B.lc, evn.new.Ap.lc)
# save(BAp.evn.BEST.out, file = "BAp-evn-BEST-out.RData")



## D. Guild long chain Fig ----


{par(mfrow = c(3, 4))
  
  plot(BAp.evp.BEST.out.lc, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [Extinct Bathonian]"]-mu["2 [Persistent Bathonian]"]),
       col= "#35B779FF")
  plot(BAp.pvp.BEST.out.lc, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [Persistent Bathonian]"]-mu["2 [Persistent Aptian]"]),
       col= "#35B779FF")
  plot(BAp.nvp.BEST.out.lc, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [New Aptian]"]-mu["2 [Persistent Aptian]"]),
       col= "gray")
  plot(BAp.evn.BEST.out.lc, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [Extinct Bathonian]"]-mu["2 [New Aptian]"]),
       col= "gray")
  
  plot(CB.evp.BEST.out.lc, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [Extinct Carnian]"]-mu["2 [Persistent Carnian]"]),
       col= "gray")
  plot(CB.pvp.BEST.out.lc, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [Persistent Carnian]"]-mu["2 [Persistent Bathonian]"]),
       col= "gray")
  plot(CB.nvp.BEST.out.lc, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [New Bathonian]"]-mu["2 [Persistent Bathonian]"]),
       col= "gray")
  plot(CB.evn.BEST.out.lc, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [Extinct Carnian]"]-mu["2 [New Bathonian]"]),
       col= "gray")
  
  plot(AnC.evp.BEST.out.lc, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [Extinct Anisian]"]-mu["2 [Persistent Anisian]"]),
       col= "gray")
  plot(AnC.pvp.BEST.out.lc, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [Persistent Anisian]"]-mu["2 [Persistent Carnian]"]),
       col= "gray")
  plot(AnC.nvp.BEST.out.lc, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [New Carnian]"]-mu["2 [Persistent Carnian]"]),
       col= "gray")
  plot(AnC.evn.BEST.out.lc, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [Extinct Anisian]"]-mu["2 [New Carnian]"]),
       col= "gray")
  
}





# V. Guild Richness ----

## A. Anisian - Carnian ----

## Extinct vs Persistent ===

# Extinct: Anisian ntps
evp.extinct.An.r <- filter(anisian, An.C == "extinct")$no_sp_pa %>% .[.>0]
# Persistent: Anisian ntps
evp.persist.An.r <- filter(anisian, An.C == "persist")$no_sp_pa %>% .[.>0]
# Bayesian comparison of groups
AnC.evp.BEST.out.r <- BESTmcmc(evp.extinct.An.r, evp.persist.An.r)
# save(AnC.evp.BEST.out, file = "AnC-evp-BEST-out.RData")



## Persistent vs Persistent ===

# Persistent: Anisian ntps
pvp.persist.An.r <- filter(anisian, An.C == "persist")$no_sp_pa %>% .[.>0]
# Persistent: Carnian ntps
pvp.persist.C.r <- filter(carnian, An.C == "persist")$no_sp_pa %>% .[.>0]
# Bayesian comparison of groups
AnC.pvp.BEST.out.r <- BESTmcmc(pvp.persist.An.r, pvp.persist.C.r)
# save(AnC.pvp.BEST.out, file = "AnC-pvp-BEST-out.RData")



## New vs Persistent ===

# New: Carnian ntps
nvp.new.C.r <- filter(carnian, An.C == "new")$no_sp_pa %>% .[.>0]
# Persistent: Carnian ntps
nvp.persist.C.r <- filter(carnian, An.C == "persist")$no_sp_pa %>% .[.>0]
# Bayesian comparison of groups
AnC.nvp.BEST.out.r <- BESTmcmc(nvp.new.C.r, nvp.persist.C.r)
# save(AnC.nvp.BEST.out, file = "AnC-nvp-BEST-out.RData")



## Extinct vs New ===

# Extinct: Anisian ntps
evn.extinct.An.r <- filter(anisian, An.C == "extinct")$no_sp_pa %>% .[.>0]
# New: Carnian ntps
evn.new.C.r <- filter(carnian, An.C == "new")$no_sp_pa %>% .[.>0]
# Bayesian comparison of groups
AnC.evn.BEST.out.r <- BESTmcmc(evn.extinct.An.r, evn.new.C.r)
# save(AnC.evn.BEST.out, file = "AnC-evn-BEST-out.RData")





## B. Carnian - Bathonian  ----

## Extinct vs Persistent ===

# Extinct: Carnian ntps
evp.extinct.C.r <- filter(carnian, C.B == "extinct")$no_sp_pa %>% .[.>0]
# Persistent: Carnian ntps
evp.persist.C.r <- filter(carnian, C.B == "persist")$no_sp_pa %>% .[.>0]
# Bayesian comparison of groups
CB.evp.BEST.out.r <- BESTmcmc(evp.extinct.C.r, evp.persist.C.r)
# save(CB.evp.BEST.out, file = "CB-evp-BEST-out.RData")




## Persistent vs Persistent ===

# Persistent: Carnian ntps
pvp.persist.C2.r <- filter(carnian, C.B == "persist")$no_sp_pa %>% .[.>0]
# Persistent: Bathonian ntps
pvp.persist.B.r <- filter(bathonian, C.B == "persist")$no_sp_pa %>% .[.>0]
# Bayesian comparison of groups
CB.pvp.BEST.out.r <- BESTmcmc(pvp.persist.C2.r, pvp.persist.B.r)
# save(CB.pvp.BEST.out, file = "CB-pvp-BEST-out.RData")



## New vs Persistent ===

# New: Bathonian ntps
nvp.new.B.r <- filter(bathonian, C.B == "new")$no_sp_pa %>% .[.>0]
# Persistent: Bathonian ntps
nvp.persist.B.r <- filter(bathonian, C.B == "persist")$no_sp_pa %>% .[.>0]
# Bayesian comparison of groups
CB.nvp.BEST.out.r <- BESTmcmc(nvp.new.B.r, nvp.persist.B.r)
# save(CB.nvp.BEST.out, file = "CB-nvp-BEST-out.RData")


## Extinct vs New ===

# Extinct: Carnian ntps
evn.extinct.C.r <- filter(carnian, C.B == "extinct")$no_sp_pa %>% .[.>0]
# New: Bathonian ntps
evn.new.B.r <- filter(bathonian, C.B == "new")$no_sp_pa %>% .[.>0]
# Bayesian comparison of groups
CB.evn.BEST.out.r <- BESTmcmc(evn.extinct.C.r, evn.new.B.r)
# save(CB.evn.BEST.out, file = "CB-evn-BEST-out.RData")






## C. Bathonian - Aptian  ----

## Extinct vs Persistent ===

# Extinct: Bathonian ntps
evp.extinct.B.r <- filter(bathonian, B.Ap == "extinct")$no_sp_pa %>% .[.>0]
# Persistent: Bathonian ntps
evp.persist.B.r <- filter(bathonian, B.Ap == "persist")$no_sp_pa %>% .[.>0]
# Bayesian comparison of groups
BAp.evp.BEST.out.r <- BESTmcmc(evp.extinct.B.r, evp.persist.B.r)
# save(BAp.evp.BEST.out, file = "BAp-evp-BEST-out.RData")



## Persistent vs Persistent ===

# Persistent: Bathonian ntps
pvp.persist.B2.r <- filter(bathonian, B.Ap == "persist")$no_sp_pa %>% .[.>0]
# Persistent: Aptian ntps
pvp.persist.Ap.r <- filter(aptian, B.Ap == "persist")$no_sp_pa %>% .[.>0]
# Bayesian comparison of groups
BAp.pvp.BEST.out.r <- BESTmcmc(pvp.persist.B2.r, pvp.persist.Ap.r)
# save(BAp.pvp.BEST.out, file = "BAp-pvp-BEST-out.RData")



## New vs Persistent ===

# New: Aptian ntps
nvp.new.Ap.r <- filter(aptian, B.Ap == "new")$no_sp_pa %>% .[.>0]
# Persistent: Aptian ntps
nvp.persist.Ap.r <- filter(aptian, B.Ap == "persist")$no_sp_pa %>% .[.>0]
# Bayesian comparison of groups
BAp.nvp.BEST.out.r <- BESTmcmc(nvp.new.Ap.r, nvp.persist.Ap.r)
# save(BAp.nvp.BEST.out, file = "BAp-nvp-BEST-out.RData")




## Extinct vs New ===

# Extinct: Bathonian ntps
evn.extinct.B.r <- filter(bathonian, B.Ap == "extinct")$no_sp_pa %>% .[.>0]
# New: Aptian ntps
evn.new.Ap.r <- filter(aptian, B.Ap == "new")$no_sp_pa %>% .[.>0]
# Bayesian comparison of groups
BAp.evn.BEST.out.r <- BESTmcmc(evn.extinct.B.r, evn.new.Ap.r)
# save(BAp.evn.BEST.out, file = "BAp-evn-BEST-out.RData")

## D. Guild percent richness Fig ----


{par(mfrow = c(3, 4))
  
  plot(BAp.evp.BEST.out.r, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [Extinct Bathonian]"]-mu["2 [Persistent Bathonian]"]),
       col= "#35B779FF")
  plot(BAp.pvp.BEST.out.r, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [Persistent Bathonian]"]-mu["2 [Persistent Aptian]"]),
       col= "gray")
  plot(BAp.nvp.BEST.out.r, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [New Aptian]"]-mu["2 [Persistent Aptian]"]),
       col= "#35B779FF")
  plot(BAp.evn.BEST.out.r, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [Extinct Bathonian]"]-mu["2 [New Aptian]"]),
       col= "#35B779FF")
  
  plot(CB.evp.BEST.out.r, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [Extinct Carnian]"]-mu["2 [Persistent Carnian]"]),
       col= "#35B779FF")
  plot(CB.pvp.BEST.out.r, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [Persistent Carnian]"]-mu["2 [Persistent Bathonian]"]),
       col= "gray")
  plot(CB.nvp.BEST.out.r, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [New Bathonian]"]-mu["2 [Persistent Bathonian]"]),
       col= "#35B779FF")
  plot(CB.evn.BEST.out.r, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [Extinct Carnian]"]-mu["2 [New Bathonian]"]),
       col= "gray")
  
  plot(AnC.evp.BEST.out.r, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [Extinct Anisian]"]-mu["2 [Persistent Anisian]"]),
       col= "#35B779FF")
  plot(AnC.pvp.BEST.out.r, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [Persistent Anisian]"]-mu["2 [Persistent Carnian]"]),
       col= "gray")
  plot(AnC.nvp.BEST.out.r, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [New Carnian]"]-mu["2 [Persistent Carnian]"]),
       col= "#35B779FF")
  plot(AnC.evn.BEST.out.r, 
       showCurve = TRUE,
       credMass = 0.95,
       compVal = 0,
       comparisonColor = "black",
       xlab = expression(mu["1 [Extinct Anisian]"]-mu["2 [New Carnian]"]),
       col= "gray")
  
}




#{par(mfrow = c(3, 4))
  
  # Anisian:Carnian
  plot(AnC.evp.BEST.out, 
       main = "Anisian:Carnian \nDifference of Means", 
       sub = "Extinct [Anisian] - Persistent [Anisian]", 
       col="pink")
  plot(AnC.pvp.BEST.out, 
       main = "", 
       sub = "Persistent [Anisian] - Persistent [Carnian]")
  plot(AnC.nvp.BEST.out, 
       main = "", 
       sub = "New [Carnian] - Persistent [Carnian]",
       col = "pink")
  plot(AnC.evn.BEST.out, 
       main = "", 
       sub = "Extinct [Anisian] - New [Carnian]")
  
  #Carnian:Bathonian
  plot(CB.evp.BEST.out,
       main = "Carnian:Bathonian \nDifference of Means",
       sub = "Extinct [Carnian] - Persistent [Carnian]",
       col="pink")
  plot(CB.pvp.BEST.out,
       main = "",
       sub = "Persistent [Carnian] - Persistent [Bathonian]")
  plot(CB.nvp.BEST.out,
       main = "",
       sub = "New [Bathonian] - Persistent [Bathonian]",
       col="pink")
  plot(CB.evn.BEST.out,
       main = "",
       sub = "Extinct [Carnian] - New [Bathonian]")
  
  #Bathonian:Aptian
  plot(BAp.evp.BEST.out, 
       main = "Bathonian:Aptian \nDifference of Means", 
       sub = "Extinct [Bathonian] - Persistent [Bathonian]",
       col="pink")
  plot(BAp.pvp.BEST.out, 
       main = "", 
       sub = "Persistent [Bathonian] - Persistent [Aptian]")
  plot(BAp.nvp.BEST.out, 
       main = "", 
       sub = "New [Aptian] - Persistent [Aptian]",
       col="pink")
  plot(BAp.evn.BEST.out, 
       main = "", 
       sub = "Extinct [Bathonian] - New [Aptian]",
       col = "pink")
#}





















