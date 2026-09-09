#====================================================================
# deepsea like OM based on Lophius from WKLIFE
# 20260423EJ
#====================================================================

#install.packages(c("tidyr", "dplyr", "ggplot2"))
#install.packages(pkgs = c("FLCore", "FLasher", "FLife"),
#  repos = c(FLR="https://flr.r-universe.dev",
#  CRAN = "https://cloud.r-project.org"))

# load packages
library(FLasher)
library(ggplot2)
load("")

#====================================================================
# deepsea like OM
# 3 areas:
#  2 submarine mounts subject to fishing (a and b);
#  one area not subject to fishing (c),
#  common S/R
#====================================================================

#--------------------------------------------------------------------
# scenario 01: minimum movement (by age and time)
#--------------------------------------------------------------------
args_om <- list(
  smpf=0.5,   # seamount partial f
  bab = 0.05, # from a to b
  bba = 0.05, # from b to a
  bac = 0,    # from a to c
  bbc = 0,    # from b to c
  smpra = 0.333, # seamount a partial R
  smprb = 0.333, # seamount b partial R
  ftrga = 0.05,  # seamount a f target
  ftrgb = 0.05  # seamount b f target
)

# objects
stka_s01 <- stkb_s01 <- stkc_s01 <- stk_s01 <- stk0
sra.gm <- srb.gm <- src.gm <- as.FLSR(stk_s01, model="geomean")
harvest(stka_s01) <- harvest(stk_s01)*args_om$smpf
harvest(stkb_s01) <- harvest(stk_s01)*(1-args_om$smpf)
harvest(stkc_s01)[] <- 0

# movement
flq0 <- stock.n(stk_s01)
flq0[] <- range(stk_s01)["min"]:range(stk_s01)["max"]
# from a to b
bab <- flq0
bab[] <- args_om$bab
# from b to a
bba <- flq0
bba[] <- args_om$bba
# from a to c
bac <- FLife:::logisticFn(age = flq0, params = FLPar(a50 = a50, asym = 1, ato95 = 1))*args_om$bac
# from b to c
bbc <- FLife:::logisticFn(age = flq0, params = FLPar(a50 = a50, asym = 1, ato95 = 1))*args_om$bbc
mov <- FLQuants(bab = bab, bba = bba, bac = bac, bbc = bbc)

# projection loop
for(i in 50:100){
  # common recruitment
  rfut <- predict(sr, ssb=ssb(stk_s01)[,ac(i - 1)])

  # movement
  # from a to b
  stknab <- stock.n(stka_s01)[,ac(i)]*mov[["bab"]][,ac(i)]
  # from b to a
  stknba <- stock.n(stkb_s01)[,ac(i)]*mov[["bba"]][,ac(i)]
  # from a to c
  stknac <- stock.n(stka_s01)[,ac(i)]*mov[["bac"]][,ac(i)]
  # from b to c
  stknbc <- stock.n(stkb_s01)[,ac(i)]*mov[["bbc"]][,ac(i)]

  stock.n(stka_s01)[,ac(i)] <- stock.n(stka_s01)[,ac(i)] - stknab - stknac + stknba
  stock.n(stkb_s01)[,ac(i)] <- stock.n(stkb_s01)[,ac(i)] - stknba - stknbc + stknab
  stock.n(stkc_s01)[,ac(i)] <- stock.n(stkc_s01)[,ac(i)] + stknac + stknbc

  # distribute recruitment evenly
  params(sra.gm)[] <- rfut*args_om$smpra
  params(srb.gm)[] <- rfut*args_om$smprb
  params(src.gm)[] <- rfut*(1-args_om$smpra-args_om$smprb)

  # set projection at F level of ftrg in the seamounts, area c has
  # harvest = 0, so no F will happen
  ctrl <- fwdControl(year = c(i + 1), quant = "f", value = args_om$ftrga)
  stka_s01 <- fwd(stka_s01, control = ctrl, sr = sra.gm)
  ctrl <- fwdControl(year = c(i + 1), quant = "f", value = args_om$ftrgb)
  stkb_s01 <- fwd(stkb_s01, control = ctrl, sr = srb.gm)
  ctrl <- fwdControl(year = c(i + 1), quant = "f", value = 0)
  stkc_s01 <- fwd(stkc_s01, control = ctrl, sr = src.gm)

  # update stock
  stk_s01 <- stka_s01 + stkb_s01 + stkc_s01

}

plot(window(FLStocks(a=stka_s01, b=stkb_s01, c=stkc_s01), 50, 100))
plot(window(stk_s01, 50, 100))

#--------------------------------------------------------------------
# scenario 02: strong vertical migration (by age and time)
#--------------------------------------------------------------------
args_om <- list(
  smpf=0.5,   # seamount partial f
  bab = 0.05, # from a to b
  bba = 0.05, # from b to a
  bac = 1,    # from a to c
  bbc = 1,    # from b to c
  smpra = 0.333, # seamount a partial R
  smprb = 0.333, # seamount b partial R
  ftrga = 0.05,  # seamount a f target
  ftrgb = 0.05  # seamount b f target
)

# objects
stka_s02 <- stkb_s02 <- stkc_s02 <- stk_s02 <- stk0
sra.gm <- srb.gm <- src.gm <- as.FLSR(stk_s02, model="geomean")
harvest(stka_s02) <- harvest(stk_s02)*args_om$smpf
harvest(stkb_s02) <- harvest(stk_s02)*(1-args_om$smpf)
harvest(stkc_s02)[] <- 0

# movement
flq0 <- stock.n(stk_s02)
flq0[] <- range(stk_s02)["min"]:range(stk_s02)["max"]
# from a to b
bab <- flq0
bab[] <- args_om$bab
# from b to a
bba <- flq0
bba[] <- args_om$bba
# from a to c
bac <- FLife:::logisticFn(age = flq0, params = FLPar(a50 = a50, asym = 1, ato95 = 1))*args_om$bac
# from b to c
bbc <- FLife:::logisticFn(age = flq0, params = FLPar(a50 = a50, asym = 1, ato95 = 1))*args_om$bbc
mov <- FLQuants(bab = bab, bba = bba, bac = bac, bbc = bbc)

# projection loop
for(i in 50:100){
  # common recruitment
  rfut <- predict(sr, ssb=ssb(stk_s02)[,ac(i - 1)])

  # movement
  # from a to b
  stknab <- stock.n(stka_s02)[,ac(i)]*mov[["bab"]][,ac(i)]
  # from b to a
  stknba <- stock.n(stkb_s02)[,ac(i)]*mov[["bba"]][,ac(i)]
  # from a to c
  stknac <- stock.n(stka_s02)[,ac(i)]*mov[["bac"]][,ac(i)]
  # from b to c
  stknbc <- stock.n(stkb_s02)[,ac(i)]*mov[["bbc"]][,ac(i)]

  stock.n(stka_s02)[,ac(i)] <- stock.n(stka_s02)[,ac(i)] - stknab - stknac + stknba
  stock.n(stkb_s02)[,ac(i)] <- stock.n(stkb_s02)[,ac(i)] - stknba - stknbc + stknab
  stock.n(stkc_s02)[,ac(i)] <- stock.n(stkc_s02)[,ac(i)] + stknac + stknbc

  # distribute recruitment evenly
  params(sra.gm)[] <- rfut*args_om$smpra
  params(srb.gm)[] <- rfut*args_om$smprb
  params(src.gm)[] <- rfut*(1-args_om$smpra-args_om$smprb)

  # set projection at F level of ftrg in the seamounts, area c has
  # harvest = 0, so no F will happen
  ctrl <- fwdControl(year = c(i + 1), quant = "f", value = args_om$ftrga)
  stka_s02 <- fwd(stka_s02, control = ctrl, sr = sra.gm)
  ctrl <- fwdControl(year = c(i + 1), quant = "f", value = args_om$ftrgb)
  stkb_s02 <- fwd(stkb_s02, control = ctrl, sr = srb.gm)
  ctrl <- fwdControl(year = c(i + 1), quant = "f", value = 0)
  stkc_s02 <- fwd(stkc_s02, control = ctrl, sr = src.gm)

  # update stock
  stk_s02 <- stka_s02 + stkb_s02 + stkc_s02

}

plot(window(FLStocks(a=stka_s02, b=stkb_s02, c=stkc_s02), 50, 100))
plot(window(stk_s02,50, 100))

#--------------------------------------------------------------------
# scenario 03: strong lateral migration (by age and time)
#--------------------------------------------------------------------
args_om <- list(
  smpf=0.5,   # seamount partial f
  bab = 0.75, # from a to b
  bba = 0.05, # from b to a
  bac = 0,    # from a to c
  bbc = 0,    # from b to c
  smpra = 0.333, # seamount a partial R
  smprb = 0.333, # seamount b partial R
  ftrga = 0.05,  # seamount a f target
  ftrgb = 0.05  # seamount b f target
)

# objects
stka_s03 <- stkb_s03 <- stkc_s03 <- stk_s03 <- stk0
sra.gm <- srb.gm <- src.gm <- as.FLSR(stk_s03, model="geomean")
harvest(stka_s03) <- harvest(stk_s03)*args_om$smpf
harvest(stkb_s03) <- harvest(stk_s03)*(1-args_om$smpf)
harvest(stkc_s03)[] <- 0

# movement
flq0 <- stock.n(stk_s03)
flq0[] <- range(stk_s03)["min"]:range(stk_s03)["max"]
# from a to b
bab <- flq0
bab[] <- args_om$bab
# from b to a
bba <- flq0
bba[] <- args_om$bba
# from a to c
bac <- FLife:::logisticFn(age = flq0, params = FLPar(a50 = a50, asym = 1, ato95 = 1))*args_om$bac
# from b to c
bbc <- FLife:::logisticFn(age = flq0, params = FLPar(a50 = a50, asym = 1, ato95 = 1))*args_om$bbc
mov <- FLQuants(bab = bab, bba = bba, bac = bac, bbc = bbc)

# projection loop
for(i in 50:100){
  # common recruitment
  rfut <- predict(sr, ssb=ssb(stk_s03)[,ac(i - 1)])

  # movement
  # from a to b
  stknab <- stock.n(stka_s03)[,ac(i)]*mov[["bab"]][,ac(i)]
  # from b to a
  stknba <- stock.n(stkb_s03)[,ac(i)]*mov[["bba"]][,ac(i)]
  # from a to c
  stknac <- stock.n(stka_s03)[,ac(i)]*mov[["bac"]][,ac(i)]
  # from b to c
  stknbc <- stock.n(stkb_s03)[,ac(i)]*mov[["bbc"]][,ac(i)]

  stock.n(stka_s03)[,ac(i)] <- stock.n(stka_s03)[,ac(i)] - stknab - stknac + stknba
  stock.n(stkb_s03)[,ac(i)] <- stock.n(stkb_s03)[,ac(i)] - stknba - stknbc + stknab
  stock.n(stkc_s03)[,ac(i)] <- stock.n(stkc_s03)[,ac(i)] + stknac + stknbc

  # distribute recruitment evenly
  params(sra.gm)[] <- rfut*args_om$smpra
  params(srb.gm)[] <- rfut*args_om$smprb
  params(src.gm)[] <- rfut*(1-args_om$smpra-args_om$smprb)

  # set projection at F level of ftrg in the seamounts, area c has
  # harvest = 0, so no F will happen
  ctrl <- fwdControl(year = c(i + 1), quant = "f", value = args_om$ftrga)
  stka_s03 <- fwd(stka_s03, control = ctrl, sr = sra.gm)
  ctrl <- fwdControl(year = c(i + 1), quant = "f", value = args_om$ftrgb)
  stkb_s03 <- fwd(stkb_s03, control = ctrl, sr = srb.gm)
  ctrl <- fwdControl(year = c(i + 1), quant = "f", value = 0)
  stkc_s03 <- fwd(stkc_s03, control = ctrl, sr = src.gm)

  # update stock
  stk_s03 <- stka_s03 + stkb_s03 + stkc_s03

}

plot(window(FLStocks(a=stka_s03, b=stkb_s03, c=stkc_s03), 50, 100))
plot(window(stk_s03, 50, 100))

