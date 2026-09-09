#====================================================================
# deepsea like OM based on Lophius from WKLIFE
# 20260423EJ
#====================================================================

#install.packages(c("tidyr", "dplyr", "ggplot2"))
#install.packages(pkgs = c("FLCore", "FLasher", "FLife"),
#  repos = c(FLR="https://flr.r-universe.dev",
#  CRAN = "https://cloud.r-project.org"))

# load packages
library(FLife) 
library(FLasher)
library(ggplot2)
library(tidyr)
library(dplyr)

#====================================================================
# Generate a stock based on life history parameters (FLife))
#====================================================================

# stocks and life-history parameters
stocks_lh <- read.csv("stocks.csv")

# select lophius
lh <- stocks_lh[1, ]
agevec <- seq(0, 90, 2)

# von Bertalanffy parameters
lh[c("linf", "k", "t0")]
tmax <- lh$t0 - log(0.05)/lh$k

# plot growth
df_length <- data.frame(age = agevec)
df_length$length <- vonB(df_length$age, params = FLPar(lh))
df_length %>%
  ggplot(aes(x = age, y = length)) +
  geom_line() +
  geom_point(data = . %>% filter(age %in% agevec)) +
  geom_hline(yintercept = lh$linf, colour = "red", linetype = "dashed") +
  geom_vline(xintercept = tmax, colour = "red") +
  labs(x = "Age (years)", y = "Length (cm)") +
  theme_bw()

# length-weight parameters
lh[c("a", "b")]

# plot length-weight relationship
data.frame(length = seq(0, lh$linf),
           weight = c(lh$a * seq(0, lh$linf)^lh$b)) %>%
  ggplot(aes(x = length, y = weight)) +
  geom_line() +
  labs(x = "Length (cm)", y = "Weight (g)") +
  theme_bw()

# plot weight at age
df_length$weight <- lh$a* df_length$length^lh$b
df_length %>%
  ggplot(aes(x = age, y = weight)) +
  geom_line() +
  geom_point(data = . %>% filter(age %in% agevec)) +
  geom_vline(xintercept = tmax, colour = "red") +
  labs(x = "Age (years)", y = "Weight (g)") +
  theme_bw()

# plot M
df_length$M <- rep(1.5*lh["k"], nrow(df_length))

df_length %>%
  ggplot(aes(x = age, y = M)) +
  geom_line() +
  geom_point(data = . %>% filter(age %in% agevec)) +
  geom_vline(xintercept = tmax, colour = "red") +
  labs(x = "Age (years)", y = "Natural mortality") +
  theme_bw() +
  coord_cartesian(ylim = c(0, 3))

# calculate a50
a50 <- vonB(len = FLQuant(lh$l50), params = FLPar(lh))

# create template for ages
ages_tmp <- FLQuant(agevec, dimnames = list(age = agevec))
mat_tmp <- FLife:::logisticFn(age = ages_tmp, params = FLPar(a50 = a50, asym = 1, ato95 = 1))
as.data.frame(mat_tmp) %>%
  ggplot(aes(x = age, y = data)) +
  geom_line() +
  geom_point(data = . %>% filter(age %in% 0:90)) +
  geom_vline(xintercept = a50, colour = "red", linetype = "1111") +
  labs(x = "Age (years)", y = "Maturity (proportion)") +
  theme_bw()

# create template for ages
ages_tmp <- FLQuant(agevec, dimnames = list(age = agevec))
sel_tmp <- FLife:::dnormalFn(age = ages_tmp,
                             params = FLPar(sel1 = a50 + 1, sel2 = 1, sel3 = 5000))
as.data.frame(sel_tmp) %>%
  ggplot(aes(x = age, y = data)) +
  geom_line() +
  geom_point(data = . %>% filter(age %in% agevec)) +
  geom_vline(xintercept = a50, colour = "red", linetype = "1111") +
  labs(x = "Age (years)", y = "Selectivity") +
  theme_bw()

# define unfished biomass and recruitment steepness
lh$v <- 1000 ### unfished SSB (arbitrary)
lh$s <- 0.75 ### recruitment steepness (arbitrary)

# plot stock-recruitment relationship
# (relative to unfished recruitment, R0)
R0 <- 1 ### for plotting only
df_r <- data.frame(SSB = 0:1000)
df_r$R <- (0.8*R0*lh$s*df_r$SSB)/(0.2*lh$v*(1-lh$s)+(lh$s-0.2)*df_r$SSB)

df_r %>%
  ggplot(aes(x = SSB, y = R)) +
  geom_line() +
  geom_vline(xintercept = 0.2*lh$v, colour = "red", linetype = "dashed") +
  labs(x = "SSB", y = "Recruitment (R/R0)") +
  theme_bw()

# generate FLPar object with life-history parameters
lh_pars <- lhPar(lh[c("linf", "k", "t0", "a", "b", "l50", "s", "v")])

# create FLBRP  with reference points
brp <- lhEql(lh_pars, range = c(min = 1, max = ceiling(tmax),
                                minfbar = 1, maxfbar = ceiling(tmax),
                                plusgroup = floor(tmax)))

# show reference points
refpts(brp)

# for plotting, remove Fmax
refpts(brp)["fmax", ] <- NA
# plot
plot(brp) + theme_bw()

# get stock
stk0 <- as(brp, "FLStock")
plot(stk0) + theme_bw()

# get stock-recruitment model object
sr <- FLSR(params = params(brp), model = model(brp))
model(sr)
params(sr)

# constant f at Fmsy
fmsy <- c(refpts(brp)["msy", "harvest"])
# apply fishing history to stock
ctrl <- fwdControl(year = 25:101, quant = "f", value = fmsy)
stk0 <- fwd(stk0, control = ctrl, sr = sr)
plot(stk0) + labs(x = "Year") + theme_bw()

#====================================================================
# deepsea like OM
# 3 areas:
#  2 submarine mounts subject to fishing (a and b);
#  one area not subject to fishing (c),
#  common S/R
#====================================================================

stka_s01 <- stkb_s01 <- stkc_s01 <- stk_s01 <- stk0
sra.gm <- srb.gm <- src.gm <- as.FLSR(stk_s01, model="geomean")
harvest(stka_s01) <- harvest(stkb_s01) <- harvest(stk_s01)*0.5
harvest(stkc_s01)[] <- 0

#--------------------------------------------------------------------
# minimum movement (by age and time)
#--------------------------------------------------------------------
flq0 <- stock.n(stk_s01)
flq0[] <- range(stk_s01)["min"]:range(stk_s01)["max"]
# from a to b
bab <- flq0
bab[] <- 0.05
# from b to a
bba <- flq0
bba[] <- 0.05
# from a to c
bac <- 0
bac <- FLife:::logisticFn(age = flq0, params = FLPar(a50 = a50, asym = 1, ato95 = 1))*bac
# from b to c
bbc <- 0
bbc <- FLife:::logisticFn(age = flq0, params = FLPar(a50 = a50, asym = 1, ato95 = 1))*bbc
mov <- FLQuants(bab = bab, bba = bba, bac = bac, bbc = bbc)

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
  params(sra.gm)[] <- rfut*0.33
  params(srb.gm)[] <- rfut*0.33
  params(src.gm)[] <- rfut*0.33

  # set projection at F level of 0.05 in the seamounts, area c has
  # harvest = 0, so no F will happen
  ctrl <- fwdControl(year = c(i + 1), quant = "f", value = 0.05)
  stka_s01 <- fwd(stka_s01, control = ctrl, sr = sra.gm)
  stkb_s01 <- fwd(stkb_s01, control = ctrl, sr = srb.gm)
  stkc_s01 <- fwd(stkc_s01, control = ctrl, sr = src.gm)

  # update stock
  stk_s01 <- stka_s01 + stkb_s01 + stkc_s01

}

plot(window(FLStocks(a=stka_s01, b=stkb_s01, c=stkc_s01), 51))
plot(window(stk_s01, 51))

#--------------------------------------------------------------------
# strong vertical migration (by age and time)
#--------------------------------------------------------------------
stka_s02 <- stkb_s02 <- stkc_s02 <- stk_s02 <- stk0
sra.gm <- srb.gm <- src.gm <- as.FLSR(stka_s02, model="geomean")
harvest(stka_s02) <- harvest(stkb_s02) <- harvest(stk_s02)*0.5
harvest(stkc_s02)[] <- 0

flq0 <- stock.n(stk_s02)
flq0[] <- range(stk_s02)["min"]:range(stk_s02)["max"]
# from a to b
bab <- flq0
bab[] <- 0.05
# from b to a
bba <- flq0
bba[] <- 0.05
# from a to c
bac <- 1
bac <- FLife:::logisticFn(age = flq0, params = FLPar(a50 = a50, asym = 1, ato95 = 1))*bac
# from b to c
bbc <- 1
bbc <- FLife:::logisticFn(age = flq0, params = FLPar(a50 = a50, asym = 1, ato95 = 1))*bbc
mov <- FLQuants(bab = bab, bba = bba, bac = bac, bbc = bbc)

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
  params(sra.gm)[] <- rfut*0.33
  params(srb.gm)[] <- rfut*0.33
  params(src.gm)[] <- rfut*0.33

  # set projection at F level of 0.05 in the seamounts, area c has
  # harvest = 0, so no F will happen
  ctrl <- fwdControl(year = c(i + 1), quant = "f", value = 0.05)
  stka_s02 <- fwd(stka_s02, control = ctrl, sr = sra.gm)
  stkb_s02 <- fwd(stkb_s02, control = ctrl, sr = srb.gm)
  stkc_s02 <- fwd(stkc_s02, control = ctrl, sr = src.gm)

  # update stock
  stk_s02 <- stka_s02 + stkb_s02 + stkc_s02
}

plot(window(FLStocks(a=stka_s02, b=stkb_s02, c=stkc_s02), 51))
plot(window(stk_s02, 51))

#--------------------------------------------------------------------
# strong lateral migration (by age and time)
#--------------------------------------------------------------------
stka_s03 <- stkb_s03 <- stkc_s03 <- stk_s03 <- stk0
sra.gm <- srb.gm <- src.gm <- as.FLSR(stka_s03, model="geomean")
harvest(stka_s03) <- harvest(stkb_s03) <- harvest(stk_s03)*0.5
harvest(stkc_s03)[] <- 0

flq0 <- stock.n(stk_s03)
flq0[] <- range(stk_s03)["min"]:range(stk_s03)["max"]
# from a to b
bab <- flq0
bab[] <- 0.5
# from b to a
bba <- flq0
bba[] <- 0.05
# from a to c
bac <- 0
bac <- FLife:::logisticFn(age = flq0, params = FLPar(a50 = a50, asym = 1, ato95 = 1))*bac
# from b to c
bbc <- 0
bbc <- FLife:::logisticFn(age = flq0, params = FLPar(a50 = a50, asym = 1, ato95 = 1))*bbc
mov <- FLQuants(bab = bab, bba = bba, bac = bac, bbc = bbc)

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
  params(sra.gm)[] <- rfut*0.33
  params(srb.gm)[] <- rfut*0.33
  params(src.gm)[] <- rfut*0.33

  # set projection at F level of 0.05 in the seamounts, area c has
  # harvest = 0, so no F will happen
  ctrl <- fwdControl(year = c(i + 1), quant = "f", value = 0.05)
  stka_s03 <- fwd(stka_s03, control = ctrl, sr = sra.gm)
  stkb_s03 <- fwd(stkb_s03, control = ctrl, sr = srb.gm)
  stkc_s03 <- fwd(stkc_s03, control = ctrl, sr = src.gm)

  # update stock
  stk_s03 <- stka_s03 + stkb_s03 + stkc_s03
}

plot(window(FLStocks(a=stka_s03, b=stkb_s03, c=stkc_s03), 51))
plot(window(stk_s03, 51))


