### ------------------------------------------------------------------------ ###
### 3-area OM and using FLR's mse framework ####
### ------------------------------------------------------------------------ ###
### based on Ernesto's code in https://github.com/ejardim/deepsea_om/blob/main/tests/deepsea_like_om.R
### adapted to work with FLR's mse framework

### structure:
### the 3 different areas (a, b, c) are stored in the unit dimension
### the OM stock object has 3 units:
### - "a"/"b"/"c" for the different areas
### (removed) - "combined" for the total stocks (a+b+c)
### the areas are stored in the unit dimension because the mse functions
### mp() and goFish() can't handle areas

### the control object with the fishing targets uses the standard fwdControl()
### but is adapted to include the unit dimension (manual fix)

### NOTE: for this to work, a different branch of the mse package is required
### install with:
### remotes::install_github("shfischer/mse", ref = "dev_units")
### The only change is that the tracking of the fwdControl() object ignores
### the unit dimension. The standard mse version will fail.

### ------------------------------------------------------------------------ ###
### load packages and function ####
### ------------------------------------------------------------------------ ###

library(FLasher)
library(ggplot2)
library(tidyr)
library(dplyr)
library(mse)
library(FLife)
library(FLBRP)

source("utilities.R")

### ------------------------------------------------------------------------ ###
### create life-history OM first ####
### ------------------------------------------------------------------------ ###
### 1st: lets generate a life-history informed OM based on blackbellied angler
### code taken from ICES 2025 MSE training course material:
### https://github.com/shfischer/ICES_MSE_course_2025_public/blob/main/FLR_examples/101_generic_OMs/101_generic_OMs.md

### Life-history data from blackbellied angler fish
lh <- data.frame(a = 6.56E-7, b = 3.15, linf = 158.6, k = 0.119, t0 = -2.282, l50 = 102.8)

# check
len_at_0 <- lh$linf*(1-exp(-lh$k*(-lh$t0)))

### Max age: age at L = 0.95 * linf
tmax <- ceiling(log(0.05)/(-(lh$k)) + (lh$t0)) ### max age

### recruitment
### define unfished biomass and recruitment steepness
lh$v <- 1000 ### unfished SSB (arbitrary)
lh$s <- 0.7 ### recruitment steepness (arbitrary)
lh$sel1 <- lh$t0 - 1/lh$k*log(1-108/lh$linf)
lh$sel2 <- 3
lh$m <- 0.3

mac <- lh$t0 - 1/lh$k*log(1-60/lh$linf) # min age caught
pgr <- ceiling(lh$t0 - 1/lh$k*log(1-130/lh$linf))

### generate some more missing life-history parameters
lh_pars <- lhPar(lh[c("linf", "k", "t0", "a", "b", "l50", "s", "v", "sel1", "sel2")])

# need to think about selectivity parameters
# create FLBRP  with reference points
brp <- lhEql(lh_pars, range = c(min = 1, max = pgr, minfbar = 1, maxfbar = pgr-1, plusgroup = pgr), m = lh$m)
brp@mat <- brp@mat*0.7 # skip spawning
brp <- brp(brp)

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

### ------------------------------------------------------------------------ ###
### 3-area scenario specification ####
### ------------------------------------------------------------------------ ###
args_om <- list(
  bab = 1, # from a to b
  bba = 0, # from b to a
  bac = 0,    # from a to c
  bbc = 0,    # from b to c
  smpra = 0.8, # seamount a partial R
  smprb = 0.2, # seamount b partial R
  ftrga = 0.02,  # seamount a f target
  ftrgb = 0.09  # seamount b f target
)

### combined stock with units (saved in units)
stk0
stks <- FLCore:::expand(stk0, unit = c("a", "b", "c"))
dim(stks)
dimnames(stks)$unit
#knife edge in a
harvest(stks)[ac(5:10), ,"a"][]<-0

### recruitment model with units
srs <- FLCore:::expand(sr, unit = c("a", "b", "c"))
dim(srs)

# movement
flq0 <- stock.n(stk0)
flq0[] <- range(stk0)["min"]:range(stk0)["max"]
# from a to b
bab <- FLife:::logisticFn(age = flq0, params = FLPar(a50 = 3, asym = 1, ato95 = 1))*args_om$bab
# from b to a
bba <- bac <- bbc <- flq0
bba[] <- args_om$bba
bac[] <- args_om$bac
bbc[] <- args_om$bbc

mov <- FLQuants(bab = bab, bba = bba, bac = bac, bbc = bbc)
saveRDS(mov, file = "tests/mov.rds")

### ------------------------------------------------------------------------ ###
### include uncertainty (iterations) - recruitment SD ####
### ------------------------------------------------------------------------ ###

n_iter <- 10
sigmaR <- 0.6

### recruitment residuals
set.seed(123)
rec_res <- FLQuant(NA, dimnames = list(year = 1:100, unit = c("a", "b", "c")))
rec_res <- rlnoise(n = n_iter, len = rec_res %=% 0 - (sigmaR^2)/2,
                   sd = sigmaR, b = 0)
plot(rec_res)

### add iter dimension to recruitment model
srs10 <- propagate(srs, n_iter)
residuals(srs10) <- rec_res

### add iter dimensation to stock
stks10 <- propagate(stks, n_iter)

### FLom
### projection - FLom
om10 <- FLom(stock = stks10,
             sr = srs10,
             projection = mseCtrl(
               method = proj_units,
               args = list(mov = mov,
                           rec_split = list("a" = args_om$smpra,
                                            "b" = args_om$smprb,
                                            "c" = (1 - args_om$smpra - args_om$smprb)))))

### ------------------------------------------------------------------------ ###
### scenario: Ftrg
###
### ------------------------------------------------------------------------ ###

### HCR - mpCtrl object - constant F
mp_ctrl <- mpCtrl(list(
  est = mseCtrl(method = perfect.sa),
  hcr = mseCtrl(method = hcr_units_constF,
                args = list(
                  ftrg = c(a = args_om$ftrga,
                           b = args_om$ftrgb,
                           c = 0)))
))

args <- list(fy = 101, y0 = 1, iy = 2)

### run MSE
res10_Ftrg <- mp(om = om10, control = mp_ctrl, args = args, parallel = FALSE)

plot(FLStocks(a = stock(res10_Ftrg)[, ac(1:100), "a"],
              b = stock(res10_Ftrg)[, ac(1:100), "b"],
              c = stock(res10_Ftrg)[, ac(1:100), "c"]))

