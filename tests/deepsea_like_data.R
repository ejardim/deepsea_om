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

save(list=ls(), file="lophius_simdata.RData")
