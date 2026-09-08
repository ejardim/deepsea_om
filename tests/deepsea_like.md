---
title: "Deep-sea like OM based on Lophius (WKLIFE)"
author: "Ernesto Jardim"
date: "2026-09-08"
output:
  html_document:
    toc: true
    toc_float: true
    number_sections: true
    code_folding: show
---



# Introduction

This report builds a deep-sea-like operating model (OM) based on the
life-history parameters of *Lophius* (monkfish), as used in WKLIFE. The OM
represents three areas: two submarine mounts subject to fishing (`a` and
`b`) and one unfished refuge area (`c`), all sharing a common
stock-recruitment relationship.


``` r
install.packages(c("tidyr", "dplyr", "ggplot2"))
install.packages(pkgs = c("FLCore", "FLasher", "FLife"),
  repos = c(FLR = "https://flr.r-universe.dev",
            CRAN = "https://cloud.r-project.org"))
```


``` r
library(FLife)
library(FLasher)
library(ggplot2)
library(tidyr)
library(dplyr)
```

# Generate a stock from life-history parameters (FLife)

## Life-history parameters


``` r
# stocks and life-history parameters
stocks_lh <- read.csv("stocks.csv")

# select lophius
lh <- stocks_lh[1, ]
agevec <- seq(0, 90, 2)
```

## Growth


``` r
# von Bertalanffy parameters
lh[c("linf", "k", "t0")]
```

```
##    linf    k   t0
## 1 110.1 0.08 0.39
```

``` r
tmax <- lh$t0 - log(0.05) / lh$k
```


``` r
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
```

![plot of chunk growth-plot](figure/growth-plot-1.png)

## Length-weight relationship


``` r
lh[c("a", "b")]
```

```
##        a     b
## 1 0.0259 2.858
```


``` r
data.frame(length = seq(0, lh$linf),
           weight = c(lh$a * seq(0, lh$linf)^lh$b)) %>%
  ggplot(aes(x = length, y = weight)) +
  geom_line() +
  labs(x = "Length (cm)", y = "Weight (g)") +
  theme_bw()
```

![plot of chunk length-weight-plot](figure/length-weight-plot-1.png)

## Weight at age


``` r
df_length$weight <- lh$a * df_length$length^lh$b

df_length %>%
  ggplot(aes(x = age, y = weight)) +
  geom_line() +
  geom_point(data = . %>% filter(age %in% agevec)) +
  geom_vline(xintercept = tmax, colour = "red") +
  labs(x = "Age (years)", y = "Weight (g)") +
  theme_bw()
```

![plot of chunk weight-at-age](figure/weight-at-age-1.png)

## Natural mortality


``` r
df_length$M <- rep(1.5 * lh["k"], nrow(df_length))

df_length %>%
  ggplot(aes(x = age, y = M)) +
  geom_line() +
  geom_point(data = . %>% filter(age %in% agevec)) +
  geom_vline(xintercept = tmax, colour = "red") +
  labs(x = "Age (years)", y = "Natural mortality") +
  theme_bw() +
  coord_cartesian(ylim = c(0, 3))
```

```
## Error in `is.finite()`:
## ! default method not implemented for type 'list'
```

## Maturity


``` r
# calculate a50
a50 <- vonB(len = FLQuant(lh$l50), params = FLPar(lh))

# create template for ages
ages_tmp <- FLQuant(agevec, dimnames = list(age = agevec))
mat_tmp <- FLife:::logisticFn(age = ages_tmp,
                               params = FLPar(a50 = a50, asym = 1, ato95 = 1))

as.data.frame(mat_tmp) %>%
  ggplot(aes(x = age, y = data)) +
  geom_line() +
  geom_point(data = . %>% filter(age %in% 0:90)) +
  geom_vline(xintercept = a50, colour = "red", linetype = "1111") +
  labs(x = "Age (years)", y = "Maturity (proportion)") +
  theme_bw()
```

![plot of chunk maturity](figure/maturity-1.png)

## Selectivity


``` r
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
```

![plot of chunk selectivity](figure/selectivity-1.png)

## Stock-recruitment relationship


``` r
# define unfished biomass and recruitment steepness
lh$v <- 1000 ### unfished SSB (arbitrary)
lh$s <- 0.75 ### recruitment steepness (arbitrary)

# plot stock-recruitment relationship
# (relative to unfished recruitment, R0)
R0 <- 1 ### for plotting only
df_r <- data.frame(SSB = 0:1000)
df_r$R <- (0.8 * R0 * lh$s * df_r$SSB) /
  (0.2 * lh$v * (1 - lh$s) + (lh$s - 0.2) * df_r$SSB)

df_r %>%
  ggplot(aes(x = SSB, y = R)) +
  geom_line() +
  geom_vline(xintercept = 0.2 * lh$v, colour = "red", linetype = "dashed") +
  labs(x = "SSB", y = "Recruitment (R/R0)") +
  theme_bw()
```

![plot of chunk stock-recruitment](figure/stock-recruitment-1.png)

## Building the FLBRP and FLStock objects


``` r
# generate FLPar object with life-history parameters
lh_pars <- lhPar(lh[c("linf", "k", "t0", "a", "b", "l50", "s", "v")])

# create FLBRP with reference points
brp <- lhEql(lh_pars, range = c(min = 1, max = ceiling(tmax),
                                 minfbar = 1, maxfbar = ceiling(tmax),
                                 plusgroup = floor(tmax)))

# show reference points
refpts(brp)
```

```
## An object of class "FLPar"
##         quant
## refpt    harvest   yield     rec       ssb       biomass   revenue   cost     
##   virgin  0.00e+00  0.00e+00  5.22e+03  1.00e+03  1.66e+04        NA        NA
##   msy     6.18e-02  2.02e+01  4.20e+03  2.56e+02  1.28e+04        NA        NA
##   crash   2.80e-01  1.72e-06  2.97e-04  4.74e-06  8.94e-04        NA        NA
##   f0.1    5.59e-02  2.01e+01  4.31e+03  2.83e+02  1.32e+04        NA        NA
##   fmax    5.00e-01 -2.63e+01 -4.52e+03 -4.03e+01 -1.36e+04        NA        NA
##         quant
## refpt    profit   
##   virgin        NA
##   msy           NA
##   crash         NA
##   f0.1          NA
##   fmax          NA
## units:  NA
```


``` r
# for plotting, remove Fmax
refpts(brp)["fmax", ] <- NA
plot(brp) + theme_bw()
```

![plot of chunk flbrp-plot](figure/flbrp-plot-1.png)


``` r
# get stock
stk <- as(brp, "FLStock")
plot(stk) + theme_bw()
```

![plot of chunk flstock](figure/flstock-1.png)

``` r
# get stock-recruitment model object
sr <- FLSR(params = params(brp), model = model(brp))
model(sr)
```

```
## rec ~ a * ssb/(b + ssb)
## <environment: 0x1c39d380>
```

``` r
params(sr)
```

```
## An object of class "FLPar"
## params
##      a      b 
## 5690.7   90.9 
## units:  NA
```

## Fishing history at Fmsy


``` r
# constant f at Fmsy
fmsy <- c(refpts(brp)["msy", "harvest"])

# apply fishing history to stock
ctrl <- fwdControl(year = 25:101, quant = "f", value = fmsy)
stk <- fwd(stk, control = ctrl, sr = sr)
plot(stk) + labs(x = "Year") + theme_bw()
```

![plot of chunk fmsy-history](figure/fmsy-history-1.png)

# Deep-sea-like operating model

The OM represents three areas: two submarine mounts subject to fishing
(`a` and `b`) and one refuge area not subject to fishing (`c`), sharing a
common stock-recruitment relationship.


``` r
stka <- stkb <- stkc <- stk
sra.gm <- srb.gm <- src.gm <- as.FLSR(stka, model = "geomean")
harvest(stka) <- harvest(stkb) <- harvest(stk) * 0.5
harvest(stkc)[] <- 0

yinit <- 50
```

## Scenario 1: minimum movement (by age and time)


``` r
flq0 <- stock.n(stk)
flq0[] <- range(stk)["min"]:range(stk)["max"]

# from a to b
bab <- flq0
bab[] <- 0.05
# from b to a
bba <- flq0
bba[] <- 0.05
# from a to c
bac <- 0.05
bac <- FLife:::logisticFn(age = flq0, params = FLPar(a50 = a50, asym = 1, ato95 = 1)) * bac
# from b to c
bbc <- 0.05
bbc <- FLife:::logisticFn(age = flq0, params = FLPar(a50 = a50, asym = 1, ato95 = 1)) * bbc

mov <- FLQuants(bab = bab, bba = bba, bac = bac, bbc = bbc)
```


``` r
for (i in 1:50) {
  # common recruitment
  rfut <- predict(sr, ssb = ssb(stk)[, ac(yinit + i - 1)])

  # movement
  # from a to b
  stknab <- stock.n(stka)[, ac(yinit + i)] * mov[["bab"]][, ac(yinit + i)]
  # from b to a
  stknba <- stock.n(stkb)[, ac(yinit + i)] * mov[["bba"]][, ac(yinit + i)]
  # from a to c
  stknac <- stock.n(stka)[, ac(yinit + i)] * mov[["bac"]][, ac(yinit + i)]
  # from b to c
  stknbc <- stock.n(stkb)[, ac(yinit + i)] * mov[["bbc"]][, ac(yinit + i)]

  stock.n(stka)[, ac(yinit + i)] <- stock.n(stka)[, ac(yinit + i)] - stknab - stknac + stknba
  stock.n(stkb)[, ac(yinit + i)] <- stock.n(stkb)[, ac(yinit + i)] - stknba - stknbc + stknab
  stock.n(stkc)[, ac(yinit + i)] <- stock.n(stkc)[, ac(yinit + i)] + stknac + stknbc

  # distribute recruitment evenly
  params(sra.gm)[] <- rfut * 0.33
  params(srb.gm)[] <- rfut * 0.33
  params(src.gm)[] <- rfut * 0.33

  # set projection at F level of 0.05 in the seamounts, area c has
  # harvest = 0, so no F will happen
  ctrl <- fwdControl(year = c(yinit + i, yinit + i + 1), quant = "f", value = 0.05)
  stka <- fwd(stka, control = ctrl, sr = sra.gm)
  stkb <- fwd(stkb, control = ctrl, sr = srb.gm)
  stkc <- fwd(stkc, control = ctrl, sr = src.gm)

  # update stock
  stk <- stka + stkb + stkc
}
```


``` r
plot(window(FLStocks(a = stka, b = stkb, c = stkc), 51))
```

![plot of chunk plot-minimum](figure/plot-minimum-1.png)

``` r
plot(window(stk, 51))
```

![plot of chunk plot-minimum](figure/plot-minimum-2.png)

## Scenario 2: strong vertical migration (by age and time)


``` r
flq0 <- stock.n(stk)
flq0[] <- range(stk)["min"]:range(stk)["max"]

# from a to b
bab <- flq0
bab[] <- 0.05
# from b to a
bba <- flq0
bba[] <- 0.05
# from a to c
bac <- 0.5
bac <- FLife:::logisticFn(age = flq0, params = FLPar(a50 = a50, asym = 1, ato95 = 1)) * bac
# from b to c
bbc <- 0.5
bbc <- FLife:::logisticFn(age = flq0, params = FLPar(a50 = a50, asym = 1, ato95 = 1)) * bbc

mov <- FLQuants(bab = bab, bba = bba, bac = bac, bbc = bbc)
```


``` r
for (i in 1:50) {
  # common recruitment
  rfut <- predict(sr, ssb = ssb(stk)[, ac(yinit + i - 1)])

  # movement
  # from a to b
  stknab <- stock.n(stka)[, ac(yinit + i)] * mov[["bab"]][, ac(yinit + i)]
  # from b to a
  stknba <- stock.n(stkb)[, ac(yinit + i)] * mov[["bba"]][, ac(yinit + i)]
  # from a to c
  stknac <- stock.n(stka)[, ac(yinit + i)] * mov[["bac"]][, ac(yinit + i)]
  # from b to c
  stknbc <- stock.n(stkb)[, ac(yinit + i)] * mov[["bbc"]][, ac(yinit + i)]

  stock.n(stka)[, ac(yinit + i)] <- stock.n(stka)[, ac(yinit + i)] - stknab - stknac + stknba
  stock.n(stkb)[, ac(yinit + i)] <- stock.n(stkb)[, ac(yinit + i)] - stknba - stknbc + stknab
  stock.n(stkc)[, ac(yinit + i)] <- stock.n(stkc)[, ac(yinit + i)] + stknac + stknbc

  # distribute recruitment evenly
  params(sra.gm)[] <- rfut * 0.33
  params(srb.gm)[] <- rfut * 0.33
  params(src.gm)[] <- rfut * 0.33

  # set projection at F level of 0.05 in the seamounts, area c has
  # harvest = 0, so no F will happen
  ctrl <- fwdControl(year = c(yinit + i, yinit + i + 1), quant = "f", value = 0.05)
  stka <- fwd(stka, control = ctrl, sr = sra.gm)
  stkb <- fwd(stkb, control = ctrl, sr = srb.gm)
  stkc <- fwd(stkc, control = ctrl, sr = src.gm)

  # update stock
  stk <- stka + stkb + stkc
}
```


``` r
plot(window(FLStocks(a = stka, b = stkb, c = stkc), 51))
```

![plot of chunk plot-strong](figure/plot-strong-1.png)

``` r
plot(window(stk, 51))
```

![plot of chunk plot-strong](figure/plot-strong-2.png)

# Session info


``` r
sessionInfo()
```

```
## R version 4.5.2 (2025-10-31)
## Platform: x86_64-pc-linux-gnu
## Running under: openSUSE Leap 16.0
## 
## Matrix products: default
## BLAS:   /usr/local/R452/lib64/R/lib/libRblas.so 
## LAPACK: /usr/local/R452/lib64/R/lib/libRlapack.so;  LAPACK version 3.12.1
## 
## locale:
##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
## 
## time zone: Europe/Lisbon
## tzcode source: system (glibc)
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] dplyr_1.2.0     tidyr_1.3.2     FLasher_0.7.3   FLFishery_0.4.0
##  [5] FLife_3.4.0     ggplotFL_2.7.2  patchwork_1.3.2 ggrepel_0.9.8  
##  [9] ggplot2_4.0.3   FLCore_2.6.33   lattice_0.22-7  knitr_1.51     
## 
## loaded via a namespace (and not attached):
##  [1] Matrix_1.7-4       gtable_0.3.6       compiler_4.5.2     tidyselect_1.2.1  
##  [5] Rcpp_1.1.2         gridExtra_2.3      scales_1.4.0       R6_2.6.1          
##  [9] labeling_0.4.3     generics_0.1.4     iterators_1.0.14   MASS_7.3-65       
## [13] tibble_3.3.1       pillar_1.11.1      RColorBrewer_1.1-3 rlang_1.2.0       
## [17] xfun_0.56          S7_0.2.2           otel_0.2.0         cli_3.6.6         
## [21] withr_3.0.3        magrittr_2.0.5     FLBRP_2.5.10       grid_4.5.2        
## [25] cowplot_1.2.0      lifecycle_1.0.5    vctrs_0.7.3        evaluate_1.0.5    
## [29] glue_1.8.1         data.table_1.18.4  farver_2.1.2       stats4_4.5.2      
## [33] purrr_1.2.1        tools_4.5.2        pkgconfig_2.0.3
```
