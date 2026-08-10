# test_fwd_flstock_units.R - DESC
# /home/mosqu003/Sandbox/test_fwd_flstock_units.R

# Copyright (c) WMR, 2026.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(FLasher)

data(ple4sex)

# SRR rec ~ 225000
srr <- predictModel(model=rec~a, params=FLPar(a=225000))

# RENAME unit to c('N', 'S')
x <- ple4sex
x <- qapply(x, function(s) {
  dimnames(s)$unit <- c("N", "S")
  s})

# SET  both units maturity
mat(x)[,,'S'] <- mat(x)[,,'N']

# Fbar target
target <- unitSums(fbar(x)[, ac(1991:2000)]) %=% 0.2

# --- RUN fwd with 2 sexes, unit=c('F', 'M'), and one SRR
te1 <- fwd(ple4sex, sr=srr, fbar=target)

# REC is split in 2, sratio = 0.5
rec(te1)[, ac(1990:2000)]


# --- RUN fwd with 2 units, units=c('N', 'S'), and  one SRR
te2 <- fwd(x, sr=srr, fbar=target)

# REC is repeated
rec(te2)[, ac(1990:2000)]


# --- RUN fwd with 2 units, separate SRRs

sr3 <- predictModel(model=rec~a, params=FLPar(c(150000, 300000),
  dimnames=list(params='a', unit=c("N", "S"), iter=1)))

te3 <- fwd(x, sr=sr3, fbar=target)

# REC is taken from each unit in SRR
rec(te3)[, ac(1990:2000)]


# ---  RUN fwd with 2 units, mat = 0 in one unit, and one SRR

y <- x
mat(y)[,,'N'] <- 0

te4 <- fwd(x, sr=srr, fbar=target)

# REC is repeated
rec(te4)[, ac(1990:2000)]


# ---  RUN fwd with 2 units, mat = 0 in one unit, and one SRR w/SSB input

sr5 <- predictModel(model=segreg()$model, params=FLPar(a=9, b=2200))

te5 <- fwd(x, sr=sr5, fbar=target)

# REC is repeated, SSB is total SSB
rec(te5)[, ac(1990:2000)]

# CHECK predictions
predict(sr5, ssb=unitSums(ssb(x)[, ac(1991:2000)]))

predict(sr5, ssb=quantSums(stock.n(x) * stock.wt(x))[, ac(1991:2000), 'S'])

predict(sr5, ssb=quantSums(stock.n(x) * stock.wt(x))[, ac(1991:2000), 'N'])


# --- RUN fwd with 2 units, separate SRRs and different selectivity

z <- x
harvest(z)[, ac(1991:2000), 'S'] <- 0

te6 <- fwd(z, sr=sr3, fbar=target)

# NO F in 'S'
fbar(te6)[, ac(1990:2000)]

rec(te6)[, ac(1990:2000)]

