### ------------------------------------------------------------------------ ###
### OM - projection module ####
### ------------------------------------------------------------------------ ###
### set up projection module
### projection - method
proj_units <- function(om, ctrl,
                       mov,
                       rec_split,
                       ...) {
                      ### get some dimensions
  # ay and etc needs to be coded in a more mse compliant way
  # this variables should be computed using information in the args
  # list, like:
  # mseargs <- list(iy=iy, fy=fy, data_lag=dl, management_lag=ml, frq=af)
  ay <- unique(ctrl@target$year)
  units <- dimnames(om)$unit
  units <- setdiff(units, "combined")
  names(units) <- units ### named vector needed for subsetting later

  ### get combined SSB
  ### default ssb() function sums up units
  ssb_ay <- ssb(stock(om)[, ac(ay - 1)])
  # ssb_ay <- lapply(units, function(x) {
  #
  # })
  # ssb_ay <- Reduce(x = ssb_ay, f = "+")
  # rfut <- predict(sr(om)[,, "combined"], ssb = ssb_ay)
  # common recruitment
  rfut <- predict(sr(om), ssb = ssb_ay)

  # movement
  # from a to b
  stknab <- stock.n(om)[, ac(ay-1), "a"]*mov[["bab"]][, ac(ay-1)]
  # from b to a
  stknba <- stock.n(om)[, ac(ay-1), "b"]*mov[["bba"]][, ac(ay-1)]
  # from a to c
  stknac <- stock.n(om)[, ac(ay-1), "a"]*mov[["bac"]][, ac(ay-1)]
  # from b to c
  stknbc <- stock.n(om)[, ac(ay-1), "b"]*mov[["bbc"]][, ac(ay-1)]

  stock.n(om)[, ac(ay-1), "a"] <- stock.n(om)[, ac(ay-1), "a"] - stknab - stknac + stknba
  stock.n(om)[, ac(ay-1), "b"] <- stock.n(om)[, ac(ay-1), "b"] - stknba - stknbc + stknab
  stock.n(om)[, ac(ay-1), "c"] <- stock.n(om)[, ac(ay-1), "c"] + stknac + stknbc

  ### recruitment model
  ### extract models for units a/b/c
  ### need to duplicate SR models because params slot (FLPar) cannot handle units...
  sr_list <- lapply(units, function(x) {
    sr_i <- sr(om)[,, ac(x)]
    ### adjust model type to geomean to define rec value manually
    model(sr_i) <- "geomean"
    ### bugs in FLSR after setting "geomean" -iterations deleted -> add again
    params(sr_i) <- propagate(params(sr_i), dim(sr_i)[6])
    residuals(sr_i) <- residuals(sr(om)[,, ac(x)])
    ### distribute recruitment between units
    params(sr_i)[] <- rfut * rec_split[[x]]
    return(sr_i)
  })

  ### ctrl (fwdControl) object contains all units
  ### we need to split them into units because FLasher can't handle units
  ctrl_list <- lapply(units, function(x) ctrl)
  ### subset each control object
  ctrl_list <- lapply(units, function(x) {
    ctrl_i <- ctrl
    ctrl_i@target <- ctrl_i@target[ctrl_i@target$unit == ac(x),, drop = FALSE]
    ctrl_i@iters <- ctrl_i@iters[ac(x),,, drop = FALSE]
    return(ctrl_i)
  })

  ### project stocks forward - by unit
  ### use for-loop to adjust values in original object to avoid more objects
  for (unit_i in units) {
    #debugonce(fwd, signature = c("FLStock", "missing", "fwdControl"))
    stock(om)[,, unit_i] <- fwd(stock(om)[,, unit_i],
                                control = ctrl_list[[unit_i]],
                                sr = sr_list[[unit_i]],
                                residuals = residuals(sr_list[[unit_i]]))
  }

  ### combine units
  ### this should sum catch and stock numbers etc but keep M etc untouched
  # stock(om)[, ac(ay), "combined"] <- Reduce(x = lapply(units, function(x) {
  #   stock(om)[, ac(ay), x]
  # }), f = "+")
  # catch(om)[, ac(ay), "combined"]
  # catch(om)[, ac(ay), "a"]
  # catch(om)[, ac(ay), "b"]
  # catch(om)[, ac(ay), "c"]

  return(list(om = om))

}

### ------------------------------------------------------------------------ ###
### HCR - constant F ####
### ------------------------------------------------------------------------ ###

### HCR - function
### constant F by unit
hcr_units_constF <- function(stk, ftrg, args, tracking,
                             ...) {
  ### get some dimensions
  units <- dimnames(stk)$unit
  units <- setdiff(units, "combined")
  names(units) <- units ### named vector needed for subsetting later

  ### create fwdControl object with units
  ctrl <- fwdControl(list(year = args$ay + args$management_lag,
                          quant = "fbar",
                          unit = units,
                          value = NA))
  n_units <- length(units)
  n_its <- dims(stk)$iter
  iters_array <- array(NA, dim = c(n_units, 3, n_its),
                       dimnames = list(row = units,
                                       value = c("min", "value", "max"),
                                       iter = seq(n_its)))
  iters_array[, "value", ] <- unlist(ftrg)
  ctrl@iters <- iters_array

  list(ctrl = ctrl, tracking = tracking)
}



