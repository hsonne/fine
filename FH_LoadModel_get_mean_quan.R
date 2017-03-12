library(kwb.utils)

### run FH_LoadModel first ###

# following lists are required:
# x_annual_loads_rain
# x_annual_loads_sew

## SOURCE ;)

# Define file types
types <- c(
  NEU_meanln_sdln = "annual mean concentrations of rainwater",
  Vol_rain= "rain runoff",
  Vol_sewage = "sewage runoff"
)

# load data
name <- "NEU_meanln_sdln"
x_conc_NEU <- readTableOrStop(data.dir, name, types[name])

name <- "Vol_rain"
vol_rain <- readTableOrStop(data.dir, name, types[name], dec = ",")

name <- "Vol_sewage"
vol_sewage <- readTableOrStop(data.dir, name, types[name])

# get Names for Variables, SUWs and paths

variables <- selectColumns(x_conc_NEU, "VariableName")
SUW_Names_rain <- unique(selectColumns(vol_rain, "SUW"))
SUW_Names_sew <- unique(selectColumns(vol_sewage, "SUW"))

paths_rain <- c("SEP", "CSO", "WWTP", "TOT")
paths_sew <- c("CSO", "WWTP", "TOT")

# MAIN
if (FALSE)
{
  # Provide loads from list "x_annual_loads_rain"
  load_rain_cso      <- x_annual_loads_rain$load_rain_cso
  load_rain_sep      <- x_annual_loads_rain$load_rain_sep
  load_rain_wwtp     <- x_annual_loads_rain$load_rain_wwtp
  load_rain_sum_path <- x_annual_loads_rain$load_rain_sum_paths
  
  # Provide loads from list "x_annual_loads_sew"
  load_sew_cso      <- x_annual_loads_sew$load_sew_cso
  load_sew_wwtp     <- x_annual_loads_sew$load_sew_wwtp
  load_sew_sum_path <- x_annual_loads_sew$load_sew_sum_paths
  
  # Define function arguments for meanQuantiles
  args1 <- list(
    
    rain_cso  = list(offset = 1, SUW_Names_rain, variables, load_rain_cso),
    rain_sep  = list(offset = 1, SUW_Names_rain, variables, load_rain_sep),
    rain_wwtp = list(offset = 1, SUW_Names_rain, variables, load_rain_wwtp),
    rain_sum  = list(offset = 0, SUW_Names_rain, variables, load_rain_sum_path),
    
    sew_cso  = list(offset = 1, SUW_Names_sew, variables, load_sew_cso),
    sew_wwtp = list(offset = 1, SUW_Names_sew, variables, load_sew_wwtp),
    sew_sum  = list(offset = 0, SUW_Names_sew, variables, load_sew_sum_path)
  )
  
  ## get mean and quantiles for loads in rainwater and all pathways
  load_rain_cso_mean_quan  <- do.call(meanQuantiles, args1$rain_cso)
  load_rain_sep_mean_quan  <- do.call(meanQuantiles, args1$rain_sep)
  load_rain_wwtp_mean_quan <- do.call(meanQuantiles, args1$rain_wwtp)

  ## get mean and quantiles for loads in sewage and all pathways
  load_sew_cso_mean_quan  <- do.call(meanQuantiles, args1$sew_cso)
  load_sew_wwtp_mean_quan <- do.call(meanQuantiles, args1$sew_wwtp)
  
  ## get mean and quantiles for load_rain_sum_path, load_sew_sum_path
  load_rain_sum_path_mean_quan <- do.call(meanQuantiles, args1$rain_sum)
  load_sew_sum_path_mean_quan  <- do.call(meanQuantiles, args1$sew_sum)
  
  # Define function arguments for combineLoads
  args2 <- list(
    
    cso  = list(variables, load_rain_cso, load_sew_cso),
    wwtp = list(variables, load_rain_wwtp, load_sew_wwtp),
    tot  = list(variables, load_rain_sum_path, load_sew_sum_path)
  )
  
  ## combine loads rainwater and sewage, load_rain_wwtp + load_sew_wwtp,
  ## rainwater + sewage
  load_cso_comb  <- do.call(combineLoads, args2$cso)
  load_wwtp_comb <- do.call(combineLoads, args2$wwtp)
  load_TOT       <- do.call(combineLoads, args2$tot)
  
  # Define function arguments for meanQuantiles (now that combined loads are 
  # available)
  args3 <- list(
    
    cso  = list(offset = 1, SUW_Names_rain, variables, load_cso_comb),
    wwtp = list(offset = 1, SUW_Names_rain, variables, load_wwtp_comb),
    tot  = list(offset = 0, SUW_Names_rain, variables, load_TOT)
  )
  
  # get mean and quantiles for load_cso, load_wwtp
  load_cso_mean_quan  <- do.call(meanQuantiles, args3$cso)
  load_wwtp_mean_quan <- do.call(meanQuantiles, args3$wwtp)
  
  # get mean and quantiles for total load
  load_TOT_mean_quan <- do.call(meanQuantiles, args3$tot)
  
  ### get mean and quantiles for volumes of rainwater and sewage
  ## rainwater volumes
  
  MC_vol_rai <- x_annual_loads_rain$MC_vol_rain[, -c(1, 2)]
  MC_vol_rain <- t(MC_vol_rai)
  
  vol_rain_mean_quan <- vol_rain[, 1:2]
  vol_rain_mean_quan[, 3:5] <- 0
  colnames(vol_rain_mean_quan)[3:5] <- c("mean", "Quan 5", "Quan 95")
  
  for (i in seq_len(3 * length(SUW_Names_rain))) {
    
    vol_rain_mean_quan[i, 3] <- mean(MC_vol_rain[, i])
    vol_rain_mean_quan[i, 4] <- quantile(MC_vol_rain[, i], probs = 0.05)
    vol_rain_mean_quan[i, 5] <- quantile(MC_vol_rain[, i], probs = 0.95)
  }
  
  # mean and quantiles for total rain volumes by SUW
  
  MC_vol_rain <- x_annual_loads_rain$MC_vol_rain
  
  vol_rain_list <- hsMatrixToListForm(
    MC_vol_rain, 
    keyFields = c("SUW", "Parameter"), 
    colNamePar = "run", 
    colNameVal = "volume"
  )
  
  sum_SUW <- aggregate(volume ~ SUW + run, data = vol_rain_list, FUN = sum)
  
  mean_SUW  <- aggregate(volume ~ SUW, data = sum_SUW, FUN = mean)
  mean_SUWt <- t(mean_SUW)
  colnames(mean_SUWt) <- mean_SUWt[1, ]

  quan5_SUW <- aggregate(volume ~ SUW, data = sum_SUW, FUN = quantile, 
                         probs = 0.05)
  quan5_SUWt <- t(quan5_SUW)
  colnames(quan5_SUWt) <- quan5_SUWt[1, ]
  
  quan95_SUW <- aggregate(volume ~ SUW, data = sum_SUW, FUN = quantile, 
                          probs = 0.95)
  quan95_SUWt <- t(quan95_SUW)
  colnames(quan95_SUWt) <- quan95_SUWt[1, ]

  # one data frame for mean and quantiles of rain volumes
  vol_rain_TOT <- toOverview(SUW_Names_rain, mean_SUWt, quant5_SUWt, quan95_SUWt)
  
  ## sewage volumes
  
  MC_vol_sew <- x_annual_loads_sew$MC_vol_sewage[, -c(1, 2)]
  MC_vol_sewage <- t(MC_vol_sew)
  
  vol_sewage_mean_quan <- vol_sewage[, 1:2]
  vol_sewage_mean_quan[, 3:5] <- 0
  colnames(vol_sewage_mean_quan)[3:5] <- c("mean", "Quan 5", "Quan 95")
  
  for (i in seq_len(2 * length(SUW_Names_sew))) {
    
    vol_sewage_mean_quan[i, 3] <- mean(MC_vol_sewage[, i])
    vol_sewage_mean_quan[i, 4] <- quantile(MC_vol_sewage[, i], probs = 0.05)
    vol_sewage_mean_quan[i, 5] <- quantile(MC_vol_sewage[, i], probs = 0.95)
  }
  
  # mean and quantiles for total volume of sewage by SUW
  
  MC_vol_sew <- x_annual_loads_sew$MC_vol_sewage
  
  vol_sew_list <- hsMatrixToListForm(
    MC_vol_sew, 
    keyFields = c("SUW", "Parameter"), 
    colNamePar = "run", 
    colNameVal = "volume"
  )
  
  sum_SUW <- aggregate(volume ~ SUW + run, data = vol_sew_list, FUN = sum)
  
  mean_SUW  <- aggregate(volume ~ SUW, data = sum_SUW, FUN = mean)
  mean_SUWt <- t(mean_SUW)
  colnames(mean_SUWt) <- mean_SUWt[1, ]
  
  quan5_SUW <- aggregate(volume ~ SUW, data = sum_SUW, FUN = quantile, 
                         probs = 0.05)
  quan5_SUWt <- t(quan5_SUW)
  colnames(quan5_SUWt) <- quan5_SUWt[1, ]
  
  quan95_SUW <- aggregate(volume ~ SUW, data = sum_SUW, FUN = quantile, 
                          probs = 0.95)
  quan95_SUWt <- t(quan95_SUW)
  colnames(quan95_SUWt) <- quan95_SUWt[1, ]
  
  # one data frame for mean and quantiles of sewage volumes
  vol_sew_TOT <- toOverview(SUW_Names_sew, mean_SUWt, quan5_SUWt, quan95_SUWt)
  
  ## get OgRe-dataframe-structure for plotting----------------------------------
  
  ## for loads in rainwater
  # by path
  
  loads_rain_by_path_mean_quan <- list()
  
  for (e in seq_along(SUW_Names_rain)) {
    
    loads_rain_by_path_mean_quan[[e]] <- data.frame(
      VariableName = x_conc_NEU$VariableName,
      SEP_mean = 0, 
      SEP_5 = 0, 
      SEP_95 = 0, 
      CSO_mean = 0, 
      CSO_5 = 0, 
      CSO_95 = 0,
      WWTP_mean = 0, 
      WWTP_5 = 0, 
      WWTP_95 = 0, 
      TOT_mean = 0, 
      TOT_5 = 0, 
      TOT_95 = 0
    )
    
    names(loads_rain_by_path_mean_quan)[e] <- SUW_Names_rain[e]
    
    for (f in seq_along(variables)) {
      
      loads_rain_by_path_mean_quan[[e]][f, 2] <- load_rain_sep_mean_quan[[f]][1, 1 + e]
      loads_rain_by_path_mean_quan[[e]][f, 3] <- load_rain_sep_mean_quan[[f]][2, 1 + e]
      loads_rain_by_path_mean_quan[[e]][f, 4] <- load_rain_sep_mean_quan[[f]][3, 1 + e]
      loads_rain_by_path_mean_quan[[e]][f, 5] <- load_rain_cso_mean_quan[[f]][1, 1 + e]
      loads_rain_by_path_mean_quan[[e]][f, 6] <- load_rain_cso_mean_quan[[f]][2, 1 + e]
      loads_rain_by_path_mean_quan[[e]][f, 7] <- load_rain_cso_mean_quan[[f]][3, 1 + e]
      loads_rain_by_path_mean_quan[[e]][f, 8] <- load_rain_wwtp_mean_quan[[f]][1, 1 + e]
      loads_rain_by_path_mean_quan[[e]][f, 9] <- load_rain_wwtp_mean_quan[[f]][2, 1 + e]
      loads_rain_by_path_mean_quan[[e]][f, 10] <- load_rain_wwtp_mean_quan [[f]][3, 1 + e]
      loads_rain_by_path_mean_quan[[e]][f, 11] <- load_rain_sum_path_mean_quan[[f]][1, 1 + e]
      loads_rain_by_path_mean_quan[[e]][f, 12] <- load_rain_sum_path_mean_quan[[f]][2, 1 + e]
      loads_rain_by_path_mean_quan[[e]][f, 13] <- load_rain_sum_path_mean_quan[[f]][3, 1 + e]
    }
  }
  
  ## for loads in sewage
  # by path
  
  loads_sew_by_path_mean_quan <- list()
  
  for (e in seq_along(SUW_Names_sew))
  {
    loads_sew_by_path_mean_quan[[e]] <- data.frame(
      VariableName = x_conc_NEU$VariableName,
      CSO_mean = 0, 
      CSO_5 = 0, 
      CSO_95 = 0, 
      WWTP_mean = 0, 
      WWTP_5 = 0, 
      WWTP_95 = 0,
      TOT_mean = 0, 
      TOT_5 = 0, 
      TOT_95 = 0
    )
    
    names(loads_sew_by_path_mean_quan)[e] <- SUW_Names_sew[e]
    
    for (f in seq_along(variables)) {
      
      loads_sew_by_path_mean_quan[[e]][f, 2] <- load_sew_cso_mean_quan[[f]][1, 1 + e]
      loads_sew_by_path_mean_quan[[e]][f, 3] <- load_sew_cso_mean_quan[[f]][2, 1 + e]
      loads_sew_by_path_mean_quan[[e]][f, 4] <- load_sew_cso_mean_quan[[f]][3, 1 + e]
      loads_sew_by_path_mean_quan[[e]][f, 5] <- load_sew_wwtp_mean_quan[[f]][1, 1 + e]
      loads_sew_by_path_mean_quan[[e]][f, 6] <- load_sew_wwtp_mean_quan[[f]][2, 1 + e]
      loads_sew_by_path_mean_quan[[e]][f, 7] <- load_sew_wwtp_mean_quan[[f]][3, 1 + e]
      loads_sew_by_path_mean_quan[[e]][f, 8] <- load_sew_sum_path_mean_quan[[f]][1, 1 + e]
      loads_sew_by_path_mean_quan[[e]][f, 9] <- load_sew_sum_path_mean_quan[[f]][2, 1 + e]
      loads_sew_by_path_mean_quan[[e]][f, 10] <- load_sew_sum_path_mean_quan[[f]][3, 1 + e]
    }
  }
  
  ## for summary of loads via CSO and WWTP
  
  loads_cso_wwtp_quan <- list()
  
  for (e in seq_along(SUW_Names_rain)) {
    
    loads_cso_wwtp_quan[[e]] <- data.frame(
      VariableName = x_conc_NEU$VariableName,
      CSO_5 = 0, 
      CSO_95 = 0, 
      WWTP_5 = 0, 
      WWTP_95 = 0
    )
    
    names(loads_cso_wwtp_quan)[e] <- SUW_Names_rain[e]
    
    for (f in seq_along(variables)) {
      
      loads_cso_wwtp_quan[[e]][f, 2] <- load_cso_mean_quan[[f]][2, 1 + e]
      loads_cso_wwtp_quan[[e]][f, 3] <- load_cso_mean_quan[[f]][3, 1 + e]
      loads_cso_wwtp_quan[[e]][f, 4] <- load_wwtp_mean_quan[[f]][2, 1 + e]
      loads_cso_wwtp_quan[[e]][f, 5] <- load_wwtp_mean_quan[[f]][3, 1 + e]
    }
  }
  
  ## for total loads
  loads_TOT_mean_quan <- summarise_loads(
    suwNames = SUW_Names_rain, 
    variables = x_conc_NEU$VariableName,
    totals = load_TOT_mean_quan
  )
}

### FUNCTIONS ###

# summarise_loads --------------------------------------------------------------
summarise_loads <- function(suwNames, variables, totals) 
{
  result <- list()
  
  for (i in seq_along(suwNames)) {
    
    x <- data.frame(
      VariableName = variables,
      TOT_mean = 0, 
      TOT_5 = 0, 
      TOT_95 = 0
    )
    
    for (row in seq_along(variables)) {
      
      x[row, 2] <- totals[[row]][1, 1 + i]
      x[row, 3] <- totals[[row]][2, 1 + i]
      x[row, 4] <- totals[[row]][3, 1 + i]
    }
    
    result[[suwNames[i]]] <- x
  }
  
  result
}

# combineLoads -----------------------------------------------------------------
combineLoads <- function(variables, x, y)
{
  result <- lapply(seq_along(variables), function(i) {
    
    xi <- x[[i]]
    yi <- y[[i]]
    
    for (column in SUW_Names_sew) {
      xi[, column] <- xi[, column] + yi[, column]
    }
    
    xi
  })
  
  # name the list elements and return
  structure(result, names = variables)
}

# meanQuantiles ----------------------------------------------------------------
# loads: loads in rainwater or sewage
meanQuantiles <- function(offset, suwNames, variables, loads)
{
  result <- list()
  
  for (a in seq_along(variables)) {
    
    result[[a]] <- data.frame(matrix(ncol = (1 + length(suwNames)), nrow = 3))
    
    result[[a]][, 1] <- c("mean", "Quan 5", "Quan 95")
    names(result)[a] <- variables[a]
    
    for (b in seq_along(suwNames)) {
      
      colnames(result[[a]]) <- c("Value", suwNames)
      
      x <- loads[[a]][, offset + b]
      
      result[[a]][1, 1 + b] <- mean(x)
      result[[a]][2, 1 + b] <- quantile(x, probs = 0.05)
      result[[a]][3, 1 + b] <- quantile(x, probs = 0.95) 
    }
  }
  
  result
}

# toOverview -------------------------------------------------------------------
toOverview <- function(suwNames, means, quantiles5, quantiles95)
{
  result <- data.frame(matrix(ncol = 1 + length(suwNames), nrow = 3))

  colnames(result) <- c("Values", suwNames)
  result[, 1] <- c("mean", "Quan 5", "Quan 95")
  
  for (column in suwNames) {
    
    result[1, column] <- means[2, column]
    result[2, column] <- quantiles5[2, column]
    result[3, column] <- quantiles95[2, column]
  }
  
  result
}
