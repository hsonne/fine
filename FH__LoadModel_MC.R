library(kwb.odm)
library(kwb.utils)
library(kwb.ogre)
require(dplyr)

# directories to be adapted by user:

# 1. directory of scripts
# script.dir <- "C:/Users/Fine/Desktop/Diplomarbeit/LoadModel"

# 2. directory of additional data
data.dir <- "C:/Users/Fine/Desktop/Diplomarbeit/LoadModel/data_LoadModel"

# MAIN
if (FALSE){
  
  # number of Monte Carlo simulations
  runs <- 1000 
  
  # 1. calculate loads of rainwater-based substances (for three pathways) 
  
  # proportion of wrong connections in seperate sewer system
  x = 0
  y = 1 - x
  
  x_annual_loads_rain <- annual_load_rain(data.dir = data.dir)
  
  # 2. calculate loads of sewage based substances 
  
  x_annual_loads_sew <- annual_load_sewage(data.dir = data.dir)
}

### FUNCTIONS ###

# readTableOrStop --------------------------------------------------------------
readTableOrStop <- function(data.dir, filename, type, csv2 = FALSE)
{
  file <- file.path(data.dir, filename)
  
  if (file.exists(file)) {
    
    if (csv2) {
      
      read.csv2(file = file, stringsAsFactors = FALSE)
      
    } else {
      
      read.table(file = file, sep = ";", dec = ".", stringsAsFactors = FALSE, 
                 header = TRUE)
    }
  } else stop(sprintf(
    "File with %s (%s) not found in data.dir", type, filename
  ))
}

# initMonteCarlo ---------------------------------------------------------------
initMonteCarlo <- function(x, runs, log = TRUE, set.names = TRUE, 
                           column.mean = "mean", column.sd = "sd", seed = NULL)
{
  result <- data.frame(matrix(ncol = nrow(x), nrow = runs))
  
  if (set.names) {
    colnames(result) <- x$VariableName
  }
  
  if (! is.null(seed)) {
    set.seed(seed)
  }
  
  for (i in seq_len(ncol(result))) {
    
    result[[i]] <- if (log) {
      
      rlnorm(n = runs, meanlog = x[i, column.mean], sdlog = x[i, column.sd])
      
    } else {
      
      rnorm(n = runs, mean = x[i, column.mean], sd = x[i, column.sd])
    }
  }
  
  result
}

# toNumeric --------------------------------------------------------------------
toNumeric <- function(x, columns)
{
  for (column in columns) {
    x[, column] <- as.numeric(x[, column])
  }
  
  x
}

# annual_load_rain -------------------------------------------------------------
annual_load_rain <- function # calculates the load for each substance
### separates pathways (rain runoff, CSO and WWTP)
(
  data.dir
  ### path of model data (annual mean concentrations "NEU_meanln_sdln.csv",
  #### mean concentrations with wrong connections "BKE_meanln_sdln.csv"
  ### rain runoff volumes "Vol_rain.csv, 
  ### removal at WWTP "substance_info.csv")
) 
{
  #load data
  x_conc_NEU <- readTableOrStop(
    data.dir, filename = "NEU_meanln_sdln.csv", 
    type = "annual mean concentrations of rainwater"
  )
  
  x_conc_BKE <- readTableOrStop(
    data.dir, filename = "BKE_meanln_sdln.csv",
    type = paste("annual mean concentrations of rainwater with wrong",
                 "connections")
  )
  
  vol_rain <- readTableOrStop(
    data.dir, filename = "Vol_rain.csv", 
    type = "rain runoff", csv2 = TRUE
  )
  
  removal_rates <- readTableOrStop(
    data.dir, filename = "substance_info.csv", 
    type = "removal rates at WWTP", csv2 = TRUE  
  )
  
  ### loads of rainwater based substances via separate sewer system and CSO
  
  # Step 1: Monte Carlo simulations to get concentrations in rainwater with 
  # proportion of wrong connections (x)
  
  MC_conc_rain <- initMonteCarlo(x = x_conc_NEU, runs = runs, seed = 0)
  
  MC_conc_rain_wrongcon <- initMonteCarlo(x = x_conc_BKE, runs = runs, seed = 1)
  
  MC_conc_rain_sep <- MC_conc_rain * y + x * MC_conc_rain_wrongcon
  
  #hist(log10(MC_conc_rain_wrongcon$`Escherichia coli`), breaks = 100)
  
  # Step 2: Monte Carlo simulations to get rain volumes
  
  vol_rain <- toNumeric(vol_rain, "sd")
  
  MC_vol_rain_1 <- initMonteCarlo(
    x = vol_rain, runs = runs, log = FALSE, set.names = FALSE, seed = 3
  )
  
  MC_vol_rain_t <- t(MC_vol_rain_1)
  MC_vol_rain_1 <- vol_rain[, 1:2]
  MC_vol_rain <- cbind(MC_vol_rain_1, MC_vol_rain_t)
  
  # Step 3: Calculation of loads in list, sep + CSO
  
  load_rain_sep <- getloadsforCSOorSEP(
    concentration = MC_conc_rain_sep,
    units = x_conc_NEU$UnitsAbbreviation, 
    volume = MC_vol_rain,
    parameter = "ROWvol, Trennsystem [m3/a]"
  )
  
  load_rain_cso <- getloadsforCSOorSEP(
    concentration = MC_conc_rain,
    units = x_conc_NEU$UnitsAbbreviation,
    volume = MC_vol_rain,
    parameter = "ROWvol, CSO [m3/a]"
  )
  
  ### loads of rainwater based substances via WWTP
  
  # Step 1: MC to get concentration in rainwater and rain volume calculation is 
  # already done (MC_conc_rain, MC_vol_rain)
  # Step 2: Monte Carlo simulations to get removal rates
  
  # missing removal rates (mean and sd) are set = 0
  removal_rates <- toNumeric(removal_rates, c("Retention_.", "Retention_sd"))

  removal_rates[is.na(removal_rates$Retention_.), 2] <- 0
  removal_rates[is.na(removal_rates$Retention_sd), 3] <- 0
  
  # get removal rates for substances in x_conc_NEU only (and in same order)
  removal_rates_red <- x_conc_NEU[, 1:2]
  
  indices <- match(removal_rates_red$VariableName, removal_rates$VariableName)
  
  removal_rates_red$mean <- as.numeric(removal_rates$Retention_.[indices])
  removal_rates_red$sd <- as.numeric(removal_rates$Retention_sd[indices])
  
  MC_removal_rates <- initMonteCarlo(
    x = removal_rates_red, runs = runs, log = FALSE, set.names = TRUE, seed = 4
  )
  
  # Step 3: Calculation of loads in list, WWTP
  
  load_rain_wwtp <- getloadsforWWTP(
    concentration = MC_conc_rain,
    units = x_conc_NEU$UnitsAbbreviation,
    volume = MC_vol_rain,
    parameter = "ROWvol, WWTP [m3/a]",
    removal = MC_removal_rates
  )
  
  # sum paths (in list)
  
  load_rain_sum_paths <- list()
  VariableNames <- x_conc_NEU$VariableName
  SUW_Names_rain <- unique(vol_rain$SUW)
  
  for (i in seq_along(VariableNames)) {
    
    load_rain_sum_paths[[i]] <- data.frame(matrix(
      ncol = length(SUW_Names_rain), 
      nrow = runs
    ))
    
    colnames(load_rain_sum_paths[[i]]) <- SUW_Names_rain
    names(load_rain_sum_paths)[i] <- colnames(MC_conc_rain)[i]
    
    for (j in seq_along(SUW_Names_rain)) {
      load_rain_sum_paths[[i]][j] <- (load_rain_cso[[i]][1 + j] 
                                      + load_rain_sep[[i]][1 + j] 
                                      + load_rain_wwtp[[i]][1 + j])
    }
  }
  
  # output
  list(
    load_rain_sep = load_rain_sep, 
    load_rain_cso = load_rain_cso, 
    load_rain_wwtp = load_rain_wwtp, 
    load_rain_sum_paths = load_rain_sum_paths,
    MC_vol_rain = MC_vol_rain
  )
}

# annual_load_sewage ----------------------------------------------------------- 
annual_load_sewage <- function # calculates the load for each substance
### separates pathways (CSO and WWTP)
(
  data.dir
  ### path of model data("Vol_sewage.csv",
  ### removal at WWTP "substance_info.csv",
  ### just for names "x_conc_NEU")
) 
  
{
  #load data
  x_conc_NEU <- readTableOrStop(
    data.dir, filename = "NEU_meanln_sdln.csv", 
    type = "annual mean concentrations of rainwater"
  )
  
  vol_sewage <- readTableOrStop(
    data.dir, filename = "Vol_sewage.csv", 
    type = "sewage runoff")
  
  sub_sew_info <- readTableOrStop(
    data.dir, filename = "substance_info.csv", 
    type = "substance information WWTP", csv2 = TRUE
  )
  
  ### loads of sewage based substances via CSO and WWTP
  
  # Step 1: Monte Carlo Simulations to get sewage volume
  
  MC_vol_sew <- initMonteCarlo(
    x = vol_sewage, runs = runs, log = FALSE, set.names = FALSE, seed = 2
  )

  MC_vol_sew_t <- t(MC_vol_sew)
  MC_vol_sew <- vol_sewage[, 1:2]
  MC_vol_sewage <- cbind(MC_vol_sew, MC_vol_sew_t)
  
  # Step 2: Monte Carlo simulations to get removal rates, concentrations of 
  # influent/effluent WWTP
  
  # read substance information
  columns <- c("Retention_.", "Retention_sd", "CinWWTP_calculated", 
               "CinWWTP_sd", "CoutWWTP", "CoutWWTP_sd")
  
  sub_sew_info <- toNumeric(sub_sew_info, columns)

  # set retention to zero, where information is lacking
  sub_sew_info$Retention_.[is.na(sub_sew_info$Retention_.)] <- 0
  sub_sew_info$Retention_sd[is.na(sub_sew_info$Retention_sd)] <- 0
  
  # MC to get retention
  MC_retention <- initMonteCarlo(
    x = sub_sew_info, runs = runs, log = FALSE, set.names = TRUE,
    column.mean = "Retention_.", column.sd = "Retention_sd", seed = 6
  )

  # MC to get concentrations in wastewater
  MC_conc_sew <- initMonteCarlo(
    x = sub_sew_info, runs = runs, log = TRUE, set.names = TRUE, 
    column.mean = "CinWWTP_calculated", column.sd = "CinWWTP_sd", seed = 7
  )
  
  # MC to get effluent concentrations - WWTP
  MC_Cout_WWTP <- initMonteCarlo(
    x = sub_sew_info, runs = runs, log = TRUE, set.names = TRUE, 
    column.mean = "CoutWWTP", column.sd = "CoutWWTP_sd", seed = 8
  )
  
  # Step 3: Calculation of loads in list, CSO + WWTP
  
  ## CSO
  load_sew_cso <- getloadsforCSOorSEP(
    concentration = MC_conc_sew,
    units = x_conc_NEU$UnitsAbbreviation, 
    volume = MC_vol_sewage,
    parameter = "ROWvol, CSO [m3/a]"
  )
  
  ##WWTP
  load_sew_wwtp <- getloadsforWWTP(
    concentration = MC_conc_sew, 
    units = x_conc_NEU$UnitsAbbreviation,
    volume = MC_vol_sewage,
    parameter = "ROWvol, WWTP [m3/a]",
    removal = MC_retention
  )
  
  # sum paths (in list)
  load_sew_sum_paths <- list()
  VariableNames <- colnames(MC_conc_sew)
  SUW_Names_sew = unique(vol_sewage$SUW)
  
  for (i in seq_along(VariableNames)) {
    
    load_sew_sum_paths[[i]] <- data.frame(matrix(
      ncol = length(SUW_Names_sew), 
      nrow = runs
    ))
    
    colnames(load_sew_sum_paths[[i]]) <- SUW_Names_sew
    names(load_sew_sum_paths)[i] <- colnames(MC_conc_sew)[i]
    
    for (j in seq_along(SUW_Names_sew)) {
      load_sew_sum_paths[[i]][j] <- (load_sew_cso[[i]][1 + j] 
                                     + load_sew_wwtp[[i]][1 + j])
    }  
  }
  
  # output
  list(
    load_sew_cso = load_sew_cso,
    load_sew_wwtp = load_sew_wwtp,
    load_sew_sum_paths = load_sew_sum_paths,
    MC_vol_sewage = MC_vol_sewage
  )
}

# changeunit--------------------------------------------------------------------
changeunit <- function(mytable)
{
  #mytable <-   load_x
  
  # if unit is mg / l
  if (unique(mytable[,"unit"]) == "mg/L") {
    
    # apply convert value to all column except for "unit"
    mytable[, -which(names(mytable) == "unit")] <- mytable[, -which(names(mytable) == "unit")] / 1000
    
  }
  
  # if unit is MPN/100 mL
  if (unique(mytable[,"unit"]) == "MPN/100 mL") {
    
    # apply convert value to all column except for "unit"
    mytable[, -which(names(mytable) == "unit")] <- mytable[, -which(names(mytable) == "unit")] * 10000
    
  }
  
  # if unit is PFU/100 mL
  if (unique(mytable[,"unit"]) == "PFU/100 mL") {
    
    # apply convert value to all column except for "unit"
    mytable[, -which(names(mytable) == "unit")] <- mytable[, -which(names(mytable) == "unit")] * 10000
    
  }
  
  mytable
}

# getLoads ---------------------------------------------------------------------
getLoads <- function
(
  concentration, 
  # concentration of rain (for CSO) or rain with wrongcons (for sep) or sewage (CSO)
  units = x_conc_NEU$UnitsAbbreviation, # just for names
  volume, # rainwater or sewage volume
  parameter,
  removal = NULL
)
{
  # calculate loads in list
  load_x <- list()
  
  for (e in seq_len(ncol(concentration))) {
    
    load_x[[e]] <- data.frame(unit = rep(units[e], times = nrow(concentration)))
    
    names(load_x)[e] <- colnames(concentration)[e]
    
    # get volume for each SUW
    SUW_Names <- unique(volume$SUW)
    
    for (f in seq_along(SUW_Names)) {
      
      indices <- which(volume$SUW == SUW_Names[f])
      vol_x_SUW <- volume[indices, ]
      
      load_x[[e]][[1 + f]] <- NA
      
      for (run in seq_len(runs)) {
        
        condition <- vol_x_SUW$Parameter == parameter
  
        load <- concentration[run, e] * vol_x_SUW[condition, 2 + run]
        
        if (! is.null(removal)) {
          load <- load * (1 - removal[run, e] / 100)
        }
        
        load_x[[e]][[1 + f]][run] <- load
          
        colnames(load_x[[e]])[1 + f] <- SUW_Names[f]
      }
    }
    
    changeunit(load_x[[e]])
  }
  
  load_x
}

# getloadsforCSOorSEP ----------------------------------------------------------
getloadsforCSOorSEP <- function
(
  concentration, # concentration of rain (for CSO) or 
  # rain with wrongcons (for sep) or sewage (CSO)
  units, # just for names
  volume, # rainwater or sewage volume
  parameter
)
{
  # concentration = MC_conc_rain
  # concentration = MC_conc_rain_sep
  # concentration = MC_conc_sew
  # volume = MC_vol_rain
  # volume = MC_vol_sewage
  # parameter = "ROWvol, Trennsystem [m3/a]"
  # parameter = "ROWvol, CSO [m3/a]"
  
  getLoads(concentration, units, volume, parameter)
}

# getloadsforWWTP --------------------------------------------------------------
getloadsforWWTP <- function
(
  concentration, # concentration of rainwater or sewage
  units, # just for names
  volume, # rainwater or sewage volume
  parameter,
  removal
)
{
  # concentration = MC_conc_rain
  # concentration = MC_Cin
  # volume = MC_vol_rain
  # volume = MC_vol_sewage
  # parameter = "ROWvol, WWTP [m3/a]"
  # removal = MC_removal_rates or MC_retention
  
  getLoads(concentration, units, volume, parameter, removal)
}
