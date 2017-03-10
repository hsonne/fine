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
  file <- file.path(data.dir, "NEU_meanln_sdln.csv")
  if (file.exists(file)) {
    x_conc_NEU <- read.table(file = file, 
                             sep = ";", dec = ".", stringsAsFactors=FALSE, 
                             header = TRUE)
  } else stop("File with annual mean concentrations of rainwater
              (NEU_meanln_sdln.csv) not found in data.dir")
  
  file <- file.path(data.dir, "BKE_meanln_sdln.csv")
  if (file.exists(file)) {
    x_conc_BKE <- read.table(file = file, 
                             sep = ";", dec = ".", stringsAsFactors=FALSE, 
                             header = TRUE)
  } else stop("File with annual mean concentrations of rainwater with wrong 
              connections(BKE_meanln_sdln.csv) not found in data.dir")
  
  file <- file.path(data.dir, "Vol_rain.csv")
  if (file.exists(file)) {
    vol_rain <- read.csv2(file=file, 
                          stringsAsFactors=FALSE)
  } else stop("File with rain runoff (Vol_rain.csv) not found in data.dir")
  
  file <- file.path(data.dir, "substance_info.csv")
  if (file.exists(file)) {
    removal_rates <- read.csv2(file=file, 
                               header=TRUE, stringsAsFactors=FALSE)
  } else stop("File with removal rates at WWTP (substance_info.csv) not found 
              in data.dir")
  
  
  ### loads of rainwater based substances via separate sewer system and CSO
  
  # Step 1: Monte Carlo simulations to get concentrations in rainwater with 
  # proportion of wrong connections (x)
  ##set.seed(1) # reproducible result
  
  MC_conc_rain <- data.frame(matrix(ncol = nrow(x_conc_NEU), nrow = runs))
  colnames(MC_conc_rain) <- x_conc_NEU$VariableName
  
  set.seed(0)
  for(i in 1:ncol(MC_conc_rain)){
    MC_conc_rain[[i]] <- rlnorm(n = runs, 
                                meanlog = x_conc_NEU$mean[i], 
                                sdlog = x_conc_NEU$sd[i])
  }
  
  MC_conc_rain_wrongcon <- data.frame(matrix(ncol = nrow(x_conc_BKE), 
                                             nrow = runs))
  colnames(MC_conc_rain_wrongcon) <- x_conc_BKE$VariableName
  
  set.seed(1)
  for(i in 1:ncol(MC_conc_rain_wrongcon)){
    MC_conc_rain_wrongcon[[i]] <- rlnorm(n = runs, 
                                         meanlog = x_conc_BKE$mean[i], 
                                         sdlog = x_conc_BKE$sd[i])
  }
  
  MC_conc_rain_sep <- MC_conc_rain * y + x * MC_conc_rain_wrongcon
  
  #hist(log10(MC_conc_rain_wrongcon$`Escherichia coli`), breaks = 100)
  
  # Step 2: Monte Carlo simulations to get rain volumes
  
  MC_vol_rain_1 <- data.frame(matrix(ncol = nrow(vol_rain), nrow = runs))
  vol_rain$sd <- as.numeric(vol_rain$sd)
  
  set.seed(3)
  for(i in 1:ncol(MC_vol_rain_1)){
    MC_vol_rain_1[[i]] <- rnorm(runs, 
                                vol_rain$mean[i], 
                                vol_rain$sd[i])
  }
  
  MC_vol_rain_t <- t(MC_vol_rain_1)
  MC_vol_rain_1 <- vol_rain[, 1:2]
  MC_vol_rain <- cbind(MC_vol_rain_1, MC_vol_rain_t)
  
  # Step 3: Calculation of loads in list, sep + CSO
  
  load_rain_sep <- getloadsforCSOorSEP(myconc_MC = MC_conc_rain_sep,
                                       x_conc_NEU, 
                                       myvol_MC = MC_vol_rain,
                                       myrowname = "ROWvol, Trennsystem [m3/a]")
  
  load_rain_cso <- getloadsforCSOorSEP(myconc_MC = MC_conc_rain,
                                       x_conc_NEU,
                                       myvol_MC = MC_vol_rain,
                                       myrowname = "ROWvol, CSO [m3/a]")
  
  ### loads of rainwater based substances via WWTP
  
  # Step 1: MC to get concentration in rainwater and rain volume calculation is 
  # already done (MC_conc_rain, MC_vol_rain)
  # Step 2: Monte Carlo simulations to get removal rates
  
  # missing removal rates (mean and sd) are set = 0
  removal_rates$Retention_. <- as.numeric(removal_rates$Retention_.)
  removal_rates$Retention_sd <- as.numeric(removal_rates$Retention_sd)
  removal_rates[which(is.na(removal_rates$Retention_.)), 2] <- 0
  removal_rates[which(is.na(removal_rates$Retention_sd)), 3] <- 0
  
  # get removal rates for substances in x_conc_NEU only (and in same order)
  removal_rates_red <- x_conc_NEU[, 1:2]
  indices <- match(removal_rates_red$VariableName, removal_rates$VariableName)
  removal_rates_red$mean <- as.numeric(removal_rates$Retention_.[indices])
  removal_rates_red$sd <- as.numeric(removal_rates$Retention_sd[indices])
  
  MC_removal_rates <- data.frame(matrix(ncol = nrow(removal_rates_red), 
                                        nrow = runs))
  colnames(MC_removal_rates) <- removal_rates_red$VariableName
  
  set.seed(4)
  for(i in 1:ncol(MC_removal_rates)){
    MC_removal_rates[[i]] <- rnorm(n = runs, 
                                   mean = removal_rates_red$mean[i], 
                                   sd = removal_rates_red$sd[i])
  }
  
  # Step 3: Calculation of loads in list, WWTP
  
  load_rain_wwtp <- getloadsforWWTP(myconc_MC = MC_conc_rain,
                                    x_conc_NEU,
                                    myvol_MC = MC_vol_rain,
                                    myrowname = "ROWvol, WWTP [m3/a]",
                                    myremoval_MC = MC_removal_rates)
  
  # sum paths (in list)
  
  load_rain_sum_paths <- list()
  VariableNames <- x_conc_NEU$VariableName
  SUW_Names_rain <- unique(vol_rain$SUW)
  
  for (i in 1:length(VariableNames)){
    
    load_rain_sum_paths[[i]] <- data.frame((matrix(ncol = length(SUW_Names_rain), 
                                                   nrow = runs)))
    colnames(load_rain_sum_paths[[i]]) <- SUW_Names_rain
    names(load_rain_sum_paths)[i] <- colnames(MC_conc_rain)[i]
    
    for (j in 1:length(SUW_Names_rain))
    {
      load_rain_sum_paths[[i]][j] <- (load_rain_cso[[i]][1 + j] 
                                      + load_rain_sep[[i]][1 + j] 
                                      + load_rain_wwtp[[i]][1 + j])
    }  
    
  }
  
  # output
  list(load_rain_sep = load_rain_sep, 
       load_rain_cso = load_rain_cso, 
       load_rain_wwtp = load_rain_wwtp, 
       load_rain_sum_paths = load_rain_sum_paths,
       MC_vol_rain = MC_vol_rain)
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
  file <- file.path(data.dir, "NEU_meanln_sdln.csv")
  if (file.exists(file)) {
    x_conc_NEU <- read.table(file = file, 
                             sep = ";", dec = ".", stringsAsFactors=FALSE, 
                             header = TRUE)
  } else stop("File with annual mean concentrations of rainwater
              (NEU_meanln_sdln.csv) not found in data.dir")
  
  file <- file.path(data.dir, "Vol_sewage.csv")
  if (file.exists(file)) {
    vol_sewage <- read.table(file = file, 
                             sep = ";", dec = ".", stringsAsFactors=FALSE, 
                             header = TRUE)
  } else stop("File with sewage runoff (Vol_sewage.csv) not found in data.dir")
  
  file <- file.path(data.dir, "substance_info.csv")
  if (file.exists(file)) {
    sub_sew_info <- read.csv2(file=file, 
                              header=TRUE, stringsAsFactors=FALSE)
  } else stop("File with substance information WWTP (substance_info.csv) 
              not found in data.dir")
  
  ### loads of sewage based substances via CSO and WWTP
  
  # Step 1: Monte Carlo Simulations to get sewage volume
  
  MC_vol_sew <- data.frame(matrix(ncol = nrow(vol_sewage), nrow = runs))
  
  set.seed(2)
  for(i in 1:ncol(MC_vol_sew)){
    MC_vol_sew[[i]] <- rnorm(runs, 
                             vol_sewage$mean[i], 
                             vol_sewage$sd[i])
  }
  
  MC_vol_sew_t <- t(MC_vol_sew)
  MC_vol_sew <- vol_sewage[, 1:2]
  MC_vol_sewage <- cbind(MC_vol_sew, MC_vol_sew_t)
  
  # Step 2: Monte Carlo simulations to get removal rates, concentrations of 
  # influent/effluent WWTP
  
  # read substance information
  sub_sew_info$Retention_. <- as.numeric(sub_sew_info$Retention_.)
  sub_sew_info$Retention_sd <- as.numeric(sub_sew_info$Retention_sd) 
  sub_sew_info$CinWWTP_calculated <- as.numeric(sub_sew_info$CinWWTP_calculated)
  sub_sew_info$CinWWTP_sd <- as.numeric((sub_sew_info$CinWWTP_sd))
  sub_sew_info$CoutWWTP <- as.numeric(sub_sew_info$CoutWWTP)
  sub_sew_info$CoutWWTP_sd <- as.numeric(sub_sew_info$CoutWWTP_sd)
  
  # set retention to zero, where information is lacking
  indices <- which(is.na(sub_sew_info$Retention_.))
  sub_sew_info$Retention_.[indices] <- 0
  indices2 <- which(is.na(sub_sew_info$Retention_sd))
  sub_sew_info$Retention_sd[indices2] <- 0
  
  # MC to get retention
  MC_retention <- data.frame(matrix(ncol = nrow(sub_sew_info), nrow = runs))
  colnames(MC_retention) <- sub_sew_info$VariableName
  
  set.seed(6)
  for(i in 1:ncol(MC_retention)){
    MC_retention[[i]] <- rnorm(n = runs, 
                               mean = sub_sew_info$Retention_.[i], 
                               sd = sub_sew_info$Retention_sd[i])
  }
  
  # MC to get concentrations in wastewater
  MC_conc_sew <- data.frame(matrix(ncol = nrow(sub_sew_info), nrow = runs))
  colnames(MC_conc_sew) <- sub_sew_info$VariableName
  
  set.seed(7)
  for(i in 1:ncol(MC_conc_sew)){
    
    MC_conc_sew[[i]] <- rlnorm(n = runs,
                               meanlog = sub_sew_info$CinWWTP_calculated[i],
                               sdlog = sub_sew_info$CinWWTP_sd[i]) 
  }
  
  # MC to get effluent concentrations - WWTP
  MC_Cout_WWTP <- data.frame(matrix(ncol = nrow(sub_sew_info), nrow = runs))
  colnames(MC_Cout_WWTP) <- sub_sew_info$VariableName
  
  set.seed(8)
  for (i in 1:ncol(MC_Cout_WWTP)){
    MC_Cout_WWTP[[i]] <- rlnorm(n = runs,
                                meanlog = sub_sew_info$CoutWWTP[i],
                                sdlog = sub_sew_info$CoutWWTP_sd[i])
  }
  
  # Step 3: Calculation of loads in list, CSO + WWTP
  
  ## CSO
  load_sew_cso <- getloadsforCSOorSEP(myconc_MC = MC_conc_sew,
                                      x_conc_NEU, 
                                      myvol_MC = MC_vol_sewage,
                                      myrowname = "ROWvol, CSO [m3/a]")
  
  ##WWTP
  load_sew_wwtp <- getloadsforWWTP(myconc_MC = MC_conc_sew, 
                                   x_conc_NEU,
                                   myvol_MC = MC_vol_sewage,
                                   myrowname = "ROWvol, WWTP [m3/a]",
                                   myremoval_MC = MC_retention)
  
  # sum paths (in list)
  load_sew_sum_paths <- list()
  VariableNames <- colnames(MC_conc_sew)
  SUW_Names_sew = unique(vol_sewage$SUW)
  
  for (i in 1:length(VariableNames))
  {
    load_sew_sum_paths[[i]] <- data.frame((matrix(ncol = length(SUW_Names_sew), 
                                                  nrow = runs)))
    
    colnames(load_sew_sum_paths[[i]]) <- SUW_Names_sew
    names(load_sew_sum_paths)[i] <- colnames(MC_conc_sew)[i]
    
    for (j in 1:length(SUW_Names_sew))
    {
      load_sew_sum_paths[[i]][j] <- (load_sew_cso[[i]][1 + j] 
                                     + load_sew_wwtp[[i]][1 + j])
    }  
  }
  
  # output
  list(load_sew_cso = load_sew_cso,
       load_sew_wwtp = load_sew_wwtp,
       load_sew_sum_paths = load_sew_sum_paths,
       MC_vol_sewage = MC_vol_sewage)
  
}

# changeunit--------------------------------------------------------------------
changeunit <- function(mytable){
  
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


# getloadsforCSOorSEP ----------------------------------------------------------
getloadsforCSOorSEP <- function(myconc_MC, # concentration of rain (for CSO) or 
                                # rain with wrongcons (for sep) or sewage (CSO)
                                x_conc_NEU, # just for names
                                myvol_MC, # rainwater or sewage volume
                                myrowname)
{
  
  # myconc_MC = MC_conc_rain
  # myconc_MC = MC_conc_rain_sep
  # myconc_MC = MC_conc_sew
  # myvol_MC = MC_vol_rain
  # myvol_MC = MC_vol_sewage
  # myrowname = "ROWvol, Trennsystem [m3/a]"
  # myrowname = "ROWvol, CSO [m3/a]"
  
  # calculate loads in list
  load_x <- list()
  
  for(e in 1:ncol(myconc_MC)){
    
    load_x[[e]] <- data.frame("unit" = rep(x_conc_NEU$UnitsAbbreviation[e], 
                                           times = nrow(myconc_MC)))
    names(load_x)[e] <- colnames(myconc_MC)[e]
    
    # get volume for each SUW
    SUW_Names <- unique(myvol_MC$SUW)
    for (f in 1:length(SUW_Names)){
      
      indices <- which(myvol_MC$SUW == SUW_Names[f])
      vol_x_SUW <- myvol_MC[indices, ]
      load_x[[e]][[1 + f]] <- NA
      
      for(run in 1:runs){ 
        
        condition <- vol_x_SUW$Parameter == myrowname 
        load_x[[e]][[1 + f]][run] <- myconc_MC[run, e] * vol_x_SUW[condition, 
                                                                   2 + run]
        colnames(load_x[[e]])[1 + f] <- SUW_Names[f]
      }
    }
    changeunit(load_x[[e]])
  }
  
  load_x
  
}

# getloadsforWWTP --------------------------------------------------------------
getloadsforWWTP <- function(myconc_MC, # concentration of rainwater or sewage
                            x_conc_NEU, # just for names
                            myvol_MC, # rainwater or sewage volume
                            myrowname,
                            myremoval_MC)
{
  
  # myconc_MC = MC_conc_rain
  # myconc_MC = MC_Cin
  # myvol_MC = MC_vol_rain
  # myvol_MC = MC_vol_sewage
  # myrowname = "ROWvol, WWTP [m3/a]"
  # myremoval_MC = MC_removal_rates or MC_retention
  
  # calculate loads in list
  load_x <- list()
  
  for(e in 1:ncol(myconc_MC)){
    
    load_x[[e]] <- data.frame("unit" = rep(x_conc_NEU$UnitsAbbreviation[e], 
                                           times = nrow(myconc_MC)))
    names(load_x)[e] <- colnames(myconc_MC)[e]
    
    # get volume for each SUW
    SUW_Names <- unique(myvol_MC$SUW)
    
    for (f in 1:length(SUW_Names)){
      
      indices <- which(myvol_MC$SUW == SUW_Names[f])
      vol_x_SUW <- myvol_MC[indices, ]
      load_x[[e]][[1 + f]] <- NA
      
      for(run in 1:runs){ 
        
        condition <- vol_x_SUW$Parameter == myrowname 
        load_x[[e]][[1 + f]][run] <- (myconc_MC[run, e] * vol_x_SUW[condition, 
                                                                    2 + run] 
                                    * (1 - myremoval_MC[run, e] / 100))
        
        colnames(load_x[[e]])[1 + f] <- SUW_Names[f]
      }
    }
    
    changeunit(load_x[[e]])
    
  }
  
  load_x
  
}
