
### run FH_LoadModel first ###

# following lists are required:
# x_annual_loads_rain
# x_annual_loads_sew

## SOURCE ;)

# Define file types
types <- c(
  NEU_meanln_sdln.csv = "annual mean concentrations of rainwater",
  Vol_rain.csv = "rain runoff",
  Vol_sewage.csv = "sewage runoff"
)

# load data
filename <- "NEU_meanln_sdln.csv"
x_conc_NEU <- readTableOrStop(data.dir, filename, types[filename])

filename <- "Vol_rain.csv"
vol_rain <- readTableOrStop(data.dir, filename, types[filename], csv2 = TRUE)

filename <- "Vol_sewage.csv"
vol_sewage <- readTableOrStop(data.dir, filename, types[filename])

# get Names for Variables, SUWs and paths

VariableNames <- x_conc_NEU$VariableName
SUW_Names_rain <- unique(vol_rain$SUW)
SUW_Names_sew <- unique(vol_sewage$SUW)
paths_rain <- c("SEP", "CSO", "WWTP", "TOT")
paths_sew <- c("CSO", "WWTP", "TOT")

#MAIN
if(FALSE)
{
  ## get mean and quantiles for loads in rainwater and all pathways-------------
  # CSO
  
  load_rain_cso <- x_annual_loads_rain$load_rain_cso
  load_rain_cso_mean_quan <- MeanQuantilesforloads(myvolume = vol_rain,
                                                   x_conc_NEU,
                                                   myloads = load_rain_cso)
  
  # SEP
  load_rain_sep <- x_annual_loads_rain$load_rain_sep
  load_rain_sep_mean_quan <- MeanQuantilesforloads(myvolume = vol_rain,
                                                   x_conc_NEU,
                                                   myloads = load_rain_sep)
  # WWTP
  load_rain_wwtp <- x_annual_loads_rain$load_rain_wwtp
  load_rain_wwtp_mean_quan <- MeanQuantilesforloads(myvolume = vol_rain,
                                                    x_conc_NEU,
                                                    myloads = load_rain_wwtp)
  
  
  ## get mean and quantiles for load_rain_sum_path
  
  load_rain_sum_path <- x_annual_loads_rain$load_rain_sum_paths
  load_rain_sum_path_mean_quan <- MeanQuantilesforloadscomb(vol_rain,
                                                            x_conc_NEU,
                                                            myloads = load_rain_sum_path)
  
  ## get mean and quantiles for loads in sewage and all pathways----------------
  
  # CSO
  load_sew_cso <- x_annual_loads_sew$load_sew_cso
  load_sew_cso_mean_quan <- MeanQuantilesforloads(myvolume = vol_sewage,
                                                  x_conc_NEU,
                                                  myloads = load_sew_cso)
  
  # WWTP
  load_sew_wwtp <- x_annual_loads_sew$load_sew_wwtp
  load_sew_wwtp_mean_quan <- MeanQuantilesforloads(myvolume = vol_sewage,
                                                   x_conc_NEU,
                                                   myloads = load_sew_wwtp)
  
  ## get mean and quantiles for load_sew_sum_path
  
  load_sew_sum_path <- x_annual_loads_sew$load_sew_sum_paths
  load_sew_sum_path_mean_quan <- MeanQuantilesforloadscomb(vol_rain = vol_sewage,
                                                           x_conc_NEU,
                                                           myloads = load_sew_sum_path)
  
  ## combine loads rainwater and sewage
  
  # load_rain_cso + load_sew_cso
  
  load_cso_comb <- list()
  
  for (i in seq_along(VariableNames)){
    
    x <- load_rain_cso[[i]]
    
    for (column in SUW_Names_sew) {
      
      x[, column] <- x[, column] + load_sew_cso[[i]][, column]
    }
    
    element <- VariableNames[i]
    load_cso_comb[[element]] <- x
  }
  
  # get mean and quantiles for load_cso
  
  load_cso_mean_quan <- MeanQuantilesforloads(myvolume = vol_rain,
                                              x_conc_NEU,
                                              myloads = load_cso_comb)
  
  # load_rain_wwtp + load_sew_wwtp
  
  load_wwtp_comb <- list()
  
  for (i in seq_along(VariableNames)){
    
    x <- load_rain_wwtp[[i]]
    
    for (column in SUW_Names_sew) {
      
      x[, column] <- x[, column] + load_sew_wwtp[[i]][, column]
    }
    
    element <- VariableNames[i]
    load_wwtp_comb[[element]] <- x
  }
  
  # get mean and quantiles for load_wwtp
  
  load_wwtp_mean_quan <- MeanQuantilesforloads(myvolume = vol_rain,
                                               x_conc_NEU,
                                               myloads = load_wwtp_comb)
  
  # total load in rainwater + sewage
  
  load_TOT <- list()
  
  for (i in seq_along(VariableNames)){
    
    x <- load_rain_sum_path[[i]]
    
    for (column in SUW_Names_sew) {
      
      x[, column] <- x[, column] + load_sew_sum_path[[i]][, column]
    }
    
    element <- VariableNames[i]
    load_TOT[[element]] <- x
  }
  
  # get mean and quantiles for total load
  
  load_TOT_mean_quan <- MeanQuantilesforloadscomb(vol_rain,
                                                  x_conc_NEU,
                                                  myloads = load_TOT)
  
  ### get mean and quantiles for volumes of rainwater and sewage
  ## rainwater volumes
  
  MC_vol_rai <- x_annual_loads_rain$MC_vol_rain[, -c(1, 2)]
  MC_vol_rain <- t(MC_vol_rai)
  
  vol_rain_mean_quan <- vol_rain[, 1:2]
  vol_rain_mean_quan[, 3:5] <- 0
  colnames(vol_rain_mean_quan)[3:5] <- c("mean", "Quan 5", "Quan 95")
  
  for (i in 1:(length(SUW_Names_rain) * 3)){
    
    vol_rain_mean_quan[i, 3] <- mean(MC_vol_rain[, i])
    vol_rain_mean_quan[i, 4] <- quantile(MC_vol_rain[, i], probs = 0.05)
    vol_rain_mean_quan[i, 5] <- quantile(MC_vol_rain[, i], probs = 0.95)
  }
  
  # mean and quantiles for total rain volumes by SUW
  
  MC_vol_rain <- x_annual_loads_rain$MC_vol_rain
  
  vol_rain_list <- kwb.utils::hsMatrixToListForm(
    MC_vol_rain, 
    keyFields = c("SUW", "Parameter"), 
    colNamePar = "run", 
    colNameVal = "volume"
  )
  
  sum_SUW <- aggregate(volume ~ SUW + run, data = vol_rain_list, FUN = sum)
  
  mean_SUW  <- aggregate(volume ~ SUW, data = sum_SUW, FUN = mean)
  mean_SUWt <- t(mean_SUW)
  colnames(mean_SUWt) <- mean_SUWt[1, ]

  quan5_SUW <- aggregate(volume ~ SUW, data = sum_SUW, 
                         FUN = quantile, probs = 0.05)
  quan5_SUWt <- t(quan5_SUW)
  colnames(quan5_SUWt) <- quan5_SUWt[1, ]
  
  quan95_SUW <- aggregate(volume ~ SUW, data = sum_SUW, 
                          FUN = quantile, probs = 0.95)
  quan95_SUWt <- t(quan95_SUW)
  colnames(quan95_SUWt) <- quan95_SUWt[1, ]

  # one data frame for mean and quantiles of rain volumes
  
  vol_rain_TOT <- data.frame(matrix(ncol = 1 + length(SUW_Names_rain), 
                                    nrow = 3))
  colnames(vol_rain_TOT) <- c("Values", SUW_Names_rain)
  vol_rain_TOT[, 1] <- c("mean", "Quan 5", "Quan 95")
  
  for (column in SUW_Names_rain){
    
    vol_rain_TOT[1, column] <- mean_SUWt[2, column]
    vol_rain_TOT[2, column] <- quan5_SUWt[2, column]
    vol_rain_TOT[3, column] <- quan95_SUWt[2, column]
  }
  
  ## sewage volumes
  
  MC_vol_sew <- x_annual_loads_sew$MC_vol_sewage[, -c(1, 2)]
  MC_vol_sewage <- t(MC_vol_sew)
  
  vol_sewage_mean_quan <- vol_sewage[, 1:2]
  vol_sewage_mean_quan[, 3:5] <- 0
  colnames(vol_sewage_mean_quan)[3:5] <- c("mean", "Quan 5", "Quan 95")
  
  for (i in 1:(length(SUW_Names_sew) * 2)){
    
    vol_sewage_mean_quan[i, 3] <- mean(MC_vol_sewage[, i])
    vol_sewage_mean_quan[i, 4] <- quantile(MC_vol_sewage[, i], probs = 0.05)
    vol_sewage_mean_quan[i, 5] <- quantile(MC_vol_sewage[, i], probs = 0.95)
  }
  
  # mean and quantiles for total volume of sewage by SUW
  
  MC_vol_sew <- x_annual_loads_sew$MC_vol_sewage
  
  vol_sew_list <- kwb.utils::hsMatrixToListForm(
    MC_vol_sew, 
    keyFields = c("SUW", "Parameter"), 
    colNamePar = "run", 
    colNameVal = "volume"
  )
  
  sum_SUW <- aggregate(volume ~ SUW + run, data = vol_sew_list, FUN = sum)
  
  mean_SUW  <- aggregate(volume ~ SUW, data = sum_SUW, FUN = mean)
  mean_SUWt <- t(mean_SUW)
  colnames(mean_SUWt) <- mean_SUWt[1, ]
  
  quan5_SUW <- aggregate(volume ~ SUW, data = sum_SUW, 
                         FUN = quantile, probs = 0.05)
  quan5_SUWt <- t(quan5_SUW)
  colnames(quan5_SUWt) <- quan5_SUWt[1, ]
  
  quan95_SUW <- aggregate(volume ~ SUW, data = sum_SUW, 
                          FUN = quantile, probs = 0.95)
  quan95_SUWt <- t(quan95_SUW)
  colnames(quan95_SUWt) <- quan95_SUWt[1, ]
  
  # one data frame for mean and quantiles of sewage volumes
  
  vol_sew_TOT <- data.frame(matrix(ncol = 1 + length(SUW_Names_sew), 
                                    nrow = 3))
  colnames(vol_sew_TOT) <- c("Values", SUW_Names_sew)
  vol_sew_TOT[, 1] <- c("mean", "Quan 5", "Quan 95")
  
  for (column in SUW_Names_sew){
    
    vol_sew_TOT[1, column] <- mean_SUWt[2, column]
    vol_sew_TOT[2, column] <- quan5_SUWt[2, column]
    vol_sew_TOT[3, column] <- quan95_SUWt[2, column]
  }
  
  ## get OgRe-dataframe-structure for plotting----------------------------------
  
  ## for loads in rainwater
  # by path
  
  loads_rain_by_path_mean_quan <- list()
  
  for (e in seq_along(SUW_Names_rain))
  {
    loads_rain_by_path_mean_quan[[e]] <- data.frame("VariableName" = x_conc_NEU$VariableName,
                                                    SEP_mean = 0, SEP_5 = 0, 
                                                    SEP_95 = 0, CSO_mean = 0, 
                                                    CSO_5 = 0, CSO_95 = 0,
                                                    WWTP_mean = 0, WWTP_5 = 0, 
                                                    WWTP_95 = 0, TOT_mean = 0, 
                                                    TOT_5 = 0, TOT_95 = 0)
    
    names(loads_rain_by_path_mean_quan)[e] <- SUW_Names_rain[e]
    
    for (f in seq_along(VariableNames)){
      
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
    loads_sew_by_path_mean_quan[[e]] <- data.frame("VariableName" = x_conc_NEU$VariableName,
                                                   CSO_mean = 0, CSO_5 = 0, 
                                                   CSO_95 = 0, WWTP_mean = 0, 
                                                   WWTP_5 = 0, WWTP_95 = 0,
                                                   TOT_mean = 0, TOT_5 = 0, 
                                                   TOT_95 = 0)
    
    names(loads_sew_by_path_mean_quan)[e] <- SUW_Names_sew[e]
    
    for (f in seq_along(VariableNames)){
      
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
  
  for (e in seq_along(SUW_Names_rain))
  {
    loads_cso_wwtp_quan[[e]] <- data.frame("VariableName" = x_conc_NEU$VariableName,
                                           CSO_5 = 0, CSO_95 = 0, 
                                           WWTP_5 = 0, WWTP_95 = 0)
    
    names(loads_cso_wwtp_quan)[e] <- SUW_Names_rain[e]
    
    for (f in seq_along(VariableNames)){
      
      loads_cso_wwtp_quan[[e]][f, 2] <- load_cso_mean_quan[[f]][2, 1 + e]
      loads_cso_wwtp_quan[[e]][f, 3] <- load_cso_mean_quan[[f]][3, 1 + e]
      loads_cso_wwtp_quan[[e]][f, 4] <- load_wwtp_mean_quan[[f]][2, 1 + e]
      loads_cso_wwtp_quan[[e]][f, 5] <- load_wwtp_mean_quan[[f]][3, 1 + e]
    }
  }
  
  ## for total loads
  
  loads_TOT_mean_quan <- list()
  
  for (e in seq_along(SUW_Names_rain))
  {
    loads_TOT_mean_quan[[e]] <- data.frame("VariableName" = x_conc_NEU$VariableName,
                                           TOT_mean = 0, TOT_5 = 0, 
                                           TOT_95 = 0)
    
    names(loads_TOT_mean_quan)[e] <- SUW_Names_rain[e]
    
    for (f in seq_along(VariableNames)){
      
      loads_TOT_mean_quan[[e]][f, 2] <- load_TOT_mean_quan[[f]][1, 1 + e]
      loads_TOT_mean_quan[[e]][f, 3] <- load_TOT_mean_quan[[f]][2, 1 + e]
      loads_TOT_mean_quan[[e]][f, 4] <- load_TOT_mean_quan[[f]][3, 1 + e]
    }
  }
  
}

### FUNCTIONS###

### Get mean and quantiles for loads in rainwater or sewage and all pathways----
MeanQuantilesforloads <- function(myvolume, # just for names
                                  x_conc_NEU, # just for names
                                  myloads) # loads in rainwater or sewage
{
  mylist <- list()
  
  SUW_Names <- unique(myvolume$SUW)
  VariableNames <- x_conc_NEU$VariableName
  
  for (a in seq_along(VariableNames))
  {
    mylist[[a]] <- data.frame(matrix(ncol = (1 + length(SUW_Names)), 
                                     nrow = 3))
    
    mylist[[a]][, 1] <- c("mean", "Quan 5", "Quan 95")
    names(mylist)[a] <- VariableNames[a]
    
    for (b in seq_along(SUW_Names)){
      
      colnames(mylist[[a]]) <- c("Value", SUW_Names)
      
      x <- myloads[[a]][, 1 + b]
      
      mylist[[a]][1, 1 + b] <- mean(x)
      mylist[[a]][2, 1 + b] <- quantile(x, probs = 0.05)
      mylist[[a]][3, 1 + b] <- quantile(x, probs = 0.95) 
      
    }
  }
  
  mylist
  
}

### Get mean and quantiles for loads in different combinations----------------
MeanQuantilesforloadscomb <- function(vol_rain, # just for names
                                      x_conc_NEU, # just for names
                                      myloads) # combined loads
{
  mylist <- list()
  
  SUW_Names_rain <- unique(vol_rain$SUW)
  VariableNames <- x_conc_NEU$VariableName
  
  for (a in seq_along(VariableNames))
  {
    mylist[[a]] <- data.frame(matrix(ncol = (1 + length(SUW_Names_rain)), 
                                     nrow = 3))
    
    mylist[[a]][,1] <- c("mean", "Quan 5", "Quan 95")
    names(mylist)[a] <- VariableNames[a]
    
    for (b in seq_along(SUW_Names_rain)){
      
      colnames(mylist[[a]]) <- c("Value", SUW_Names_rain)
      
      x <- (myloads[[a]][, b])
      
      mylist[[a]][1, 1 + b] <- mean(x)
      mylist[[a]][2, 1 + b] <- quantile(x, probs = 0.05)
      mylist[[a]][3, 1 + b] <- quantile(x, probs = 0.95)
    }
  }
  
  mylist
  
}
