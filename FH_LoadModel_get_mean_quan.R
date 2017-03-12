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
  vol_rain_mean_quan <- getMeanAndQuantiles(
    x = vol_rain, 
    monteCarlo = t(x_annual_loads_rain$MC_vol_rain[, -(1:2)]), 
    suwNames = SUW_Names_rain,
    multiple = 3
  )

  # mean and quantiles for total rain volumes by SUW
  vol_rain_TOT <- toTotal(
    x = x_annual_loads_rain$MC_vol_rain, 
    suwNames = SUW_Names_rain
  )

  ## sewage volumes
  vol_sewage_mean_quan <- getMeanAndQuantiles(
    x = vol_sewage[, 1:2], 
    monteCarlo = t(x_annual_loads_sew$MC_vol_sewage[, -(1:2)]), 
    suwNames = SUW_Names_sew, 
    multiple = 2
  )
  
  # mean and quantiles for total volume of sewage by SUW
  # one data frame for mean and quantiles of sewage volumes
  vol_sew_TOT <- toTotal(
    x = x_annual_loads_sew$MC_vol_sewage,
    suwNames = SUW_Names_sew
  )
  
  ## get OgRe-dataframe-structure for plotting----------------------------------
  
  ## for loads in rainwater
  # by path
  
  loads_rain_by_path_mean_quan <- summarise_loads(
    suwNames = SUW_Names_rain,
    inputs = list(
      SEP = load_rain_sep_mean_quan, CSO = load_rain_cso_mean_quan,
      WWTP = load_rain_wwtp_mean_quan, TOT = load_rain_sum_path_mean_quan
    ),
    columns = c(
      "SEP_mean", "SEP_5", "SEP_95", "CSO_mean", "CSO_5", "CSO_95",
      "WWTP_mean", "WWTP_5", "WWTP_95", "TOT_mean", "TOT_5", "TOT_95"
    )
  )

  ## for loads in sewage
  # by path
  
  loads_sew_by_path_mean_quan <- summarise_loads(
    suwNames = SUW_Names_sew,
    variables = variables,
    inputs = list(
      CSO = load_sew_cso_mean_quan, WWTP = load_sew_wwtp_mean_quan, 
      TOT = load_sew_sum_path_mean_quan
    ),
    columns = c(
      "CSO_mean", "CSO_5", "CSO_95", "WWTP_mean", "WWTP_5", "WWTP_95",
      "TOT_mean", "TOT_5", "TOT_95"
    )
  )
  
  ## for summary of loads via CSO and WWTP
  loads_cso_wwtp_quan <- summarise_loads(
    suwNames = SUW_Names_rain,
    variables = variables,
    inputs = list(CSO = load_cso_mean_quan, WWTP = load_wwtp_mean_quan),
    columns = c("CSO_5", "CSO_95", "WWTP_5", "WWTP_95")
  )

  ## for total loads
  loads_TOT_mean_quan <- summarise_loads(
    suwNames = SUW_Names_rain, 
    variables = variables,
    inputs = list(TOT = load_TOT_mean_quan),
    columns = c("TOT_mean", "TOT_5", "TOT_95")
  )
}

### FUNCTIONS ###

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

# getMeanAndQuantiles ----------------------------------------------------------
getMeanAndQuantiles <- function(x, monteCarlo, suwNames, multiple)
{
  result <- x
  result[, 3:5] <- 0
  colnames(result)[3:5] <- c("mean", "Quan 5", "Quan 95")

  for (i in seq_len(multiple * length(suwNames))) {

    values <- monteCarlo[, i]
      
    result[i, 3] <- mean(values)
    result[i, 4] <- quantile(values, probs = 0.05)
    result[i, 5] <- quantile(values, probs = 0.95)
  }
  
  result
}

# toTotal ----------------------------------------------------------------------
toTotal <- function(x, suwNames)
{
  x.long <- hsMatrixToListForm(
    x, 
    keyFields = c("SUW", "Parameter"), 
    colNamePar = "run", 
    colNameVal = "volume"
  )
  
  sum_SUW <- aggregate(volume ~ SUW + run, data = x.long, FUN = sum)
  
  # one data frame for mean and quantiles of rain volumes
  toOverview(
    suwNames = suwNames, 
    means = aggregateBySUW(sum_SUW, mean), 
    quantiles5 = aggregateBySUW(sum_SUW, quantile, probs = 0.05), 
    quantiles95 = aggregateBySUW(sum_SUW, quantile, probs = 0.95)
  )
}

# aggregateBySUW ---------------------------------------------------------------
aggregateBySUW <- function(data, FUN, ...) 
{
  result <- t(aggregate(volume ~ SUW, data = data, FUN = FUN, ...))
  colnames(result) <- result[1, ]
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

# summarise_loads --------------------------------------------------------------
summarise_loads <- function(suwNames, variables, inputs, columns) 
{
  # Initialise a data frame with one column: VariableName
  x.init <- data.frame(VariableName = variables)
  
  # Add columns filled with zeroes by calling cbind with do.call
  args <- as.list(rep(0, length(columns)))
  x.init <- do.call(cbind, c(list(x.init), structure(args, names = columns)))
  
  # Loop through the vector suwNames
  result <- lapply(seq_along(suwNames), function(i) {
    
    # Start with the initial data frame
    x <- x.init
    
    # Loop through the rows of the data frame (corresponding to the variables)
    for (row.to in seq_len(nrow(x))) {
      
      # Loop through the column names (excluding the first column)
      for (column in columns) {
        
        # Split the column name at the underscore
        parts <- strsplit(column, "_")[[1]]
        
        # Select the input data frame according to the first part of the column
        # name and select the row from which to take the value according to the
        # right part of the column name ("mean" = 1, "5" = 2, "95" = 3)
        input <- selectElements(inputs, parts[1])
        row.from <- match(parts[2], c("mean", "5", "95"))
        
        # Copy the value from the input data frame to the output data frame
        x[row.to, column] <- input[[row.to]][row.from, 1 + i]
      }
    }
    
    x
  })
  
  structure(result, names = suwNames)
}
