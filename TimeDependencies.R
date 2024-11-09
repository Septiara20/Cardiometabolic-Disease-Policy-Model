source("Biomarkers.R")

# Define the medcode IDs for each covariate based on your mappings

# Use covariate names that directly match msdata
covariateMedcodes <- list(
  bmi = medcodeBMI$medcodeid,
  hdl = medcodeHDL$medcodeid,
  ldl = medcodeLDL$medcodeid,
  triglycerides = medcodeTriglycerides$medcodeid,
  cholesterol = medcodeTotalChol$medcodeid,
  glucose = medcodeFastingGlucose$medcodeid,
  sbp = medcodeSBP$medcodeid,
  dbp = medcodeDBP$medcodeid,
  hba1c = medcodeHBA1C$medcodeid
)

ehrNamesMap <- list(
  bmi = "bmi",
  hdl = "hdl",
  ldl = "ldl",
  triglycerides = "triglyceride", # Fix name mismatch
  cholesterol = "totalcholesterol",
  glucose = "fastingglucose",
  sbp = "sbp",
  dbp = "dbp",
  hba1c = "hba1c"
)

# Updated function to fetch the most recent value of a covariate
getMostRecentCovariateValues <- function(observationDataset, covariateName, medcodeID, msdata, ehrName) {
  message(glue("Fetching most recent values for {covariateName}"))
  
  # Ensure msdata is a data.table and normalize patid types
  setDT(msdata)
  msdata <- normalizePatidType(msdata)
  msdata[, Tstop_date := as.Date("1990-01-01") + Tstop]
  
  # Load observation data for the given covariate and normalize patid
  covariateData <- observationDataset %>%
    filter(patid %in% unique(msdata$patid), medcodeid %in% medcodeID) %>%
    collect() %>%
    as.data.table()
  covariateData <- normalizePatidType(covariateData)
  covariateData <- covariateData[obsdate >= as.Date("1990-01-01")]
  
  # Clean values using EHRBiomarkr
  message(glue("Number of {ehrName} biomarkers before cleaning: {nrow(covariateData)}"))
  covariateData <- clean_biomarker_values(covariateData, value, ehrName)
  message(glue("Number of {ehrName} biomarkers after cleaning: {nrow(covariateData)}"))
  
  # Sort and key data
  setorder(covariateData, patid, obsdate)
  setkey(covariateData, patid, obsdate)
  setorder(msdata, patid, Tstop_date)
  setkey(msdata, patid, Tstop_date)
  
  # Join and find most recent values
  result <- covariateData[msdata, on = .(patid), allow.cartesian = TRUE][
    obsdate <= Tstop_date, .SD[which.max(obsdate)], by = .(patid, Tstop_date)]
  
  # Ensure all rows from msdata are preserved
  if (covariateName %in% names(msdata)) {
    msdata[, (covariateName) := NULL]
  }
  
  msdata <- merge(
    msdata,
    result[, .(patid, Tstop_date, value)], # Keep only relevant columns
    by = c("patid", "Tstop_date"),
    all.x = TRUE
  )
  
  # Rename column to match covariateName
  setnames(msdata, "value", covariateName)
  
  # Clean up intermediate columns
  msdata[, Tstop_date := NULL]
  
  return(msdata)
}

# Updated function to apply time-dependent covariates
applyTimeDependentCovariates <- function(msdata, observationDataset, covariateMedcodes, ehrNamesMap) {
  message("Applying time-dependent covariates...")

  # Debug: Check mappings before starting
  message("CovariateMedcodes: ", paste(names(covariateMedcodes), collapse = ", "))
  message("ehrNamesMap: ", paste(names(ehrNamesMap), collapse = ", "))

  # Process each covariate
  for (covariateName in names(covariateMedcodes)) {
    # Get the corresponding EHR name for cleaning
    ehrName <- ehrNamesMap[[covariateName]]
    if (is.null(ehrName)) {
      stop(glue("Mapping for {covariateName} is missing in ehrNamesMap. Check your map!"))
    }
    
    message(glue("Processing covariate: {covariateName} (cleaning as {ehrName})"))
    msdata <- getMostRecentCovariateValues(
      observationDataset,
      covariateName,
      covariateMedcodes[[covariateName]],
      msdata,
      ehrName
    )
  }

  message("Completed applying time-dependent covariates.")
  
  # Convert patid back to a factor if necessary for downstream compatibility
  msdata[, patid := as.factor(patid)]
  
  return(msdata)
}

addTimeDependentAge <- function(msdata, wideTransitionTable) {
  # Merge the age from wideTransitionTable into msdata based on patid
  #msdata <- merge(msdata, wideTransitionTable[, .(patid, age)], by = "patid", all.x = TRUE)
  
  # Calculate age at each Tstart (in years)
  msdata$age_at_Tstart <- msdata$age + (msdata$Tstart / 365.25)

  # Define the age categories based on the updated age
  msdata$age <- cut(msdata$age_at_Tstart, 
                    breaks = c(-Inf, 24, 34, 44, 54, 64, Inf),
                    labels = c("18-24", "25-34", "35-44", "45-54", "55-64", ">65"),
                    right = FALSE)  # Adjust to include the lower bound correctly

  # Check for any NAs and handle them
  if (any(is.na(msdata$age))) {
    cat("Warning: NAs found in age after cutting. Consider reviewing the age ranges or input data.\n")
  }

  # Clean up intermediate columns if not needed
  msdata <- msdata[, !names(msdata) %in% c("age_at_Tstart")]

  return(msdata)
}