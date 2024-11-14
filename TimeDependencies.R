source("Biomarkers.R")
source("Continuous.R")

# Constants
################################################################################
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
  triglycerides = "triglyceride",
  cholesterol = "totalcholesterol",
  glucose = "fastingglucose",
  sbp = "sbp",
  dbp = "dbp",
  hba1c = "hba1c"
)

allCovariates <- c("gender", "deprivationIndex", "diabetesFH", "cvdFH", 
              "bmi", "hdl", "ldl", "triglycerides", "cholesterol", 
              "glucose", "sbp", "dbp", "latestSmokingStatus", 
              "atrialFib", "hyperlipidaemia", "hypertension", 
              "age", "ethnicity", "alcoholStatus")

biomarkerCovariates <- c("bmi", "hdl", "ldl", "triglycerides", "cholesterol", 
                         "glucose", "sbp", "dbp")
################################################################################

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

addTimeDependentCovariates <- function(msdata, observationDataset) {
  message("Applying time-dependent covariates...")

  message("CovariateMedcodes: ", paste(names(covariateMedcodes), collapse = ", "))
  message("ehrNamesMap: ", paste(names(ehrNamesMap), collapse = ", "))

  for (covariateName in names(covariateMedcodes)) {
    # Get the corresponding EHR name for cleaning
    ehrName <- ehrNamesMap[[covariateName]]
    if (is.null(ehrName)) {
      stop(glue("Mapping for {covariateName} is missing in ehrNamesMap. Check map!"))
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
  
  # Convert patid back to a factor to match original format
  msdata[, patid := as.factor(patid)]
  
  return(msdata)
}

buildContinuousData <- function(wideTransitionTable, observationDataset) {
    wideTransitionTable <- convertNAsToFalse(wideTransitionTable)
    wideTransitionTable <- replaceSpecialCharacters(wideTransitionTableContSub)
    wideTransitionTable <- addFactors(wideTransitionTable)
    msdata <- prepareDataForModel(transitionMatrix, wideTransitionTable, allCovariates)
    msdata <- addContinuousTimeDependentAge(msdata, wideTransitionTable)
    msdata <- addTimeDependentCovariates(msdata, observationDataset)
    return(msdata)
}

buildCategoricalData <- function(wideTransitionTable, observationDataset) {
    wideTransitionTable <- convertNAsToFalse(wideTransitionTable)
    wideTransitionTable <- replaceSpecialCharacters(wideTransitionTable)
    wideTransitionTable <- addFactors(wideTransitionTable, TRUE)
    msdata <- prepareDataForModel(transitionMatrix, wideTransitionTable, allCovariates)
    msdata <- addTimeDependentAge(msdata, wideTransitionTable)
    msdata <- addTimeDependentCovariates(msdata, observationDataset)
    msdata <- performCategorization(msdata)
    return(msdata)
}

runContinuousCoxTimeDep <- function(msData, transition) {
    message(glue("Fitting model for transition: {transition} and covariates: {paste(covariates, collapse = ', ')}"))
    setDT(msData)
    filteredData <- msData[trans == transition]
    filteredData <- as.data.frame(filteredData)
    gc() 
    rm(msData)
    gc()

    # Step 4: Construct model representing Cox survival function
    formula <- as.formula(paste("Surv(Tstart, Tstop, status) ~", paste(allCovariates, collapse = " + ")))
    
    # Remove factors (only categorical data needs factors)
    for (cov in biomarkerCovariates) {
      if (is.factor(filteredData[[cov]])) {
        filteredData[[cov]] <- as.numeric(as.character(filteredData[[cov]]))
      }
    }

    # Step 5: Fit Cox model with the expanded covariates
    filteredData$bmi <- round(filteredData$bmi, 1)
    coxModel <- coxph(formula, data = filteredData, x = TRUE, method = "breslow")
    
    return(coxModel)
}

runCategoricalCoxTimeDep <- function(msData, transition) {
    message(glue("Fitting model for transition: {transition} and covariates: {paste(covariates, collapse = ', ')}"))
    setDT(msData)
    filteredData <- msData[trans == transition]
    filteredData <- as.data.frame(filteredData)
    gc() 
    rm(msData)
    gc()

    # Step 4: Construct model representing Cox survival function
    formula <- as.formula(paste("Surv(Tstart, Tstop, status) ~", paste(allCovariates, collapse = " + ")))

    # Step 5: Fit Cox model with the expanded covariates
    coxModel <- coxph(formula, data = filteredData, x = TRUE, method = "breslow")
    
    return(coxModel)
}

normalizePatidType <- function(data) {
  if (is.factor(data$patid)) {
    data[, patid := as.numeric(as.character(patid))]
  } else if (!is.numeric(data$patid)) {
    data[, patid := as.numeric(patid)]
  }
  return(data)
}


