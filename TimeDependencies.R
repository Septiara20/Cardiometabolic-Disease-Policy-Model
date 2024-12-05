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

restoreMsdataClass <- function(msdata, transitionMatrix) {
  if (is.null(attr(msdata, "trans"))) {
    attr(msdata, "trans") <- transitionMatrix
  }
  class(msdata) <- c("msdata", class(msdata))
  return(msdata)
}

buildContinuousData <- function(wideTransitionTable, observationDataset) {
    wideTransitionTable <- convertNAsToFalse(wideTransitionTable)
    wideTransitionTable <- replaceSpecialCharacters(wideTransitionTableContSub)
    wideTransitionTable <- addFactors(wideTransitionTable)
    msdata <- prepareDataForModel(transitionMatrix, wideTransitionTable, allCovariates)
    msdata <- addContinuousTimeDependentAge(msdata, wideTransitionTable)
    msdata <- addTimeDependentCovariates(msdata, observationDataset)
    msdata <- restoreMsdataClass(msdata, transitionMatrix)
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
    msdata <- restoreMsdataClass(msdata, transitionMatrix)
    return(msdata)
}

runContinuousCoxTimeDep <- function(msData, covariates, transition) {
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
    coxModel <- coxph(formula, data = filteredData, x = TRUE, method = "breslow")
    
    return(coxModel)
}

runCategoricalCoxTimeDep <- function(msData, covariates, transition) {
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
    
    # Get summary and round coefficients and confidence intervals
    summaryModelCox <- summary(coxModel)
    coefs <- summaryModelCox$coefficients
    coefs <- round(coefs, 3)
    
    # Add confidence intervals
    confInterval <- round(confint(coxModel), 3)
    
    # Combine coefficients with confidence intervals
    resultValueCox <- data.frame(
        Estimate = round(coefs[, "coef"],3),
        `Hazard Ratio` = round(exp(coefs[, "coef"]),3),
        `Std. Error` = round(coefs[, "se(coef)"], 3),
        `z value` = round(coefs[, "z"], 3),
        `Pr(>|z|)` = round(coefs[, "Pr(>|z|)"], 3),
        `Lower CI` = round(confInterval[, 1], 3),
        `Upper CI` = round(confInterval[, 2], 3)
    )
    
  return(list(model = coxModel, results = resultValueCox))
}

normalizePatidType <- function(data) {
  if (is.factor(data$patid)) {
    data[, patid := as.numeric(as.character(patid))]
  } else if (!is.numeric(data$patid)) {
    data[, patid := as.numeric(patid)]
  }
  return(data)
}

removeFactors <- function(data, covList) {
  # Loop through each column in the list
  for (col in covList) {
    if (col %in% names(data)) { # Check if the column exists in the data
      if (is.factor(data[[col]])) { # Check if the column is a factor
        data[[col]] <- as.character(data[[col]]) # Convert factor to character
      }
    } else {
      warning(glue("Column {col} not found in the dataset.")) # Warn if column doesn't exist
    }
  }
  return(data)
}
