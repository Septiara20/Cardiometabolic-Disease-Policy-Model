########################################################
# Experimental code to run Cox for continuous covariates
########################################################

# So we don't forget to load the other files
source("Biomarkers.R")
source("SurvivalAnalysis.R")

################################################################################
# Usage : 
#runContinuousCox(wideTransitionTableContSub, c("hdl", "ldl", "atrialFib", "bmi", "gender", "age"), 1)
################################################################################


# Operates on a specified transition given a covariate list and the
# wideTransitionTableContSub
runContinuousCox <- function(wideTransitionTable, covariates, transition) {

    wideTransitionTable <- convertNAsToFalse(wideTransitionTable)
    
    # Step 1: Prepare data in long format using prepareDataForModel (from SurvivalAnalysis.R)
    preparedData <- prepareDataForModel(transitionMatrix, wideTransitionTable, covariateList = covariates)
    
    # Step 2: Add age if necessary
    if ("age" %in% covariates) {
        preparedData <- addContinuousTimeDependentAge(preparedData, wideTransitionTable)
    }
    
    # Step 3: Extract only the data for the given transition
    message(glue("Fitting model for transition: {transition} and covariates: {paste(covariates, collapse = ', ')}"))
    setDT(preparedData)
    preparedData <- preparedData[trans == transition]
    filteredData <- as.data.frame(preparedData)
    gc() 
    rm(preparedData)
    rm(wideTransitionTable)
    gc()

    # Step 4: Construct model representing Cox survival function
    formula <- as.formula(paste("Surv(Tstart, Tstop, status) ~", paste(covariates, collapse = " + ")))
    
    # Remove factors (categorical data needs factors)
    for (cov in covariates) {
      if (is.factor(filteredData[[cov]])) {
        filteredData[[cov]] <- as.numeric(as.character(filteredData[[cov]]))
      }
    }

    # Step 5: Fit Cox model with the expanded covariates
    filteredData$bmi <- round(filteredData$bmi, 1)
    coxModel <- coxph(formula, data = filteredData, x = TRUE, method = "breslow")
    
    return(coxModel)
}

addContinuousTimeDependentAge <- function(msdata, wideTransitionTable) {
  # Merge the age from wideTransitionTable into msdata based on patid
  # Assuming 'age' is the baseline age in wideTransitionTable

  # Calculate age at each Tstart (in years) and store it as the updated age
  msdata$age <- msdata$age + (msdata$Tstart / 365.25)
  
  # Clean up intermediate columns if not needed
  if ("age_at_Tstart" %in% names(msdata)) {
    msdata <- msdata[, !names(msdata) %in% c("age_at_Tstart")]
  }
  
  return(msdata)
}

#Run cox model for each transition
runContinuousCox(wideTransitionTableContSub, c("gender", "deprivationIndex", "diabetesFH", "cvdFH", 
                "bmi", "hdl", "triglycerides", "cholesterol", 
                "glucose", "sbp", "latestSmokingStatus", "alcoholStatus",
                "atrialFib", "hyperlipidaemia"), 1)

runContinuousCox(wideTransitionTableContSub, c("gender", "deprivationIndex", "diabetesFH", "cvdFH", 
                "bmi", "hdl", "triglycerides", "cholesterol", 
                "glucose", "sbp", "latestSmokingStatus", "alcoholStatus",
                "atrialFib", "hyperlipidaemia"), 2)

runContinuousCox(wideTransitionTableContSub, c("gender", "deprivationIndex", "diabetesFH", "cvdFH", 
                "bmi", "hdl", "triglycerides", "cholesterol", 
                "glucose", "sbp", "latestSmokingStatus", "alcoholStatus",
                "atrialFib", "hyperlipidaemia"), 3)

