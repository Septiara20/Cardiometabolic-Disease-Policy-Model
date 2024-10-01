################################################################################
# Cardiometabolic disease policy model in the UK setting
# Septiara Putri 
# Health Economics and Health Technology Assessment (HEHTA), 
# University of Glasgow, UK
# 2024
################################################################################
# Biomarker descriptive statistic extractor
################################################################################
# Dependencies
################################################################################
library(arrow)
library(data.table)
library(lubridate)
library(dplyr)
library(glue)
library(bit64)
library(devtools)
library(parallel)
library(foreach)
library(doParallel)
library(forcats)
#install_github("Exeter-Diabetes/EHRBiomarkr")
library(EHRBiomarkr)
observationDataset <- open_dataset('Data/Observation/')
################################################################################
# Functions
################################################################################

####################################################################
# Simple data extraction functions
# These are necessary precursors to baseline table generation
####################################################################

###############################################################################
# Biomarker code called from SurvivalAnalysis.R is in this section

# Gets the date of the first CMD event for patients who have a CMD event
getFirstCMDDates <- function(filteredPatients, sampledPatient) {
  
  # Combine medical codes into one vector for efficient filtering
  combinedMedcodes <- unique(c(medcodeStroke$medcodeid, medcodeDiabetesT2$medcodeid, medcodeMyoInf$medcodeid))
  
  # Combine ICD codes into one vector for efficient filtering
  combinedICDCodes <- unique(c(icdCodeStroke$icd, icdCodeMyoInf$icd, icdCodeT2Diab$icd))
  
  # Pre-filter relevant observations and diagnosis, adding a source column
  relevantObservations <- observationDataset %>%
    filter(patid %in% filteredPatients & medcodeid %in% medcodeDiabetesT2$medcodeid) %>%
    select(patid, obsdate) %>%
    mutate(obsdate = as.Date(obsdate, format = "%Y-%m-%d"), source = "GP") %>%
    collect()

  relevantHospEpi <- sampledDiagEpi[patid %in% filteredPatients & ICD %in% combinedICDCodes, .(patid, epistart = as.Date(epistart, format = "%Y-%m-%d"), source = "Hospital")]
  
  # Convert to data.table for fast operations
  setDT(relevantObservations)
  setDT(relevantHospEpi)
  
  # Rename columns for consistency
  setnames(relevantHospEpi, "epistart", "date")
  
  # Set keys for fast joining and filtering
  setkey(relevantObservations, patid)
  setkey(relevantHospEpi, patid)
  
  # Combine the datasets
  combinedEvents <- rbindlist(list(
    relevantObservations[, .(patid, date = obsdate, source)],
    relevantHospEpi[, .(patid, date, source)]
  ), use.names = TRUE)

  # Find the earliest event date for each patient, including the source
  eventDates <- combinedEvents[, .SD[which.min(date)], by = patid]
  
  # Optional: Handle patients with all NA dates to ensure they are included with NA dates and source
  allPatients <- unique(c(filteredPatients, relevantObservations$patid, relevantHospEpi$patid))
  allPatientsData <- data.table(patid = allPatients)
  setkey(allPatientsData, patid)
  setkey(eventDates, patid)
  
  # Join to ensure all patients are included, even those with all NA dates and determine the source
  finalEventDates <- allPatientsData[eventDates, on = "patid"]
  
  sampledPatientClone <- copy(sampledPatient)
  # Define the start date of the study
  startDate <- as.Date("1990-01-01")
  studyStartYear <- as.numeric(format(startDate, "%Y"))
  
  # Calculate age at the start of the study
  sampledPatientClone[, age := studyStartYear - yob]
  
  # Merge with sampledPatient to get age, patId pair
  finalEventDates <- merge(finalEventDates, sampledPatientClone[, .(patid, age)], by = "patid", all.x = TRUE)
  
  finalEventDates <- convertAgeToCategory(finalEventDates)
  
  return(finalEventDates)
}

getBiomarkersWithin5Years <- function(firstCMDDates) {
  # Define study start date
  studyStartDate <- as.Date("1990-01-01")
  
  # Exclude AlcoholStatus from the biomarker list
  combinedBiomarkerMedcodes <- unique(c(
    medcodeAtrialFib$medcodeid, 
    medcodeBMI$medcodeid, 
    medcodeDBP$medcodeid, 
    medcodeFastingGlucose$medcodeid, 
    medcodeHBA1C$medcodeid, 
    medcodeHDL$medcodeid, 
    medcodeHyperlipidaemia$medcodeid, 
    medcodeHypertension$medcodeid, 
    medcodeLDL$medcodeid, 
    medcodeSBP$medcodeid, 
    medcodeSmokingStatus$medcodeid, 
    medcodeTotalChol$medcodeid, 
    medcodeTriglycerides$medcodeid
  ))
  
  combinedBiomarkerMedcodes <- as.integer64(combinedBiomarkerMedcodes)
  
  # Filter observationDataset for relevant biomarkers and patients
  relevantBioObservations <- observationDataset %>%
    filter(patid %in% firstCMDDates$patid, 
           medcodeid %in% combinedBiomarkerMedcodes) %>%
    collect()

  # Convert to data.table for faster manipulation
  relevantBioObservations <- as.data.table(relevantBioObservations)
  
  # Convert firstCMDDates to data.table and merge to get first CMD date for each patient
  firstCMDDates <- as.data.table(firstCMDDates)
  relevantBioObservations <- merge(relevantBioObservations, firstCMDDates, by = "patid")
  
  # Calculate the 5-year window for each observation
  relevantBioObservations[, fiveYearsBeforeCMD := pmax(date - years(5), studyStartDate)]
  
  # Filter to only include observations within the 5-year window before CMD and after study start date
  relevantBioObservations <- relevantBioObservations[obsdate >= fiveYearsBeforeCMD & obsdate <= date]
  
  if ("source" %in% colnames(relevantBioObservations)) { relevantBioObservations[, source := NULL] }
  # Return the filtered data.table
  return(relevantBioObservations)
}

getAverageEventDateByAgeGroup <- function(filteredPatients, firstCMDDates, sampledPatient) {
  nonCMDPatients <- data.table(setdiff(filteredPatients, firstCMDDates$patid))
  setnames(nonCMDPatients, "patid")
  sampledPatientClone <- copy(sampledPatient)
  # Define the start date of the study
  startDate <- as.Date("1990-01-01")
  studyStartYear <- as.numeric(format(startDate, "%Y"))
  
  # Calculate age at the start of the study
  sampledPatientClone[, age := studyStartYear - yob]
  
  # Merge with sampledPatient to get age, patId pair
  patIdsWithAges <- merge(nonCMDPatients, sampledPatientClone[, .(patid, age)], by = "patid", all.x = TRUE)
  patIdsWithAges <- convertAgeToCategory(patIdsWithAges)
  
  # Convert firstCMDDates$cmdDate to Date type if it's not already
  firstCMDDates[, date := as.Date(date)]
  
  # Step 1: Group the first CMD dates by age group and calculate the average date
  avgCMDDateByAgeGroup <- firstCMDDates[, .(avgCMDDate = mean(date, na.rm = TRUE)), by = age]
  
  # Step 2: Merge the average CMD date for each age group into the non-CMD patients
  patIdsWithAges <- merge(patIdsWithAges, avgCMDDateByAgeGroup, by = "age", all.x = TRUE)
  setorder(patIdsWithAges, patid)
  setnames(patIdsWithAges, "avgCMDDate", "date")
  # Now, patIdsWithAges has each non-CMD patient with their age group and the average CMD date for that age group
  return(patIdsWithAges)
}

getNonCMDBiomarkersPseudoEvent <- function(filteredPatients, firstCMDDates, sampledPatient) {
  patIdsWithExpectedDate <- getAverageEventDateByAgeGroup(filteredPatients, firstCMDDates, sampledPatient)
  return (getBiomarkersWithin5Years(patIdsWithExpectedDate))
}

getBaselineTable <- function(biomarkersCMD, biomarkersNonCMD, filteredPatients) {
  setkey(biomarkersCMD, patid, obsdate)  # Set a key for faster sorting and merging
  setkey(biomarkersNonCMD, patid, obsdate)  # Set a key for faster sorting and merging

  # Combine the unique patient IDs from both datasets
  uniqueCMD <- unique(biomarkersCMD$patid)
  uniqueNonCMD <- unique(biomarkersNonCMD$patid)
  allUniquePatients <- unique(c(uniqueCMD, uniqueNonCMD))  # Union of both patient ID sets
  message(glue("{length(allUniquePatients)} unique patients who have biomarkers"))
  # Filter both data tables to include all these patients
  #biomarkersCMD <- biomarkersCMD[patid %in% allUniquePatients]
  #biomarkersNonCMD <- biomarkersNonCMD[patid %in% allUniquePatients]

  # Define observation types with corresponding medcode IDs
  observationTypes <- list(
    bmi = medcodeBMI$medcodeid,
    totalcholesterol = medcodeTotalChol$medcodeid,
    dbp = medcodeDBP$medcodeid,
    sbp = medcodeSBP$medcodeid,
    fastingglucose = medcodeFastingGlucose$medcodeid,
    hba1c = medcodeHBA1C$medcodeid,
    hdl = medcodeHDL$medcodeid,
    ldl = medcodeLDL$medcodeid,
    triglyceride = medcodeTriglycerides$medcodeid,
    alcohol = medcodeAlcoholStatus$medcodeid,
    smoking = medcodeSmokingStatus$medcodeid,
    atrialfib = medcodeAtrialFib$medcodeid,
    hyperlipidaemia = medcodeHyperlipidaemia$medcodeid,
    hypertension = medcodeHypertension$medcodeid
  )
  # Process each observation type through extractObservations
  processedData <- lapply(names(observationTypes), function(type) {
    if (type %in% c("alcohol", "smoking", "atrialfib", "hyperlipidaemia", "hypertension")) {
      observations <- extractBooleanObservations(biomarkersCMD, biomarkersNonCMD, observationTypes[[type]], type)
    } else {
      observations <- extractObservations(biomarkersCMD, biomarkersNonCMD, observationTypes[[type]], type)
    }
    observations[, type := type]  # Add a column to keep track of the type
    return(observations)
  })

  alcoholObservations <- extractBooleanObservations(biomarkersCMD, biomarkersNonCMD, medcodeAlcoholStatus$medcodeid, "alcohol")
  smokingObservations <- extractBooleanObservations(biomarkersCMD, biomarkersNonCMD, medcodeSmokingStatus$medcodeid, "smoking")
  atrialFibObservations <- extractBooleanObservations(biomarkersCMD, biomarkersNonCMD, medcodeAtrialFib$medcodeid, "atrial fib")
  hyperlipidaemiaObservations <- extractBooleanObservations(biomarkersCMD, biomarkersNonCMD, medcodeHyperlipidaemia$medcodeid, "hyperlipidaemia")
  hypertensionObservations <- extractBooleanObservations(biomarkersCMD, biomarkersNonCMD, medcodeHypertension$medcodeid, "hypertension")
  medcodeAlcohol <- as.integer64(medcodeAlcoholStatus$medcodeid)
  medcodeAtrialFib <- as.integer64(medcodeAtrialFib$medcodeid)
  medcodeSmoking <- as.integer64(medcodeSmokingStatus$medcodeid)
  medcodeHypertension <- as.integer64(medcodeHypertension$medcodeid)
  medcodeHyperlipidaemia <- as.integer64(medcodeHyperlipidaemia$medcodeid)
  
  # Combine all processed data
  combinedObservations <- rbindlist(processedData, use.names = TRUE, fill = TRUE)
  setkey(combinedObservations, patid, obsdate)
  # Debug: Check the number of patients after combining observations
  message("Number of patients after combining observations: ", length(unique(combinedObservations$patid)))
  # Create a base data.table with all patients to ensure no one is left out
  allPatientData <- data.table(patid = allUniquePatients)
  
  # Full outer join to ensure all patient IDs are included
  combinedObservations <- allPatientData[combinedObservations, on = "patid", nomatch = NA]
  
  message(glue("{length(unique(combinedObservations$patid))} patients represented in the observations"))

  message("All data have been prepared. Please wait, this may take a while...")
  # Compute required results for each patient
  results <- combinedObservations[, .(
    age = as.numeric(format(min(obsdate), "%Y")) - sampledPatient[sampledPatient$patid == .BY$patid, yob],
    gender = sampledPatient[sampledPatient$patid == .BY$patid, gender],
    ethnicity = sampledEpiHes[sampledEpiHes$patid == .BY$patid, .SD[order(-admidate)][1, ethnos]],
    bmi = mean(value[type == "bmi"], na.rm = TRUE),
    cholesterol = mean(value[type == "totalcholesterol"], na.rm = TRUE),
    dbp = mean(value[type == "dbp"], na.rm = TRUE),
    sbp = mean(value[type == "sbp"], na.rm = TRUE),
    glucose = mean(value[type == "fastingglucose"], na.rm = TRUE),
    hba1c = mean(value[type == "hba1c"], na.rm = TRUE),
    hdl = mean(value[type == "hdl"], na.rm = TRUE),
    ldl = mean(value[type == "ldl"], na.rm = TRUE),
    triglycerides = mean(value[type == "triglyceride"], na.rm = TRUE),
    latestAlcoholUse = medcodeAlcoholStatus[medcodeid == alcoholObservations[patid == .BY$patid, .SD[order(-obsdate)]][1, medcodeid], category],
    latestSmokingStatus = medcodeSmokingStatus[medcodeid == smokingObservations[patid == .BY$patid, .SD[order(-obsdate)]][1, medcodeid], category],
    atrialFib = length(atrialFibObservations[patid == .BY$patid, patid]) > 0,
    hyperlipidaemia = length(hyperlipidaemiaObservations[patid == .BY$patid, patid]) > 0,
    hypertension = length(hypertensionObservations[patid == .BY$patid, patid]) > 0, 
    imd = sampledIMD2010[patid == .BY$patid, imd2010_5]
  ), by = .(patid)]

  # Ensure all results are combined into a single data.table
  output <- as.data.table(results)
  return(output)
}

categorizeColumn <- function(data, columnName, breaks, labels = NULL) {
  # Ensure column exists in data
  if (!columnName %in% names(data)) {
    stop(glue("Column '{columnName}' not found in the data."))
  }
  
  # Use cut to categorize the column based on the breaks and labels
  data[[columnName]] <- cut(data[[columnName]], 
                            breaks = breaks, 
                            labels = labels, 
                            include.lowest = TRUE, 
                            right = FALSE)
  
  # Return the updated data
  return(data)
}

performCategorization <- function(data) {
  # BMI
  bmiBreaks <- c(0, 18.5, 25, 30, Inf)
  bmiLabels <- c("<18.5", "18.5-25", "25-30", ">30")
  workingData <- categorizeColumn(data, "bmi", bmiBreaks, bmiLabels)
  
  # HDL
  hdlBreaks <- c(0, 1.03, 1.55, Inf)
  hdlLabels <- c("<1.03", "1.03-1.55", ">1.55")
  workingData <- categorizeColumn(workingData, "hdl", hdlBreaks, hdlLabels)

  # LDL
  ldlBreaks <- c(0, 2.6, 3.4, 4, 4.9, Inf)
  ldlLabels <- c("<2.6", "2.6-3.4", "3.4-4.0", "4.0-4.9", ">4.9")
  workingData <- categorizeColumn(workingData, "ldl", ldlBreaks, ldlLabels)
  
  # Triglycerides
  triglyBreaks <- c(0, 1.7, 2.3, 5.6, Inf)
  triglyLabels <- c("<1.7", "1.7-2.3", "2.3-5.6", ">5.6")
  workingData <- categorizeColumn(workingData, "triglycerides", triglyBreaks, triglyLabels)

  # Total Cholesterol
  totalCholBreaks <- c(0, 5, 5.5, 6, Inf)
  totalCholLabels <- c("<5.0", "5.0-5.5", "5.5-6.0", ">6.0")
  workingData <- categorizeColumn(workingData, "cholesterol", totalCholBreaks, totalCholLabels)

  # Glucose
  glucoseBreaks <- c(0, 5.5, 7, Inf)
  glucoseLabels <- c("<5.5", "5.5-7.0", ">7.0")
  workingData <- categorizeColumn(workingData, "glucose", glucoseBreaks, glucoseLabels)

  # Systolic Blood Pressure (SBP)
  sbpBreaks <- c(0, 120, 140, Inf)
  sbpLabels <- c("<120", "120-140", ">140")
  workingData <- categorizeColumn(workingData, "sbp", sbpBreaks, sbpLabels)

  # Diastolic Blood Pressure (DBP)
  dbpBreaks <- c(0, 80, 90, Inf)
  dbpLabels <- c("<80", "80-90", ">90")
  workingData <- categorizeColumn(workingData, "dbp", dbpBreaks, dbpLabels)

  # HbA1c
  hba1cBreaks <- c(0, 42, 48, Inf)
  hba1cLabels <- c("<42", "42-48", ">48")
  workingData <- categorizeColumn(workingData, "hba1c", hba1cBreaks, hba1cLabels)

  return(workingData)
}


###############################################################################



# Determines which patients have more than one year of history before their first
# (First year after adult age confirmed)
getCMDPatientsWith1yHistory <- function(fullCohort, firstCMDDates, cmdPatBiomarkers, sampledPatient) {
  # Ensure input data tables have keys set for faster subsetting and joining
  setkey(firstCMDDates, patid)
  setkey(cmdPatBiomarkers, patid)
  setkey(sampledPatient, patid)
  
  # Vectorized operation to calculate the cutoff date for each patient for CMD
  firstCMDDates[, cutoffDate := date %m-% months(12)]
  
  # Join firstCMDDates with sampledPatient to add 'yob' to firstCMDDates
  firstCMDDates <- firstCMDDates[sampledPatient, .(patid, date, cutoffDate, yob = i.yob), on = .(patid)]
  
  # Calculate the age at CMD event and ensure it's at least 18 years after 'yob'
  firstCMDDates[, ageAtCMD := year(date) - as.numeric(yob)]
  eligiblePatients <- firstCMDDates[ageAtCMD >= 18]
  
  # Calculate the date the patient turns 18
  eligiblePatients[, age18Date := as.IDate(paste(as.numeric(yob) + 18, "01", "01", sep="-"))]

  # Join eligiblePatients with biomarkers, adjusting filtering conditions
  eligiblePatients <- eligiblePatients[cmdPatBiomarkers, 
                                       .(patid, firstCMDDate = date, cutoffDate, obsdate, age18Date),
                                       on = .(patid),
                                       nomatch = 0]
  
  # Adjust filtering to check that there is any valid biomarker record before cutoff and after the patient turns 18
  eligiblePatients <- eligiblePatients[obsdate < cutoffDate & obsdate >= age18Date, .(patid), by = .(patid)]

  # Extract the 'patid' column as a vector, ensuring uniqueness
  uniqueEligiblePatids <- unique(eligiblePatients$patid)
  
  # Return the vector of unique patids
  return(uniqueEligiblePatids)
}


# Gets the first year of biomarker data for the non-CMD patients
# (First year after adult age confirmed)
getNonCMDBiomarkersWithinYear <- function(fullCohort, firstCMDDates) {
  combinedBiomarkerMedcodes <- unique(c(medcodeAtrialFib$medcodeid, 
                                        medcodeAlcoholStatus$medcodeid, 
                                        medcodeBMI$medcodeid,
                                        medcodeDBP$medcodeid,
                                        medcodeFastingGlucose$medcodeid,
                                        medcodeHBA1C$medcodeid,
                                        medcodeHDL$medcodeid,
                                        medcodeHyperlipidaemia$medcodeid,
                                        medcodeHypertension$medcodeid,
                                        medcodeLDL$medcodeid,
                                        medcodeSBP$medcodeid,
                                        medcodeSmokingStatus$medcodeid,
                                        medcodeTotalChol$medcodeid,
                                        medcodeTriglycerides$medcodeid))
  
  combinedBiomarkerMedcodes <- as.integer64(combinedBiomarkerMedcodes)
  
  nonCMDPatientIds <- setdiff(fullCohort$patid, firstCMDDates$patid)

  # Convert observationDataset to a data.table
  relevantBioObservations <- observationDataset %>%
                             filter(patid %in% nonCMDPatientIds, medcodeid %in% combinedBiomarkerMedcodes) %>%
                             collect() %>%
                             as.data.table()

  # Select only 'patid' and 'yob' from sampledPatient for joining
  sampledPatientDT <- as.data.table(sampledPatient)[, .(patid, yob)]

  # Join to include 'yob'
  relevantBioObservations <- merge(relevantBioObservations, sampledPatientDT, by = "patid")

  # Calculate 'yobPlus18'
  relevantBioObservations[, yobPlus18 := as.Date(paste0(yob, '-01-01')) + years(18)]

  # Filter for first observation after turning 18 and calculate 'firstObsDateAfter18'
  firstObservationsAfter18 <- relevantBioObservations[obsdate >= yobPlus18, .(firstObsDateAfter18 = min(obsdate)), by = patid]

  # Merge this information back to include 'firstObsDateAfter18' in the observations
  relevantBioObservations <- merge(relevantBioObservations, firstObservationsAfter18, by = "patid", all.x = TRUE)

  # Further filter to include only observations within 1 year of the first observation after turning 18
  finalOutput <- relevantBioObservations[obsdate <= firstObsDateAfter18 + years(1) & obsdate >= firstObsDateAfter18]

  return(finalOutput)
}
# Simple auxiliary function that produces a set of all biomarker medcodes in use
getCombinedBiomarkerMedcodes <- function() {
    combinedBiomarkerMedcodes <- unique(c(medcodeAtrialFib$medcodeid, 
                                        medcodeAlcoholStatus$medcodeid, 
                                        medcodeBMI$medcodeid,
                                        medcodeDBP$medcodeid,
                                        medcodeFastingGlucose$medcodeid,
                                        medcodeHBA1C$medcodeid,
                                        medcodeHDL$medcodeid,
                                        medcodeHyperlipidaemia$medcodeid,
                                        medcodeHypertension$medcodeid,
                                        medcodeLDL$medcodeid,
                                        medcodeSBP$medcodeid,
                                        medcodeSmokingStatus$medcodeid,
                                        medcodeTotalChol$medcodeid,
                                        medcodeTriglycerides$medcodeid))
    combinedBiomarkerMedcodes <- as.integer64(combinedBiomarkerMedcodes)
    return(combinedBiomarkerMedcodes)
}
# Gets the first year of biomarker data for CMD patients where the patient
# has a year or more of history before their first CMD event
# getCMDBiomarkersWithinYear <- function(fullCohort, firstCMDDates, patientsWith1yHistory) {
#   # Combine medcode IDs into one vector
#   combinedBiomarkerMedcodes <- unique(c(medcodeAtrialFib$medcodeid, 
#                                         medcodeAlcoholStatus$medcodeid, 
#                                         medcodeBMI$medcodeid,
#                                         medcodeDBP$medcodeid,
#                                         medcodeFastingGlucose$medcodeid,
#                                         medcodeHBA1C$medcodeid,
#                                         medcodeHDL$medcodeid,
#                                         medcodeHyperlipidaemia$medcodeid,
#                                         medcodeHypertension$medcodeid,
#                                         medcodeLDL$medcodeid,
#                                         medcodeSBP$medcodeid,
#                                         medcodeSmokingStatus$medcodeid,
#                                         medcodeTotalChol$medcodeid))
# 
#   # Filter observationDataset using dplyr and collect the results
#   relevantBioObservations <- observationDataset %>%
#                              filter(patid %in% patientsWith1yHistory, medcodeid %in% combinedBiomarkerMedcodes) %>%
#                              collect()
# 
#   # Convert to data.table
#   relevantBioObservations <- as.data.table(relevantBioObservations)
# 
#   # Find the first observation date for each patient
#   firstObservations <- relevantBioObservations[, .(firstObsDate = min(obsdate)), by = patid]
# 
#   # Merge this information back to the original dataset
#   relevantBioObservations <- merge(relevantBioObservations, firstObservations, by = "patid")
# 
#   # Filter to include only observations within 1 year of the first observation date
#   relevantBioObservations <- relevantBioObservations[obsdate <= firstObsDate + years(1) & obsdate >= firstObsDate]
# 
#   return(relevantBioObservations)
# }
# First year after age confirmed
getCMDBiomarkersWithinYear <- function(firstCMDDates, patientsWith1yHistory) {
  # Assume medcode vectors are defined somewhere in your environment
  combinedBiomarkerMedcodes <- unique(c(medcodeAtrialFib$medcodeid, 
                                        medcodeAlcoholStatus$medcodeid, 
                                        medcodeBMI$medcodeid,
                                        medcodeDBP$medcodeid,
                                        medcodeFastingGlucose$medcodeid,
                                        medcodeHBA1C$medcodeid,
                                        medcodeHDL$medcodeid,
                                        medcodeHyperlipidaemia$medcodeid,
                                        medcodeHypertension$medcodeid,
                                        medcodeLDL$medcodeid,
                                        medcodeSBP$medcodeid,
                                        medcodeSmokingStatus$medcodeid,
                                        medcodeTotalChol$medcodeid, 
                                        medcodeTriglycerides$medcodeid))
  
  combinedBiomarkerMedcodes <- as.integer64(combinedBiomarkerMedcodes)

  # Filter observationDataset for relevant biomarkers and patients with a history before their CMD event
  relevantBioObservations <- observationDataset %>%
                             filter(patid %in% patientsWith1yHistory, medcodeid %in% combinedBiomarkerMedcodes) %>%
                             collect()

  # Convert to data.table
  relevantBioObservations <- as.data.table(relevantBioObservations)

  # Ensure sampledPatient is a data.table for joining
  sampledPatientDT <- as.data.table(sampledPatient)
  
  # Join to include 'yob'
  relevantBioObservations <- merge(relevantBioObservations, sampledPatientDT[, .(patid, yob)], by = "patid")
  
  # Calculate when each patient turns 18
  relevantBioObservations[, yobPlus18 := as.Date(paste0(yob, '-01-01')) + years(18)]
  
  # Filter to find the first observation for each patient after they turn 18
  relevantBioObservations <- relevantBioObservations[obsdate >= yobPlus18]

  # Find the first observation date for each patient after they turn 18
  firstObservationsAfter18 <- relevantBioObservations[, .(firstObsDateAfter18 = min(obsdate)), by = patid]

  # Merge this information back to the original dataset
  relevantBioObservations <- merge(relevantBioObservations, firstObservationsAfter18, by = "patid")

  # Filter to include only observations within 1 year of the first observation date after turning 18
  relevantBioObservations <- relevantBioObservations[obsdate <= firstObsDateAfter18 + years(1) & obsdate >= firstObsDateAfter18]

  return(relevantBioObservations)
}
# Gets ALL biomarker data before each patient's first CMD Event
getCMDBiomarkerObservations <- function(firstCMDEventDates) {
  
  combinedBiomarkerMedcodes <- unique(c(medcodeAtrialFib$medcodeid, 
                                        medcodeAlcoholStatus$medcodeid, 
                                        medcodeBMI$medcodeid,
                                        medcodeDBP$medcodeid,
                                        medcodeFastingGlucose$medcodeid,
                                        medcodeHBA1C$medcodeid,
                                        medcodeHDL$medcodeid,
                                        medcodeHyperlipidaemia$medcodeid,
                                        medcodeHypertension$medcodeid,
                                        medcodeLDL$medcodeid,
                                        medcodeSBP$medcodeid,
                                        medcodeSmokingStatus$medcodeid,
                                        medcodeTotalChol$medcodeid,
                                        medcodeTriglycerides$medcodeid))
  
  combinedBiomarkerMedcodes <- as.integer64(combinedBiomarkerMedcodes)
  
  cmdDatePatients <- firstCMDEventDates[, .(patid)]
  
  filteredObservationDataset <- observationDataset %>%
    filter(patid %in% cmdDatePatients$patid) %>%
    filter(medcodeid %in% combinedBiomarkerMedcodes)
  
  relevantObservations <- as.data.table(filteredObservationDataset %>% collect())
  relevantObservations[, obsdate := as.Date(obsdate)]
  
  # Ensure date is in the correct format in firstCMDEventDates
  firstCMDEventDates[, date := as.Date(date)]
  
  # Convert both data.tables to ensure 'patid' is set as a key for a proper join
  setkey(firstCMDEventDates, patid)
  setkey(relevantObservations, patid)
  
  # Perform the join with an emphasis on matching 'patid' exactly
  # And then filter by 'obsdate' being before 'date' from firstCMDEventDates
  finalObservations <- firstCMDEventDates[relevantObservations, .(patid, consid, pracid, obsid, obsdate = i.obsdate, obsdate, staffid,
                                              parentobsid, medcodeid, value, numunitid, obstypeid, numrangelow, 
                                              numrangehigh, probobsid,
                                              date),
                                              on = .(patid), nomatch = 0][obsdate < date]
  
  # Remove temporary columns if needed, and return the final result
  finalObservations[, date := NULL] # Assuming 'date' is no longer needed
  
  return(finalObservations)
}

##################################
# Baseline table generation functions
##################################

extractBooleanObservations <- function(cmdPatientBiomarkers, nonCmdPatientBiomarkers, medcodeVar, type) {
    # Convert the specified medcode variable to integer64, assuming it's available in the environment
  medcodes <- as.integer64(medcodeVar)
  
  # Filter for observations based on the specified medcodes for both CMD and Non-CMD patients
  cmdObservations <- cmdPatientBiomarkers[medcodeid %in% medcodes]
  nonCmdObservations <- nonCmdPatientBiomarkers[medcodeid %in% medcodes]
  
  # Add a new column to indicate 'CMD' or 'Non-CMD'
  cmdObservations[, patientType := 'CMD']
  nonCmdObservations[, patientType := 'Non-CMD']
  
  # Combine the two data.tables
  combinedObservations <- rbind(cmdObservations, nonCmdObservations)
  
  message(glue("Number of patients who have this type of biomarker ({type}): {length(unique(combinedObservations[,patid]))}"))
  message(glue("Overall number of biomarker observations: {nrow(combinedObservations)}"))
  
  return(combinedObservations)
}
# Extracts the subset of observations matching a medcodeid vector parameter
extractObservations <- function(cmdPatientBiomarkers, nonCmdPatientBiomarkers, medcodeVar, type) {
  # Convert the specified medcode variable to integer64, assuming it's available in the environment
  medcodes <- as.integer64(medcodeVar)
  
  # Filter for observations based on the specified medcodes for both CMD and Non-CMD patients
  cmdObservations <- cmdPatientBiomarkers[medcodeid %in% medcodes]
  nonCmdObservations <- nonCmdPatientBiomarkers[medcodeid %in% medcodes]
  
  # Add a new column to indicate 'CMD' or 'Non-CMD'
  cmdObservations[, patientType := 'CMD']
  nonCmdObservations[, patientType := 'Non-CMD']
  
  # Combine the two data.tables
  combinedObservations <- rbind(cmdObservations, nonCmdObservations)
  
  message(glue("Number of patients who have this type of biomarker ({type}): {length(unique(combinedObservations[,patid]))}"))
  message(glue("Overall number of biomarker observations: {nrow(combinedObservations)}"))

  # Use EHRBiomarkr library to remove invalid values
  combinedObservations <- clean_biomarker_values(combinedObservations, value, type)
  
  message(glue("Number of patients who have this type of biomarker after cleaning: ({type}): {length(unique(combinedObservations[,patid]))}"))
  message(glue("Overall number of biomarker observations after cleaning: {nrow(combinedObservations)}"))

  
  return(combinedObservations)
}
# This function produces the baseline biomarker table.

# This function produces additional data complementing the baseline table, without augmenting the original table.
getSurvivalTable <- function(biomarkersCMD, biomarkersNonCMD) {
  allBiomarkers <- rbind(biomarkersCMD, biomarkersNonCMD)
  for(i in seq_len(baselineTable$patid)) {
    patientId <- baselineTable$patid[i]
    patientDied <- nrow(sampledDeath[patid == patientId]) > 0
    
    patientStrokeObservations <- allBiomarkers[patid == patientId & medcodeid %in% medcodeStroke$medcodeid]
    patientStrokeHosp <- sampledDiagEpi[patid == patientId & ICD %in% icdCodeStroke$icd]
    patientHadStroke <- nrow(patientStrokeObservations) > 0 | nrow(patientStrokeHosp) > 0
    
    patientMIObservations <- allBiomarkers[patid == patientId & medcodeid %in% medcodeMyoInf$medcodeid]
    patientMIHosp <- sampledDiagEpi[patid == patientId & ICD %in% icdCodeMyoInf$icd]
    patientHadMI <- nrow(patientMIObservations) > 0 | nrow(patientMIHosp) > 0
    
    patientDiabetesObservations <- allBiomarkers[patid == patientId & medcodeid %in% medcodeDiabetesT2$medcodeid]
    patientDiabetesHosp <- sampledDiagEpi[patid == patientId & ICD %in% icdCodeT2Diab$icd]
    patientHadDiabetes <- nrow(patientDiabetesObservations) > 0 | nrow(patientDiabetesHosp) > 0
    
    if (patientDied == TRUE) {
      deathDate <- sampledDeath[patid == patientId, dod]
      firstObservation <- allBiomarkers[patid == patientId, .(firstObsDate = min(obsdate))]
      interval <- firstObservation$firstObsDate %--% deathDate
      studyLength <- time_length(interval, unit = "month")
    }
    else {
      firstObservation <- allBiomarkers[patid == patientId, .(firstObsDate = min(obsdate))]
      interval <- firstObservation$firstObsDate %--% as.Date("2020-10-15", format="%Y-%m-%d")
      studyLength = time_length(interval, unit = "month")
    }
    
    # Output a data.table with columns 'patid, 'studyLength', 'stroke', 'mi', 'diabetes', 'isDead'
  }
}

getSurvivalTable_mine <- function(baselineTable, biomarkersCMD, biomarkersNonCMD) {
  patientIds <- baselineTable$patid
  message("Preparing data...")
  relevantObservations <- observationDataset |> filter(patid %in% patientIds & (medcodeid %in% medcodeStroke$medcodeid |
                                                                                medcodeid %in% medcodeMyoInf$medcodeid |
                                                                                medcodeid %in% medcodeDiabetesT2$medcodeid)) |> collect()
  setDT(relevantObservations)
  setkey(relevantObservations, patid)
  relevantSampledDeath <- sampledDeath[patid %in% patientIds]
  message("Data prepared. Analysing...")
  
  allBiomarkers <- rbind(biomarkersCMD, biomarkersNonCMD)
  
  survivalTable <- data.table(
  patid = as.double(patientIds),           # Ensure patid is of type double
  lengthOfStudy = as.double(rep(0, length(patientIds))),  # Initialize lengthOfStudy as double
  isDead = rep(FALSE, length(patientIds)), # Initialize isDead as boolean
  stroke = rep(FALSE, length(patientIds)), # Initialize stroke as boolean
  mi = rep(FALSE, length(patientIds)),     # Initialize mi as boolean
  diabetes = rep(FALSE, length(patientIds))# Initialize diabetes as boolean
)
  
  for(i in seq_len(nrow(survivalTable)))
  {
    patientId <- survivalTable[i]$patid
    isDead <- nrow(relevantSampledDeath[patid == patientId]) > 0
    survivalTable[i, isDead := isDead]
    
    firstObservation <- allBiomarkers[patid == patientId] |> arrange(obsdate) |> tail(1)
    if (isDead == TRUE) {
      deathDate <- sampledDeath[patid == patientId, dod]
      interval <- firstObservation$obsdate %--% deathDate
      studyLength <- time_length(interval, unit = "month")
    }
    else {
      interval <- firstObservation$obsdate %--% as.Date("2020-10-15", format="%Y-%m-%d")
      studyLength <- time_length(interval, unit = "month")
    }
    survivalTable[i, lengthOfStudy := studyLength]
    hadStroke <- nrow(relevantObservations[patid == patid & medcodeid %in% medcodeStroke$medcodeid]) > 0 | nrow(sampledDiagHosp[patid == patid & ICD %in% icdCodeStroke$icd]) > 0
    hadMI <- nrow(relevantObservations[patid == patid & medcodeid %in% medcodeMyoInf$medcodeid]) > 0 | nrow(sampledDiagHosp[patid == patid & ICD %in% icdCodeMyoInf$icd]) > 0
    hadDiabetes <- nrow(relevantObservations[patid == patid & medcodeid %in% medcodeDiabetesT2$medcodeid]) > 0 | nrow(sampledDiagHosp[patid == patid & ICD %in% icdCodeT2Diab$icd]) > 0
    survivalTable[i, stroke := hadStroke]
    survivalTable[i, mi := hadMI]
    survivalTable[i, mi := hadDiabetes]
    if (i %% 1 == 0)
    {
      message(glue("Processed {i} patients... please wait..."))
    }
  }

  return(survivalTable)
}

getSurvivalTable_mine2 <- function(baselineTable, biomarkersCMD, biomarkersNonCMD) {
  patientIds <- baselineTable$patid
  message("Preparing data...")
  
  # Pre-filter and collect relevant observations
  relevantObservations <- observationDataset %>%
    filter(patid %in% patientIds & 
           (medcodeid %in% c(medcodeStroke$medcodeid, 
                             medcodeMyoInf$medcodeid, 
                             medcodeDiabetesT2$medcodeid))) %>%
    collect()
  setDT(relevantObservations)
  setkey(relevantObservations, patid)

  # Preparing other necessary datasets
  relevantSampledDeath <- sampledDeath[patid %in% patientIds, .(patid, dod)]
  setDT(relevantSampledDeath)
  setkey(relevantSampledDeath, patid)

  message("Data prepared. Analysing...")
  
  # Combine and prepare biomarker data
  allBiomarkers <- rbind(biomarkersCMD, biomarkersNonCMD)
  setDT(allBiomarkers)
  setkey(allBiomarkers, patid)

  # Initialize the survival table
  survivalTable <- data.table(
    patid = as.double(patientIds),
    lengthOfStudy = 0,
    isDead = FALSE,
    stroke = FALSE,
    mi = FALSE,
    diabetes = FALSE
  )

  # Calculate first observation dates for each patient
  firstObservations <- allBiomarkers[, .(FirstObsDate = min(obsdate)), by = patid]
  setkey(firstObservations, patid)
  
  # Merge first observation dates and death information
  survivalTable <- merge(survivalTable, firstObservations, by = "patid")
  survivalTable <- merge(survivalTable, relevantSampledDeath, by = "patid", all.x = TRUE)

  # Calculate study length
survivalTable[, lengthOfStudy := fifelse(!is.na(dod),
                                         as.numeric(difftime(dod, FirstObsDate, units = "days")) / 30.44,
                                         as.numeric(difftime(as.Date("2020-10-15"), FirstObsDate, units = "days")) / 30.44)]
  survivalTable[, isDead := !is.na(dod)]

  # Calculate conditions using more controlled operations
  survivalTable[relevantObservations[medcodeid %in% medcodeStroke$medcodeid], stroke := TRUE, on = "patid"]
  survivalTable[relevantObservations[medcodeid %in% medcodeMyoInf$medcodeid], mi := TRUE, on = "patid"]
  survivalTable[relevantObservations[medcodeid %in% medcodeDiabetesT2$medcodeid], diabetes := TRUE, on = "patid"]

  message("Analysis complete.")
  return(survivalTable)
}



##################################
# Functions for meta-analysis and aggregation of baseline biomarker table
##################################
calculateAgeGenderDistribution <- function(baselineTable) {
  # Ensure gender is treated as a factor with levels in the desired order
  baselineTable$gender <- factor(baselineTable$gender, levels = c("M", "F", "I"))

  # Create age groups
  baselineTable <- baselineTable %>%
    mutate(age_group = case_when(
      age >= 18 & age <= 24 ~ '18-24',
      age >= 25 & age <= 34 ~ '25-34',
      age >= 35 & age <= 44 ~ '35-44',
      age >= 45 & age <= 54 ~ '45-54',
      age >= 55 & age <= 64 ~ '55-64',
      age >= 65 ~ '65+',
      TRUE ~ 'Unknown'
    ))
  
  # Calculate counts per group
  group_counts <- baselineTable %>%
    group_by(age_group, gender) %>%
    summarise(count = n(), .groups = 'drop')
  
  # Calculate the total population count for overall percentage calculation
  total_population <- nrow(baselineTable)
  
  # Calculate percentages based on the total population
  percentage_table <- group_counts %>%
    mutate(percentage = (count / total_population * 100)) %>%
    mutate(percentage = format(round(percentage, 2), nsmall = 2)) %>%
    arrange(age_group, gender)  # Ensure correct ordering here

  # Calculate standard deviation (across the proportions in each age group based on total population)
  std_dev_table <- group_counts %>%
    mutate(proportion = count / total_population) %>%
    group_by(age_group) %>%
    summarise(sd = round(sd(proportion), 2), .groups = 'drop') %>%
    arrange(age_group)  # Only age_group sorting needed here

  # Return the results as a list of data frames
  return(list(percentage_table = percentage_table, std_dev_table = std_dev_table))
}

calculateStatisticsByGroup <- function(data, columnName, breaks, genderLevels = c("M", "F", "I")) {
  # Ensure gender is treated as a factor with levels in the desired order
  data$gender <- factor(data$gender, levels = genderLevels)
  
  # Create a new grouping variable based on numerical groupings
  data <- data %>%
    mutate(group = cut(get(columnName), breaks = breaks, include.lowest = TRUE, right = FALSE))

  # Calculate total population count for percentage calculation
  total_population <- nrow(data)

  # Calculate mean, standard deviation, and count for each group by gender
  stats_table <- data %>%
    group_by(group, gender) %>%
    summarise(
      count = n(), 
      mean = sprintf("%.2f", mean(get(columnName), na.rm = TRUE)),  # Format mean to 2 decimal places as string
      sd = sd(get(columnName), na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(percentage = sprintf("%.2f", (count / total_population * 100))) %>%  # Format percentage to 2 decimal places as string
    select(group, gender, mean, sd, count, percentage)  # Specifying order of columns

  return(stats_table)
}

calculateValueMatches <- function(data, columnName, matchValues, genderLevels = c("M", "F", "I")) {
  # Ensure gender is treated as a factor with levels in the desired order
  data$gender <- factor(data$gender, levels = genderLevels)
  
  # Filter data to include only rows where the column matches one of the specified values
  filtered_data <- data %>%
    filter(get(columnName) %in% matchValues)

  # Calculate counts and percentages
  results <- filtered_data %>%
    group_by(Value = get(columnName), gender) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    mutate(Total = sum(Count)) %>%
    group_by(Value) %>%
    mutate(Percentage = sprintf("%.2f", round((Count / Total * 100), 2))) %>%
    ungroup() %>%
    select(Value, gender, Count, Percentage) %>%
    arrange(Value, gender)

  return(results)
}

showMultiDiseaseStats <- function(survivalTable) {
  message(glue("Number of heart attack-only patients: {nrow(survivalTable[stroke == FALSE & mi == TRUE & diabetes == FALSE])}"))
  message(glue("Number of heart attack patients in total: {nrow(survivalTable[mi == TRUE])}"))
  message("")
  message(glue("Number of stroke-only patients: {nrow(survivalTable[stroke == TRUE & mi == FALSE & diabetes == FALSE])}"))
  message(glue("Number of stroke patients in total: {nrow(survivalTable[stroke == TRUE])}"))
  message("")
  message(glue("Number of diabetes-only patients: {nrow(survivalTable[diabetes == TRUE & stroke == FALSE & mi == FALSE])}"))
  message(glue("Number of diabetes patients total: {nrow(survivalTable[diabetes == TRUE])}"))
  message("")
  message(glue("Number of patients with stroke + heart attack: {nrow(survivalTable[stroke == TRUE & mi == TRUE & diabetes == FALSE])}"))
  message(glue("Number of patients with stroke + diabetes: {nrow(survivalTable[stroke == TRUE & mi == FALSE & diabetes == TRUE])}"))
  message(glue("Number of patients with heart attack + diabetes: {nrow(survivalTable[stroke == FALSE & mi == TRUE & diabetes == TRUE])}"))
  message(glue("Number of patients with all three diseases: {nrow(survivalTable[stroke == TRUE & mi == TRUE & diabetes == TRUE])}"))
  message(glue("Number of patients without any measured CMD disease: {nrow(survivalTable[stroke == FALSE & mi == FALSE & diabetes == FALSE])}"))
}

getSurvivalTable <- function(baselineTable, biomarkersCMD, biomarkersNonCMD, firstCMDDates) {
  patientIds <- baselineTable$patid
  message("Preparing data...")

  # Use dplyr to filter observationDataset and collect it as a data frame
  relevantObservations <- observationDataset %>%
    filter(patid %in% patientIds,
           medcodeid %in% c(medcodeStroke$medcodeid, medcodeMyoInf$medcodeid, medcodeDiabetesT2$medcodeid)) %>%
    collect() %>%
    as.data.table()

  # Sampled death and hospital diagnosis processing
  relevantSampledDeath <- sampledDeath[patid %in% patientIds, .(patid, dod)]
  relevantDiagHosp <- sampledDiagHosp[patid %in% patientIds]

  # Merge biomarkers and compute the first date
  allBiomarkers <- rbind(biomarkersCMD, biomarkersNonCMD)
  firstBiomarkerDate <- allBiomarkers[patid %in% patientIds, .(first_date = min(obsdate)), by = patid]
  
  # Prepare health event indicators using both medcodeid and ICD
  healthEventsMedcode <- relevantObservations[, .(
    stroke = any(medcodeid %in% medcodeStroke$medcodeid),
    mi = any(medcodeid %in% medcodeMyoInf$medcodeid),
    diabetes = any(medcodeid %in% medcodeDiabetesT2$medcodeid)
  ), by = patid]
  
  healthEventsICD <- relevantDiagHosp[, .(
    stroke_ICD = any(ICD %in% icdCodeStroke$icd),
    mi_ICD = any(ICD %in% icdCodeMyoInf$icd),
    diabetes_ICD = any(ICD %in% icdCodeT2Diab$icd)
  ), by = patid]

  # Merge and resolve column names
  healthEvents <- merge(healthEventsMedcode, healthEventsICD, by = "patid", all = TRUE)
  healthEvents[, `:=` (
    stroke = coalesce(stroke, FALSE) | coalesce(stroke_ICD, FALSE),
    mi = coalesce(mi, FALSE) | coalesce(mi_ICD, FALSE),
    diabetes = coalesce(diabetes, FALSE) | coalesce(diabetes_ICD, FALSE)
  )]
  # Remove ICD columns before merging to avoid duplicates
  healthEvents[, c("stroke_ICD", "mi_ICD", "diabetes_ICD") := NULL]

  # Merge all data into survivalTable
  survivalTable <- data.table(patid = as.double(patientIds))
  survivalTable[, `:=` (isDead = FALSE, lengthOfStudy = 0)]
  survivalTable <- merge(survivalTable, firstBiomarkerDate, by = "patid", all.x = TRUE)
  survivalTable <- merge(survivalTable, healthEvents, by = "patid", all.x = TRUE)
  survivalTable <- merge(survivalTable, firstCMDDates, by = "patid", all.x = TRUE)
  survivalTable <- merge(survivalTable, relevantSampledDeath, by = "patid", all.x = TRUE)
  
  # Update table based on conditions
  survivalTable[, `:=` (
    isDead = !is.na(dod),
    lengthOfStudy = ifelse(isDead, 
                           as.numeric(interval(first_date, dod), 'days'),
                           as.numeric(interval(first_date, as.Date("2020-10-15")), 'days')),
    stroke = as.logical(stroke),
    mi = as.logical(mi),
    diabetes = as.logical(diabetes),
    timeToFirstCMDEvent = ifelse(stroke | mi | diabetes,
                                 as.numeric(interval(first_date, date), 'days'),
                                 NA_real_)
  )]

  message("Analysis complete.")

  return(survivalTable)
}

#write.csv(survivalTable, file="survivalTable.csv")
# There are negative values for timeToFirstCMDEvent. Why?


###############################################################################
# Revision section following expert advice

# The core idea here is to instead take the average of the 5 years leading
# up to each patient's CMD event

# Patients who never have a CMD event will not be considered in the biomarker
# set

###############################################################################
# Retainable functions:
# GetCMDBiomarkerObservations
# Get

###############################################################################
