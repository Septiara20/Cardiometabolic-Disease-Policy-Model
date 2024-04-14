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
install_github("Exeter-Diabetes/EHRBiomarkr")
library(EHRBiomarkr)
observationDataset <- open_dataset('Data/Observation/')
################################################################################
# Functions
################################################################################
# This function returns the balance of CMD and non-CMD patients in the full cohort
getCMDBalance <- function(fullCohort) {
  setDT(fullCohort)
  data <- data.table(
    numCMDPatients = numeric(),
    numNonCMDPatients = numeric()
  )

  message("Getting CMD patients from Observation...")
  cmdPatientsInObservation <- unique(observationDataset |> filter(medcodeid %in% medcodeDiabetesT2$medcodeid |
                                                         medcodeid %in% medcodeStroke$medcodeid |
                                                         medcodeid %in% medcodeMyocardialInf$medcodeid)
                                                        |> select(patid)
                                                        |> collect())
  
  message("Getting CMD patients from episodes...")
  cmdPatientsInHosp <- unique(sampledDiagEpi[ICD %in% icdCodeStroke$icd | ICD %in% icdCodeMyoInf$icd | ICD %in% icdCodeT2Diab$icd, patid])
  
  message("Combining...")
  combined <- unique(union(cmdPatientsInObservation$patid, cmdPatientsInHosp))
  
  message("Calculating inverse...")
  inverse <- fullCohort[!patid %in% combined]$patid
  
  newData <- data.table(numCMDPatients = length(combined), 
                        numNonCMDPatients = length(inverse))
  data <- rbind(data, newData)
  
  return(data)
}
# Probably unnecessary. Checking that the non-cmd patient value in above function
# is correct
getNonCMDPatients <- function(fullCohort) {
  setDT(fullCohort)
  
  message("Getting relevant obs...")
  relevantObservations <- data.table(observationDataset |> filter(patid %in% fullCohort$patid) |> select(patid, medcodeid) |> collect())
  relevantEpisodes <- sampledDiagEpi[patid %in% fullCohort$patid, .(patid, ICD)]
  
  message("Intersecting observations and episodes patients...")
  patIdsWithBoth <- intersect(unique(relevantObservations[,patid]), unique(relevantEpisodes[,patid]))
  nonCmd = numeric()

  
  message(glue("Processing {length(patIdsWithBoth)} patients... plooz waaaat"))
  for (i in seq_len(patIdsWithBoth)) {
    observationsForThisPatient <- relevantObservations[patid == patIdsWithBoth[i]]
    if (length(observationsForThisPatient[medcodeid %in% medcodeDiabetesT2 |
                                   medcodeid %in% medcodeMyocardialInf |
                                   medcodeid %in% medcodeStroke, medcodeid]) == 0) {
      nonCmd = nonCmd + 1
    }
    if (i %% 1 == 0)
    {
      message(glue("Processed {i} patients. Keep chil'n..."))
    }
  }
  return(nonCmd)
}
# Determines which patients have more than one year of history before their first
# CMD event
getCMDPatientsWith1yHistory <- function(fullCohort, firstCMDDates, cmdPatBiomarkers) {
  # Ensure input data tables have keys set for faster subsetting and joining
  setkey(firstCMDDates, patid)
  setkey(cmdPatBiomarkers, patid)
  
  # Vectorized operation to calculate the cutoff date for each patient
  # Note: This remains unchanged as it defines the period we're interested in based on the CMD diagnosis date
  firstCMDDates[, cutoffDate := date %m-% months(12)]
  
  # Adjust the join to use 'enterdate' from cmdPatBiomarkers instead of 'obsdate'
  # and compare it against the 'date' from firstCMDDates
  eligiblePatients <- firstCMDDates[cmdPatBiomarkers, 
                                    .(patid, firstCMDDate = date, cutoffDate, enterdate),
                                    on = .(patid, date > enterdate),
                                    nomatch = 0][
                                      enterdate < cutoffDate, .(patid), by = .(patid)]

  # The above returns a data.table grouped by patid, potentially with duplicates removed by the grouping operation
  # Extract the patid column as a vector
  uniqueEligiblePatids <- eligiblePatients$patid
  
  # Return the vector of unique patids
  return(uniqueEligiblePatids)
}
# First year after age confirmed
getCMDPatientsWith1yHistory_Over18 <- function(fullCohort, firstCMDDates, cmdPatBiomarkers, sampledPatient) {
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
  
  # Perform the join between eligiblePatients and cmdPatBiomarkers, filtering based on 'enterdate' and 'cutoffDate'
  eligiblePatients <- eligiblePatients[cmdPatBiomarkers, 
                                       .(patid, firstCMDDate = date, cutoffDate, enterdate, yob),
                                       on = .(patid),
                                       nomatch = 0][enterdate < cutoffDate & enterdate >= firstCMDDate - years(18), .(patid), by = .(patid)]
  
  # Extract the 'patid' column as a vector, ensuring uniqueness
  uniqueEligiblePatids <- unique(eligiblePatients$patid)
  
  # Return the vector of unique patids
  return(uniqueEligiblePatids)
}

# Gets the date of the first CMD event for patients who have a CMD event
getFirstCMDDates <- function(filteredPatients) {
  
  # Combine medical codes into one vector for efficient filtering
  combinedMedcodes <- unique(c(medcodeStroke$medcodeid, medcodeDiabetesT2$medcodeid, medcodeMyocardialInf$medcodeid))
  
  # Combine ICD codes into one vector for efficient filtering
  combinedICDCodes <- unique(c(icdCodeStroke$icd, icdCodeMyoInf$icd, icdCodeT2Diab$icd))
  
  # Pre-filter relevant observations and diagnosis, adding a source column
  relevantObservations <- observationDataset %>%
    filter(patid %in% filteredPatients$patid & medcodeid %in% combinedMedcodes) %>%
    select(patid, obsdate) %>%
    mutate(obsdate = as.Date(obsdate, format = "%Y-%m-%d"), source = "GP") %>%
    collect()

  relevantHospEpi <- sampledDiagEpi[patid %in% filteredPatients$patid & ICD %in% combinedICDCodes, .(patid, epistart = as.Date(epistart, format = "%Y-%m-%d"), source = "Hospital")]
  
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
  allPatients <- unique(c(filteredPatients$patid, relevantObservations$patid, relevantHospEpi$patid))
  allPatientsData <- data.table(patid = allPatients)
  setkey(allPatientsData, patid)
  setkey(eventDates, patid)
  
  # Join to ensure all patients are included, even those with all NA dates and determine the source
  finalEventDates <- allPatientsData[eventDates, on = "patid"]
  
  return(finalEventDates)
}

# Gets the first year of biomarker data for the non-CMD patients
getNonCMDBiomarkersWithinYear <- function(fullCohort, firstCMDDates) {
  # Combine medcode IDs into one vector
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
                                        medcodeTotalChol$medcodeid))

  nonCMDPatientIds <- setdiff(filteredPatients$patid, firstCMDDates$patid)

  # Filter observationDataset using dplyr and collect the results
  relevantBioObservations <- observationDataset %>%
                             filter(patid %in% nonCMDPatientIds, medcodeid %in% combinedBiomarkerMedcodes) %>%
                             collect()

  # Convert to data.table
  relevantBioObservations <- as.data.table(relevantBioObservations)

  # Find the first observation date for each patient
  firstObservations <- relevantBioObservations[, .(firstObsDate = min(enterdate)), by = patid]

  # Merge this information back to the original dataset
  relevantBioObservations <- merge(relevantBioObservations, firstObservations, by = "patid")

  # Filter to include only observations within 1 year of the first observation date
  relevantBioObservations <- relevantBioObservations[enterdate <= firstObsDate + years(1) & enterdate >= firstObsDate]

  return(relevantBioObservations)
}
# First year after age confirmed
getNonCMDBiomarkersWithinYear_Over18 <- function(fullCohort, firstCMDDates) {
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
  firstObservationsAfter18 <- relevantBioObservations[enterdate >= yobPlus18, .(firstObsDateAfter18 = min(enterdate)), by = patid]

  # Merge this information back to include 'firstObsDateAfter18' in the observations
  relevantBioObservations <- merge(relevantBioObservations, firstObservationsAfter18, by = "patid", all.x = TRUE)

  # Further filter to include only observations within 1 year of the first observation after turning 18
  finalOutput <- relevantBioObservations[enterdate <= firstObsDateAfter18 + years(1) & enterdate >= firstObsDateAfter18]

  return(finalOutput)
}

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
getCMDBiomarkersWithinYear <- function(fullCohort, firstCMDDates, patientsWith1yHistory) {
  # Combine medcode IDs into one vector
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
                                        medcodeTotalChol$medcodeid))

  # Filter observationDataset using dplyr and collect the results
  relevantBioObservations <- observationDataset %>%
                             filter(patid %in% patientsWith1yHistory, medcodeid %in% combinedBiomarkerMedcodes) %>%
                             collect()

  # Convert to data.table
  relevantBioObservations <- as.data.table(relevantBioObservations)

  # Find the first observation date for each patient
  firstObservations <- relevantBioObservations[, .(firstObsDate = min(enterdate)), by = patid]

  # Merge this information back to the original dataset
  relevantBioObservations <- merge(relevantBioObservations, firstObservations, by = "patid")

  # Filter to include only observations within 1 year of the first observation date
  relevantBioObservations <- relevantBioObservations[enterdate <= firstObsDate + years(1) & enterdate >= firstObsDate]

  return(relevantBioObservations)
}
# First year after age confirmed
getCMDBiomarkersWithinYear_Over18 <- function(firstCMDDates, patientsWith1yHistory) {
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
  relevantBioObservations <- relevantBioObservations[enterdate >= yobPlus18]

  # Find the first observation date for each patient after they turn 18
  firstObservationsAfter18 <- relevantBioObservations[, .(firstObsDateAfter18 = min(enterdate)), by = patid]

  # Merge this information back to the original dataset
  relevantBioObservations <- merge(relevantBioObservations, firstObservationsAfter18, by = "patid")

  # Filter to include only observations within 1 year of the first observation date after turning 18
  relevantBioObservations <- relevantBioObservations[enterdate <= firstObsDateAfter18 + years(1) & enterdate >= firstObsDateAfter18]

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
  finalObservations <- firstCMDEventDates[relevantObservations, .(patid, consid, pracid, obsid, obsdate = i.obsdate, enterdate, staffid,
                                              parentobsid, medcodeid, value, numunitid, obstypeid, numrangelow, 
                                              numrangehigh, probobsid,
                                              date),
                                              on = .(patid), nomatch = 0][obsdate < date]
  
  # Remove temporary columns if needed, and return the final result
  finalObservations[, date := NULL] # Assuming 'date' is no longer needed
  
  return(finalObservations)
}

getAtrialFibStats <- function (cmdPatientBiomarkers, nonCmdPatientBiomarkers) {
  medcodes <- as.integer64(medcodeAtrialFib$medcodeid)
  
  cmdPatients <- unique(cmdPatientBiomarkers[,patid])
  nonCmdPatients <- unique(nonCmdPatientBiomarkers[,patid])
  
  atrialFibPatients_CMD <- unique(cmdPatientBiomarkers[medcodeid %in% medcodes, patid])
  nonAtrialFibPatients_CMD <- setdiff(cmdPatients, atrialFibPatients_CMD)
  
  nonCmdPatients <- nonCmdPatientBiomarkers$patid
  atrialFibPatients_NonCMD <- unique(nonCmdPatientBiomarkers[medcodeid %in% medcodes, patid])
  nonAtrialFibPatients_NonCMD <- setdiff(nonCmdPatients, atrialFibPatients_NonCMD)
  
  stats = data.table(atrialFibCount_CMD = length(atrialFibPatients_CMD),
                     nonAtrialFibCount_CMD = length(nonAtrialFibPatients_CMD),
                     atrialFibCount_NonCMD = length(atrialFibPatients_NonCMD),
                     nonAtrialFibCount_NonCMD = length(nonAtrialFibPatients_NonCMD))
  
  return(stats)
}

getSmokingStats <- function(cmdPatientBiomarkers, nonCmdPatientBiomarkers) {
  setDT(medcodeSmokingStatus)
  
  cmdPatients <- unique(cmdPatientBiomarkers[,patid])
  nonCmdPatients <- unique(nonCmdPatientBiomarkers[,patid])
  
  currentSmokerMedcodes <- as.integer64(medcodeSmokingStatus[category == "Active smoker", medcodeid])
  pastSmokerMedcodes <- as.integer64(medcodeSmokingStatus[category == "Ex-smoker", medcodeid])
  nonSmokerMedcodes <- as.integer64(medcodeSmokingStatus[category == "Non-smoker", medcodeid])
  
  currentSmokerPatients_CMD <- unique(cmdPatientBiomarkers[medcodeid %in% currentSmokerMedcodes, patid])
  pastSmokerPatients_CMD <- unique(cmdPatientBiomarkers[medcodeid %in% pastSmokerMedcodes, patid])
  peopleWhoGaveUpSmoking_CMD <- intersect(currentSmokerPatients_CMD, pastSmokerPatients_CMD)
  nonSmokerPatients_CMD <- unique(cmdPatientBiomarkers[medcodeid %in% nonSmokerMedcodes, patid])
  # Create a vector of all patient IDs who have a smoking status recorded (both current and past smokers)
  smokerPatientIDs_CMD <- unique(c(currentSmokerPatients_CMD, pastSmokerPatients_CMD))

  # Identify nonSmokerPatients_CMD directly by checking which patients in the overall patient list (cmdPatients)
  # do not have their IDs in the smokerPatientIDs list. 
  # We achieve this by using the negation of `%in%` operator.
  nonSmokerPatients_CMD <- cmdPatients[!cmdPatients %in% smokerPatientIDs_CMD]
  
  ###############
  # Non-CMD
  ###############
  currentSmokerPatients_NonCMD <- unique(nonCmdPatientBiomarkers[medcodeid %in% currentSmokerMedcodes, patid])
  pastSmokerPatients_NonCMD <- unique(nonCmdPatientBiomarkers[medcodeid %in% pastSmokerMedcodes, patid])
  peopleWhoGaveUpSmoking_NonCMD <- intersect(currentSmokerPatients_NonCMD, pastSmokerPatients_NonCMD)
  nonSmokerPatients_NonCMD <- unique(nonCmdPatientBiomarkers[medcodeid %in% nonSmokerMedcodes, patid])
  # Create a vector of all patient IDs who have a smoking status recorded (both current and past smokers)
  smokerPatientIDs_NonCMD <- unique(c(currentSmokerPatients_NonCMD, pastSmokerPatients_NonCMD))

  # Identify nonSmokerPatients_CMD directly by checking which patients in the overall patient list (cmdPatients)
  # do not have their IDs in the smokerPatientIDs list. 
  # We achieve this by using the negation of `%in%` operator.
  nonSmokerPatients_NonCMD <- nonCmdPatients[!nonCmdPatients %in% smokerPatientIDs_NonCMD]
  
  stats = data.table(smokerPatients_CMD = length(currentSmokerPatients_CMD),
                     pastSmokerPatients_CMD = length(pastSmokerPatients_CMD),
                     nonSmokerPatients_CMD = length(nonSmokerPatients_CMD),
                     smokersWhoGaveUp_CMD = length(peopleWhoGaveUpSmoking_CMD),
                     nonSmokers_CMD = length(nonSmokerPatients_CMD),
                     smokerPatients_NonCMD = length(currentSmokerPatients_NonCMD),
                     pastSmokerPatients_NonCMD = length(pastSmokerPatients_NonCMD),
                     nonSmokerPatients_NonCMD = length(nonSmokerPatients_NonCMD),
                     smokersWhoGaveUp_NonCMD = length(peopleWhoGaveUpSmoking_NonCMD),
                     nonSmokers_NonCMD = length(nonSmokerPatients_NonCMD))
  return(stats)
}

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

getBaselineTable <- function(biomarkersCMD, biomarkersNonCMD) {
  setkey(biomarkersCMD, patid, obsdate)  # Set a key for faster sorting and merging
  setkey(biomarkersNonCMD, patid, obsdate)  # Set a key for faster sorting and merging

  # Combine the unique patient IDs from both datasets
  uniqueCMD <- unique(biomarkersCMD$patid)
  uniqueNonCMD <- unique(biomarkersNonCMD$patid)
  allUniquePatients <- unique(c(uniqueCMD, uniqueNonCMD))  # Union of both patient ID sets
  message(glue("{length(allUniquePatients)} unique patients"))
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
