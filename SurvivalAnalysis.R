################################################################################
# Cardiometabolic disease policy model in the UK setting
# Septiara Putri 
# Health Economics and Health Technology Assessment (HEHTA), 
# University of Glasgow, UK
# 2024
################################################################################
# Survival analysis and state transition table generator
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
library(mstate)
library(colorspace)
#install_github("Exeter-Diabetes/EHRBiomarkr")
library(EHRBiomarkr) 
library(cmprsk)

observationDataset <- open_dataset('Data/Observation/')

################################################################################
# Load other code here
source("Biomarkers.R")
source("TimeDependencies.R")
################################################################################

getQualifyingPatients <- function() {
  # We find patients who have a CMD diagnosis or event pre-1990 
  # and return the patients who qualify for the study
  combinedMedcodes <- unique(c(medcodeStroke$medcodeid, medcodeDiabetesT2$medcodeid, medcodeMyoInf$medcodeid))
  combinedICDCodes <- unique(c(icdCodeStroke$icd, icdCodeMyoInf$icd, icdCodeT2Diab$icd))
  
  patIds <- sampledPatient[,patid]
      
  pre90Observations <- observationDataset %>%
    filter(patid %in% patIds & medcodeid %in% combinedMedcodes & obsdate < ymd('1990-01-01')) %>%
    collect()
  pre90Observations <- setDT(pre90Observations)
  setkey(pre90Observations, "patid")
  

  pre90Diagnoses <- sampledDiagEpi[patid %in% patIds & ICD %in% combinedICDCodes & epistart < ymd('1990-01-01')]
  pre90Diagnoses <- setDT(pre90Diagnoses)
  disqualifiedPatientIds <- unique(union(pre90Diagnoses[,patid], pre90Observations[,patid]))
  qualifyingPatientIds <- setdiff(patIds, disqualifiedPatientIds)
  
  patientsOver18In1990 <- sampledPatient[yob < 1972, patid]
  message(glue("{nrow(sampledPatient) - length(patientsOver18In1990)} patients under 18 in 1990, {length(disqualifiedPatientIds)} patients with CMD pre-90, {nrow(sampledPatient)} total patients"))

  qualifyingPatientIds <- intersect(qualifyingPatientIds, patientsOver18In1990)
  message(glue("Total number of qualifying patients: {length(qualifyingPatientIds)}"))
  return(qualifyingPatientIds)
}

getCMDObservationsPost1990 <- function(patIds) {
    # Combine medical codes into one vector for efficient filtering
  relevantObservations <- observationDataset %>%
    filter(patid %in% patIds & medcodeid %in% medcodeDiabetesT2$medcodeid & obsdate >= ymd('1990-01-01')) %>%
    collect()
  relevantObservations <- setDT(relevantObservations)
  setkey(relevantObservations, "patid")
  return(relevantObservations)
}

getHospDiagnosesPost1990 <- function(patIds) {
    # Combine ICD codes into one vector for efficient filtering
  combinedICDCodes <- unique(c(icdCodeStroke$icd, icdCodeMyoInf$icd, icdCodeT2Diab$icd))
  return(sampledDiagEpi[patid %in% patIds & ICD %in% combinedICDCodes & epistart >= ymd('1990-01-01')])
}

generateStateTransitionTable <- function(sampledPatient) {
  results <- data.table(
    patid = integer(),
    date = as.Date(character()),
    type = character(),
    spno = integer(),
    source = character()  # Add a source column
  )
  
  message("Getting the qualifying patients...")
  qualifyingPatients <- getQualifyingPatients()
  message("Getting the post-1990 CMD event data set...")
  hospitalDiagnoses <- getHospDiagnosesPost1990(qualifyingPatients)
  observations <- getCMDObservationsPost1990(qualifyingPatients)
  sampledPatient <- sampledPatient[patid %in% qualifyingPatients]
  
  message("Converting medcodes to character type...")
  # Convert all IDs to character to avoid coercion issues
  medcodeStroke$medcodeid <- as.character(medcodeStroke$medcodeid)
  medcodeMyoInf$medcodeid <- as.character(medcodeMyoInf$medcodeid)
  medcodeDiabetesT2$medcodeid <- as.character(medcodeDiabetesT2$medcodeid)
  icdCodeStroke$icd <- as.character(icdCodeStroke$icd)
  icdCodeMyoInf$icd <- as.character(icdCodeMyoInf$icd)
  icdCodeT2Diab$icd <- as.character(icdCodeT2Diab$icd)
  observations$medcodeid <- as.character(observations$medcodeid)
  hospitalDiagnoses$ICD <- as.character(hospitalDiagnoses$ICD)
  
  message("Ensuring data.tables are built...")
  # Convert to data.table for faster processing
  setDT(hospitalDiagnoses)
  setDT(observations)
  setDT(sampledPatient)
  setDT(medcodeStroke)
  setDT(medcodeMyoInf)
  setDT(medcodeDiabetesT2)
  setDT(icdCodeStroke)
  setDT(icdCodeMyoInf)
  setDT(icdCodeT2Diab)
  
  # Determine the type for observations
  observations[, type := fcase(
    medcodeid %in% medcodeDiabetesT2$medcodeid, "diabetes",
    default = NA_character_
  )]
  
  # Determine the type for hospital diagnoses
  hospitalDiagnoses[, type := fcase(
    ICD %in% icdCodeStroke$icd, "stroke",
    ICD %in% icdCodeMyoInf$icd, "mi",
    ICD %in% icdCodeT2Diab$icd, "diabetes",
    default = NA_character_
  )]
  
  message("Transforming observation and hospital tables so that the columns match...")
  # Ensure date columns are of Date type and type columns are character
  observations[, date := as.Date(obsdate)]
  hospitalDiagnoses[, date := as.Date(epistart)]
  observations[, type := as.character(type)]
  hospitalDiagnoses[, type := as.character(type)]
  
  # Add a placeholder for the 'spno' column in observations
  observations[, spno := NA_integer_]
  
  # Add a source column
  observations[, source := "GP"]
  hospitalDiagnoses[, source := "Hospital"]
  
  # Select and rename columns to align
  observations <- observations[, .(patid, date, type, spno, source)]
  hospitalDiagnoses <- hospitalDiagnoses[, .(patid, date, type, spno, source)]
  
  message("Combining observations and hospital diagnoses...")
  # Merge observations and hospital diagnoses
  allEvents <- rbindlist(list(observations, hospitalDiagnoses), use.names = TRUE, fill = TRUE)
  
  message("Adding death events...")
  # Add death events with spno and source
  sampledPatient[, date := as.Date(cprd_ddate)]
  sampledPatient[, type := "death"]
  sampledPatient[, spno := NA_integer_] # death events don't have an spno
  sampledPatient[, source := "Registry"]
  
  deathEvents <- sampledPatient[!is.na(date), .(patid, date, type, spno, source)]
  
  message("Combining all events...")
  # Combine all events
  allEvents <- rbindlist(list(allEvents, deathEvents))
  
  message("Sorting...")
  # Sort the results by patid and date
  setorder(allEvents, patid, date)
  
  message("Augmenting with time information...")
  # Calculate 'day' and 'daysSinceLastEvent'
  allEvents[, day := as.integer(difftime(date, as.Date("1990-01-01"), units = "days"))]
  allEvents[, daysSinceLastEvent := c(0, diff(day)), by = patid]
  
  # Convert daysSinceLastEvent to numeric and correct the base date calculation
  allEvents[, daysSinceLastEvent := as.numeric(daysSinceLastEvent)]
  
  message("Done! Returning table...")
  return(allEvents)
}

removeEventsWithinContCare <- function(eventsTable, sampledDiagHosp) {
  # Convert 'discharged' to Date type if it's not already
  sampledDiagHosp[, discharged := as.Date(discharged)]
  
  # Precompute a mapping of spno to discharge date
  dischargeDateMap <- sampledDiagHosp[, .(discharged = max(discharged, na.rm = TRUE)), by = spno]
  
  # Add discharge dates to eventsTable using a join
  filteredEvents <- merge(eventsTable, dischargeDateMap, by = "spno", all.x = TRUE, suffixes = c("", "_discharged"))

  # Initialize a counter for tracking progress
  totalPatients <- length(unique(filteredEvents$patid))
  counter <- 0
  
  # Process each patid separately
  processedEvents <- filteredEvents[, {
    # Update progress counter and print progress
    counter <<- counter + 1
    if (counter %% 10000 == 0 || counter == totalPatients) {
      cat(sprintf("Processing patid %d of %d (%.2f%% complete)\n", counter, totalPatients, (counter / totalPatients) * 100))
    }
    
    # Logical vector to keep retained events
    retain <- rep(TRUE, .N)
    
    # Separate processing for 'mi' and 'stroke' events
    for (eventType in c("mi", "stroke")) {
      # Initialize the latest discharge date as a very early date
      latestDischargeDate <- as.Date("1990-01-01")
      
      # Track if the first event has been found
      firstEventFound <- FALSE
      
      # Iterate through events for this patid
      for (i in seq_len(.N)) {
        row <- .SD[i]
        
        # Skip events not of the specified type (mi or stroke)
        if (row$type != eventType) {
          next
        }
        
        # Process only the specified event type (mi or stroke)
        if (row$type == eventType) {
          # Retain the first occurrence
          if (!firstEventFound) {
            retain[i] <- TRUE
            firstEventFound <- TRUE
            # Only update latest discharge date if it's not NA
            if (!is.na(row$discharged)) {
              latestDischargeDate <- row$discharged
            }
          } else {
            # For subsequent events, check the discharge date
            if (!is.na(row$date) && !is.na(latestDischargeDate) && row$date <= latestDischargeDate) {
              retain[i] <- FALSE  # Discard this event
            } else {
              # Retain this event and update latest discharge date
              retain[i] <- TRUE
              # Only update if discharged is not NA
              if (!is.na(row$discharged)) {
                latestDischargeDate <- row$discharged
              }
            }
          }
        }
      }
    }
    
    # Return only events that are retained
    .SD[retain]
  }, by = patid]
  setorder(processedEvents, patid, date)
  # Rebuild daysSinceLastEvent
  processedEvents[, daysSinceLastEvent := c(0, diff(day)), by = patid]
  return(processedEvents)
}

relabelFirstSecondCVD <- function(longStateTransitionTable) {
  
  data <- copy(longStateTransitionTable)
  # Sort the data by patid and date to ensure events are processed in order
  setorder(data, patid, date)
  
  # Create a copy of the original type column to preserve the original labels
  data[, original_type := type]
  
  # Create a flag for the first and second occurrences of 'mi' and 'stroke'
  data[, event_number := 1:.N, by = .(patid, type)]
  
  # Relabel the second occurrence of 'mi' and 'stroke'
  data[event_number == 2 & type == "mi", type := "post-mi"]
  data[event_number == 2 & type == "stroke", type := "post-stroke"]
  
  # Filter to keep only the first two occurrences of 'mi' and 'stroke'
  filtered_data <- data[(event_number <= 2) | (type == "diabetes")]
  
  # Remove the helper columns used for processing
  filtered_data[, c("original_type", "event_number") := NULL]
  
  # Return the processed data
  return(filtered_data)
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

removeCompetingCVDEvents <- function(stateTransitionTable) {
  # Add an index column to keep track of the original order
  stateTransitionTable[, index := .I]
  
  # Identify the first cardiovascular event type and its index for each patient
  stateTransitionTable[, firstCVEvent := type[which(type %in% c("mi", "stroke"))[1]], by = patid]
  stateTransitionTable[, firstCVEventIndex := which(type %in% c("mi", "stroke"))[1], by = patid]
  
  # Handle patients without any cardiovascular events
  stateTransitionTable[is.na(firstCVEvent), firstCVEvent := "none"]
  
  # Filter the table to remove invalid transitions
  stateTransitionTable <- stateTransitionTable[, {
    if (firstCVEvent[1] == "mi") {
      .SD[!(type %in% c("stroke", "post-stroke") & .I > firstCVEventIndex[1])]
    } else if (firstCVEvent[1] == "stroke") {
      .SD[!(type %in% c("mi", "post-mi") & .I > firstCVEventIndex[1])]
    } else {
      .SD
    }
  }, by = patid]
  
  # Remove helper columns
  stateTransitionTable[, `:=`(firstCVEvent = NULL, firstCVEventIndex = NULL)]
  
  # Sort by the original index to maintain order
  setorder(stateTransitionTable, index)
  
  # Remove the index column
  stateTransitionTable[, index := NULL]
  
    # Rebuild daysSinceLastEvent
  stateTransitionTable[, daysSinceLastEvent := c(0, diff(day)), by = patid]
  return(stateTransitionTable)
}

filterFirstDiabetes <- function(stateTransitionTable) {
  # Filter rows with date on or before 2020-10-14
  message("Removing any events after 2020-10-14...")
  filteredTable <- stateTransitionTable[date <= as.Date("2020-10-14")]
  # Extract first diabetes diagnosis per patient
  firstDiabetes <-
    filteredTable[type == "diabetes", .SD[1], by = patid]
  # Extract all other diagnosis types
  otherDiagnoses <- filteredTable[type != "diabetes"]
  # Combine the filtered diabetes diagnoses with the other diagnoses
  finalFilteredTable <- rbind(firstDiabetes, otherDiagnoses)
  # Optional: to keep the rows sorted
  setorder(finalFilteredTable, patid, date)
  message("Rebuilding daysSinceLastEvent column...")
  finalFilteredTable[, daysSinceLastEvent := c(0, diff(day)), by = patid]
  return(finalFilteredTable)
}


cleanAndFilter <- function(stateTransitionTable) {
  message("Removing follow-on diabetes events...")
  stateTransitionTable <- filterFirstDiabetes(stateTransitionTable)
  message("Removing secondary events occurring within a period of continuous care...")
  stateTransitionTable <- removeEventsWithinContCare(stateTransitionTable, sampledDiagHosp)
  message("Retaining first and second CVD events, and renaming the second post-{CVD}...")
  stateTransitionTable <- relabelFirstSecondCVD(stateTransitionTable)
  message("Removing competing CVD events...")
  stateTransitionTable <- removeCompetingCVDEvents(stateTransitionTable)
  return(stateTransitionTable)
}

convertToWideFormat <- function(stateTransitionTable) {
  # Define the start and end dates
  startDate <- as.Date("1990-01-01")
  endDate <- as.Date("2020-10-14")
  
  # Calculate the length of the study
  studyLength <- as.numeric(difftime(endDate, startDate, units = "days"))
  
  # Initialize the wide format table
  wideTable <- unique(stateTransitionTable[, .(patid)])
  
  # Define the columns for the wide format
  eventTypes <- c("mi", "post-mi", "stroke", "post-stroke", "diabetes", "death")
  for (event in eventTypes) {
    wideTable[, (event) := NA_integer_]
    wideTable[, paste0(event, ".s") := 0]
  }
  
  # Fill in the columns based on the events in the long format data
  stateTransitionTable[, daysIntoStudy := day]
  for (event in eventTypes) {
    wideTable[stateTransitionTable[type == event], (event) := i.daysIntoStudy, on = .(patid)]
  }
  
  # Update the censoring columns for events that occurred
  for (event in eventTypes) {
    wideTable[!is.na(get(event)), paste0(event, ".s") := 1]
  }
  
  # Update censoring columns to the death time if the patient has died
  wideTable[!is.na(death), `:=`(
    mi = ifelse(is.na(mi), death, mi),
    `post-mi` = ifelse(is.na(`post-mi`), death, `post-mi`),
    stroke = ifelse(is.na(stroke), death, stroke),
    `post-stroke` = ifelse(is.na(`post-stroke`), death, `post-stroke`),
    diabetes = ifelse(is.na(diabetes), death, diabetes),
    mi.s = ifelse(is.na(mi), 0, mi.s),
    `post-mi.s` = ifelse(is.na(`post-mi`), 0, `post-mi.s`),
    stroke.s = ifelse(is.na(stroke), 0, stroke.s),
    `post-stroke.s` = ifelse(is.na(`post-stroke`), 0, `post-stroke.s`),
    diabetes.s = ifelse(is.na(diabetes), 0, diabetes.s)
  )]
  
  # Post-process to set the length of study for non-death events
  for (event in eventTypes) {
    wideTable[is.na(get(event)), (event) := studyLength]
  }
  
  return(wideTable)
}
# Function to add back patients with no events
addBackPatientsWithNoEvents <- function(stateTransitionTable) {
  qualifyingPatientIds <- getQualifyingPatients()
  
  # Identify patients with no events
  patientsWithNoEvents <- setdiff(qualifyingPatientIds, unique(stateTransitionTable$patid))
  message(glue("{length(patientsWithNoEvents)} patients identified with no events..."))
  message(glue("Total number of patients: {length(patientsWithNoEvents) + length(unique(stateTransitionTable$patid))}"))
  # Define the start and end dates
  startDate <- as.Date("1990-01-01")
  endDate <- as.Date("2020-10-14")
  
  # Calculate the length of the study
  studyLength <- as.numeric(difftime(endDate, startDate, units = "days"))
  
  # Create a data table for patients with no events
  noEventsTable <- data.table(patid = patientsWithNoEvents)
  eventTypes <- c("mi", "post-mi", "stroke", "post-stroke", "diabetes", "death")
  for (event in eventTypes) {
    noEventsTable[, (event) := NA_integer_]
    noEventsTable[, paste0(event, ".s") := 0]
  }
  
  # Set censoring at the end of the study for all event types
  for (event in eventTypes) {
    noEventsTable[, (event) := studyLength]
    noEventsTable[, paste0(event, ".s") := 0]
  }
  
  # Combine the no events data with the existing stateTransitionTable
  combinedTable <- rbind(stateTransitionTable, noEventsTable, fill = TRUE)
  
  setorder(combinedTable, patid)
  
  return(combinedTable)
}
# Function to add extra columns: age, gender, and deprivationIndex
addExtraColumns <- function (dataTable, sampledPatient, sampledIMD2010) {
  message("Adding covariate columns...")
  sampledPatientClone <- copy(sampledPatient)
  # Define the start date of the study
  startDate <- as.Date("1990-01-01")
  studyStartYear <- as.numeric(format(startDate, "%Y"))
  
  # Calculate age at the start of the study
  sampledPatientClone[, age := studyStartYear - yob]
  
  # Map gender codes to full names
  genderMapping <- c("F" = "Female", "M" = "Male", "I" = "Intersex")
  sampledPatientClone[, gender := genderMapping[gender]]
  
  # Merge with sampledPatient to get age and gender
  dataTable <- merge(dataTable, sampledPatientClone[, .(patid, age, gender)], by = "patid", all.x = TRUE)
  
  # Merge with sampledIMD2010 to get deprivationIndex
  dataTable <- merge(dataTable, sampledIMD2010[, .(patid, imd2010_5)], by = "patid", all.x = TRUE)
  
  # Rename the deprivation index column
  setnames(dataTable, "imd2010_5", "deprivationIndex")
  
  return(dataTable)
}

getBaselineBiomarkers <- function (sampledPatient) {
  filteredPatients <- getQualifyingPatients()
  firstCMDDates <- getFirstCMDDates(filteredPatients, sampledPatient)
  biomarkersCMD <- getBiomarkersWithin5Years(firstCMDDates)
  biomarkersNonCMD <- getNonCMDBiomarkersPseudoEvent(filteredPatients, firstCMDDates, sampledPatient)
  baselineTable <- getBaselineTable(biomarkersCMD, biomarkersNonCMD, filteredPatients)
  baselineTable <- performCategorization(baselineTable)
}

mergeDataFromBaseline <- function(wideTransitionTable, baselineTable) {
  # Merge all the categorized columns from baselineTable into wideTransitionTable
  mergedData <- merge(wideTransitionTable, 
                      baselineTable[, .(patid, bmi, hdl, ldl, triglycerides, cholesterol, glucose, sbp, dbp, hba1c, 
                                        latestSmokingStatus, atrialFib, hyperlipidaemia, hypertension)], 
                      by = "patid", all.x = TRUE)

  # Return the updated wideTransitionTable with the merged data
  
  # NOTE: This leaves the FH columns unpopulated for patients with no biomarkers!!
  return(mergedData)
}

addEthnicity <- function(wideTransitionTable, sampledEpiHes) {
  # For each patient (patid), get the most recent 'ethnos' value from 'sampledEpiHes'
  ethnicityData <- sampledEpiHes[, .(ethnicity = .SD[order(-admidate)][1, ethnos]), by = patid]

  # Merge the ethnicity data into wideTransitionTable using 'patid' as the key
  wideTransitionTable <- merge(wideTransitionTable, ethnicityData, by = "patid", all.x = TRUE)
  
  return(wideTransitionTable)
}


addAlcoholStatus <- function(wideTransitionTable) {
  # Define the end date of the study (based on earlier conversations, assuming it's 2019-12-31)
  studyEndDate <- as.Date("2019-12-31")
  
  # Extract patient IDs from wideTransitionTable
  patIds <- wideTransitionTable[, patid]
  
  # Filter alcohol observations: valid medcodeid for alcohol status and within the study period
  alcoholObservations <- observationDataset |> 
    filter(medcodeid %in% medcodeAlcoholStatus$medcodeid & 
           obsdate >= ymd('1990-01-01') & 
           obsdate <= studyEndDate & 
           patid %in% patIds) |> 
    collect()
  
  # Merge alcoholObservations with medcodeAlcoholStatus to get the 'category'
  alcoholObservations <- merge(alcoholObservations, medcodeAlcoholStatus, by = "medcodeid", all.x = TRUE)
  
  # Find the earliest alcohol observation category for each patient
  earliestAlcoholStatus <- alcoholObservations %>%
    group_by(patid) %>%
    summarize(alcoholStatus = first(category[order(obsdate)]), .groups = 'drop')
  
  # Merge the earliest alcohol status (category) back into the wideTransitionTable
  wideTransitionTable <- merge(wideTransitionTable, earliestAlcoholStatus, by = "patid", all.x = TRUE)
  
  return(wideTransitionTable)
}

addFamilyHistoryCovariates <- function(wideTransitionTable) {
  message("Adding family history covariate. May take a while...")
  familyHistoryObservations <- observationDataset |> filter(medcodeid %in% medcodeCVD_FH$medcodeid | medcodeid %in% medcodeDiab_FH$medcodeid) |> collect()
  setDT(familyHistoryObservations)
  
  cvdFamilyHistory <- familyHistoryObservations[medcodeid %in% medcodeCVD_FH$medcodeid]
  diabFamilyHistory <- familyHistoryObservations[medcodeid %in% medcodeDiab_FH$medcodeid]
  
    # Create diabetesFH and cvdFH columns with default FALSE
  wideTransitionTable$diabetesFH <- FALSE
  wideTransitionTable$cvdFH <- FALSE
  
  # Set to TRUE if patid is found in family history datasets
  wideTransitionTable$diabetesFH[wideTransitionTable$patid %in% diabFamilyHistory$patid] <- TRUE
  wideTransitionTable$cvdFH[wideTransitionTable$patid %in% cvdFamilyHistory$patid] <- TRUE
  return(wideTransitionTable)
}

# Gets all of the event data and turns it into wide format
performDataPreparation <- function (sampledPatient, sampledDiagHosp) {
  longStateTransitionTable <- generateStateTransitionTable(sampledPatient)
  longStateTransitionTable <- cleanAndFilter(longStateTransitionTable)
  longStateTransitionTable <- addBackPatientsWithNoEvents(longStateTransitionTable)
  wideTransitionTable <- convertToWideFormat(longStateTransitionTable)
  wideTransitionTable <- addExtraColumns(wideTransitionTable, sampledPatient, sampledIMD2010)
  wideTransitionTable <- addFamilyHistoryCovariates(wideTransitionTable)
  wideTransitionTable <- removeIntersexPatients(wideTransitionTable)
  return(wideTransitionTable)
}

removeIntersexPatients <- function(wideTransitionTable) {
  wideTransitionTable <- wideTransitionTable[gender != "Intersex"]
}

convertAgeToCategory <- function(wideTransitionTable) {
  wideTransitionTable$age<- cut(
    wideTransitionTable$age,
    breaks = c(-Inf, 24, 34, 44, 54, 64, Inf),
    labels = c("18-24", "25-34", "35-44", "45-54", "55-64", ">65"),
    right = TRUE  # Includes the upper bound in the interval
  )
  
  # Return the modified data frame
  return(wideTransitionTable)
}

createTransitionMatrix <- function() {
  tmat <-
    transMat(
      x = list(c(2, 3, 5, 7), c(3, 5, 7), c(4, 7), c(7), c(6, 7), c(7), c()),
      names = c(
        "Disease-free",
        "Diabetes",
        "MI - 1 event",
        "MI - >1 event",
        "Stroke - 1 event",
        "Stroke - >1 event",
        "Death"
      )
    )
}

prepareDataForModel <- function (transitionMatrix, wideTransitionTable, covariateList = NULL) {
  preparedData <- suppressWarnings(msprep(time = c(NA, "diabetes", "mi", "post-mi", "stroke", "post-stroke", "death"), 
                         trans = transitionMatrix,
                         data = wideTransitionTable,
                         status = c(NA, "diabetes.s", "mi.s", "post-mi.s", "stroke.s", "post-stroke.s", "death.s"),
                         keep = covariateList,
                         id = "patid"))
  #preparedData <- addTimeDependentAge(preparedData, wideTransitionTable)
  return(preparedData)
}

# # Can generalise for different covariates, but let's worry about that later
# includeCovariates <- function(preparedMStateData, wideTransitionTable) {
#   message("Expanding time-independent covariates gender, deprivationIndex, diabetesFH and cvdFH...")
#   covariates <- c("gender", "deprivationIndex", "diabetesFH", "cvdFH")
#   
#   preparedMStateData <- expand.covs(preparedMStateData, covariates, append = TRUE, longnames = TRUE)
#   
#   message("Adding time-dependent covariate age...")
#   preparedMStateData <- addTimeDependentAge(preparedMStateData, wideTransitionTable)
#   return(preparedMStateData)
# }


includeCovariates <- function(preparedMStateData, wideTransitionTable, covariates) {
  # Separate time-dependent and time-independent covariates
  timeDependentCovariates <- "age"  # Assuming 'age' is the only time-dependent one for now
  timeIndependentCovariates <- setdiff(covariates, timeDependentCovariates)

  allTimeIndependentCovs <- c("gender", "deprivationIndex", "diabetesFH", "cvdFH", 
                "bmi", "hdl", "ldl", "triglycerides", "cholesterol", 
                "glucose", "sbp", "dbp", "hba1c", "latestSmokingStatus", "alcoholStatus",
                "atrialFib", "hyperlipidaemia", "hypertension")
  
  # Expand time-independent covariates
  if (length(timeIndependentCovariates) > 0) {
    message(glue("Expanding time-independent covariates: {paste(timeIndependentCovariates, collapse = ', ')}..."))
    preparedMStateData <- expand.covs(preparedMStateData, timeIndependentCovariates, append = TRUE, longnames = TRUE)
  } else {
    # message("No time-independent covariates to expand.")
  }
  
  # Add time-dependent covariates (e.g., 'age')
  if (timeDependentCovariates %in% covariates) {
    message("Adding time-dependent covariate 'age'...")
    preparedMStateData <- addTimeDependentAge(preparedMStateData, wideTransitionTable)
  } else {
    # Handle cases where there are no time-dependent covariates
    message("No time-dependent covariates to add.")
  }

  # Step 3: Ensure that only expanded covariates are kept, and age if it exists in the covariate list
  expandedCovariates <- setdiff(names(preparedMStateData), covariates)
  
  # If age is in the covariate list, retain it, otherwise exclude it as well
  relevantCovariates <- if ("age" %in% covariates) {
    c(expandedCovariates, "age")
  } else {
    expandedCovariates
  }
  
  setDT(preparedMStateData)
  preparedMStateData <- preparedMStateData[, relevantCovariates, with = FALSE]
  
  return(preparedMStateData)
}


covariates <- c("gender", "deprivationIndex", "diabetesFH", "cvdFH", 
                "bmi", "hdl", "ldl", "triglycerides", "hba1c", "cholesterol", 
                "glucose", "sbp", "dbp", "latestSmokingStatus", 
                "atrialFib", "hyperlipidaemia", "hypertension")


runCoxModel <- function(wideTransitionTable, transitionMatrix, covariates = NULL, selectedTransitions = NULL) {
  
  wideTransitionTable <- replaceSpecialCharacters(wideTransitionTable)
  wideTransitionTable <- convertNAsToFalse(wideTransitionTable)
  
  wideTransitionTable <- addFactors(wideTransitionTable)
  # Step 1: Prepare data in long format using prepareDataForModel
  preparedData <- prepareDataForModel(transitionMatrix, wideTransitionTable, covariateList = covariates)
  
  # Step 2: Expand covariates and add time-dependent age
  preparedData <- includeCovariates(preparedData, wideTransitionTable, covariates)
  
  # Step 3: Filter out NAs for complete-case analysis
  filteredData <- na.omit(preparedData)  # Filters out rows with missing covariates
  
  print(glue("Number of patients in analysis: pre {length(unique(filteredData$patid))}"))
  setDT(filteredData)
  
  # Step 4: Check for selected transitions and filter if necessary
  if (!is.null(selectedTransitions)) {
    message(glue("Fitting model for transitions: {paste(selectedTransitions, collapse = ', ')}"))
    filteredData <- filteredData[trans %in% selectedTransitions]
  } else {
    selectedTransitions <- c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13')
    message("Fitting model for all transitions.")
  }
  
  gc() 
  rm(preparedData)
  gc()
  filteredData <- as.data.frame(filteredData)
  
  # Step 5: Build the formula dynamically based on whether 'age' is in the covariate list
  expandedCovariates <- colnames(filteredData)[grepl(paste0("\\.(", paste(selectedTransitions, collapse = "|"), ")$"), colnames(filteredData))]
  
  # Add age to the model conditionally if it is in the covariate list
  if ("age" %in% covariates) {
    formula <- as.formula(paste("Surv(Tstart, Tstop, status) ~", paste(c("age", expandedCovariates), collapse = " + ")))
  } else {
    formula <- as.formula(paste("Surv(Tstart, Tstop, status) ~", paste(expandedCovariates, collapse = " + ")))
  }
  
  # Step 6: Fit Cox model with the expanded covariates
  coxModel <- coxph(formula, data = filteredData, x = TRUE, method = "breslow")
  
  return(coxModel)
}

runModelForAllTransitionsAndAge <- function(wideTransitionTable, transitionNames) {
  doRunModelForTransitionsAndCov(wideTransitionTable, transitionMatrix, c("age"), transitionNames)
}

getSubsetWithAllCovs <- function(wideTransitionTable, covariates) {
  wideTransitionTable <- replaceSpecialCharacters(wideTransitionTable)
  wideTransitionTable <- convertNAsToFalse(wideTransitionTable)
  
  # Remove people with Intersex gender. Should be moved to cleaning function above
  wideTransitionTable <- wideTransitionTable[gender != "Intersex"]
  
  wideTransitionTable <- addFactors(wideTransitionTable)
  # Step 1: Prepare data in long format using prepareDataForModel
  preparedData <- prepareDataForModel(transitionMatrix, wideTransitionTable, covariateList = covariates)
  
  # Step 2: Expand covariates and add time-dependent age
  preparedData <- includeCovariates(preparedData, wideTransitionTable, covariates)
  
  # Step 3: Filter out NAs for complete-case analysis
  filteredData <- na.omit(preparedData)  # Filters out rows with missing covariates
  
  returnTable <- wideTransitionTable[patid %in% unique(filteredData$patid) ]
  return(returnTable)
}

# Runs the model for all covariates AT ONCE except hba1c
runModelForAllAtSameTime <- function(wideTransitionTable, transitionNames) {
  covariates <- c("gender", "deprivationIndex", "diabetesFH", "cvdFH", 
                "bmi", "hdl", "ldl", "triglycerides", "cholesterol", 
                "glucose", "sbp", "dbp", "latestSmokingStatus", 
                "atrialFib", "hyperlipidaemia", "hypertension", "age", "ethnicity")
  
    doRunModelForTransitionsAndCov(wideTransitionTable, transitionMatrix, covariates, transitionNames)
}

runModelForAllTransitionsAndCovs <- function(wideTransitionTable, transitionNames) {
  covariates <- c("hdl", "ldl", "sbp", "dbp","bmi", "triglycerides", "cholesterol", "glucose", "hba1c", "gender", "deprivationIndex", "diabetesFH", "cvdFH", "latestSmokingStatus", "atrialFib", "hyperlipidaemia", "hypertension", "alcoholStatus")
  
  for (covariate in covariates) {
    doRunModelForTransitionsAndCov(wideTransitionTable, transitionMatrix, covariate, transitionNames)
  }
}

doRunModelForTransitionsAndCov <- function(wideTransitionTable, transitionMatrix, covariate, transitionNames) {
  for (i in seq_along(transitionNames)) {
    print(runCoxModel(wideTransitionTable, transitionMatrix, covariate, c(i)))
    cat("\n")
  }
}

addFactors <- function(wideTransitionTable) 
{
  wideTransitionTable$gender <- factor(wideTransitionTable$gender, levels = c("Male", "Female"))
  wideTransitionTable$deprivationIndex <- factor(wideTransitionTable$deprivationIndex, levels = c(1, 2, 3, 4, 5))
  wideTransitionTable$latestSmokingStatus <- factor(wideTransitionTable$latestSmokingStatus, levels = c('Non smoker', 
                                                                                                        'Active smoker', 
                                                                                                        'Ex smoker'))
  wideTransitionTable$alcoholStatus <- factor(wideTransitionTable$alcoholStatus, levels = c('AlcoholConsumptionLevel0', 
                                                                                            'AlcoholConsumptionLevel1', 
                                                                                            'AlcoholConsumptionLevel2', 
                                                                                            'AlcoholConsumptionLevel3'))
  wideTransitionTable$bmi <- factor(wideTransitionTable$bmi, levels = c('less.than.18.5', '18.5.to.25', '25.to.30', 'greater.than.30'))
  
  wideTransitionTable$bmi <- relevel(wideTransitionTable$bmi, ref = '18.5.to.25')
  wideTransitionTable$ethnicity <- factor(wideTransitionTable$ethnicity, levels = c("White", "Indian", "Unknown", "Mixed", "Bl_Afric", "Other", 
                     "Chinese", "Bl_Other", "Bangladesi", "Pakistani", 
                     "Oth_Asian", "Bl_Carib"))
  wideTransitionTable$ethnicity <- relevel(wideTransitionTable$ethnicity, ref = 'White')
  wideTransitionTable$latestSmokingStatus <- relevel(wideTransitionTable$latestSmokingStatus, ref = 'Non smoker')
  
  return(wideTransitionTable)
}

replaceSpecialCharacters <- function(wideTransitionTable) {
  
  wideTransitionTable$latestSmokingStatus <- gsub("-", " ", wideTransitionTable$latestSmokingStatus)
  # Loop over all columns
  for (col in names(wideTransitionTable)) {
    # Check if the column is of type character or factor
    if (is.character(wideTransitionTable[[col]]) || is.factor(wideTransitionTable[[col]])) {
      # Replace '<' with 'lt' and '>' with 'gt'
      wideTransitionTable[[col]] <- gsub("<", "less.than.", wideTransitionTable[[col]])
      wideTransitionTable[[col]] <- gsub(">", "greater.than.", wideTransitionTable[[col]])
      wideTransitionTable[[col]] <- gsub("-", ".to.", wideTransitionTable[[col]])
    }
  }
  return(wideTransitionTable)
}

############################################################## Cox Analysis

# Use the mapping in your function to print the transition names
analyseLifestyleFactors <- function() {
  coxModelResults <- runCoxModel(wideTransitionTable, transitionMatrix, c('latestSmokingStatus', 'alcoholStatus', 'bmi'), c(1, 2, 4))
  
  # Loop through each transition and print the results with meaningful names
  for (transition in c(1, 2, 4)) {
    cat("Transition:", transitionNames[[as.character(transition)]], "\n")
    print(coxModelResults[[transition]])
    cat("\n")
  }
}

replaceCovariates <- function(wideTransitionTable, baselineTable) {
  # Define the covariates to replace
  covariates <- c("bmi", "hdl", "ldl", "triglycerides", "cholesterol", 
                  "glucose", "sbp", "dbp", "hba1c")
  
  # Extract the subset of wideTransitionTable that matches baselineTable
  matchedRows <- wideTransitionTable[patid %in% baselineTable$patid]
  
  # Perform a full join between matched rows from wideTransitionTable and baselineTable
  updatedTransitionTable <- merge(matchedRows, baselineTable[, c("patid", covariates), with = FALSE], 
                                  by = "patid", 
                                  all.x = TRUE, 
                                  suffixes = c("", "_baseline"))
  
  # Replace covariate values in wideTransitionTable only where non-missing in baselineTable
  for (cov in covariates) {
    updatedTransitionTable[[cov]] <- ifelse(!is.na(updatedTransitionTable[[paste0(cov, "_baseline")]]),
                                            updatedTransitionTable[[paste0(cov, "_baseline")]],
                                            updatedTransitionTable[[cov]])
    updatedTransitionTable[[paste0(cov, "_baseline")]] <- NULL  # Remove the temporary columns
  }
  
  # Find the unmatched rows from the original wideTransitionTable
  unmatchedRows <- wideTransitionTable[!(patid %in% baselineTable$patid)]
  
  # Combine the updated rows with the unmatched rows
  finalTransitionTable <- rbind(updatedTransitionTable, unmatchedRows)
  
  return(finalTransitionTable)
}

################################################################################
############################################################@@@@@@@ Utilities

convertNAsToFalse <- function(data) {
  cols_to_convert <- c("hypertension", "hyperlipidaemia", "atrialFib")
  
  for (col in cols_to_convert) {
    if (col %in% colnames(data)) {
      data[[col]][is.na(data[[col]])] <- FALSE
    }
  }
  
  return(data)
}

# Define a list to map transitions to names
transitionNames <- list(
  '1' = "Disease-free to Diabetes",
  '2' = "Disease-free to MI - 1 event",
  '3' = "Disease-free to Stroke - 1 event",
  '4' = "Disease-free to Death",
  '5' = "Diabetes to MI - 1 event",
  '6' = "Diabetes to Stroke - 1 event",
  '7' = "Diabetes to Death",
  '8' = "MI - 1 event to MI - >1 event",
  '9' = "MI - 1 event to Death",
  '10' = "MI - >1 event to Death",
  '11' = "Stroke - 1 event to Stroke - >1 event",
  '12' = "Stroke - 1 event to Death",
  '13' = "Stroke - >1 event to Death"
)



