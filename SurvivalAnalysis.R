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


library(data.table)

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
  #wideTransitionTable <- convertAgeToCategory(wideTransitionTable) # Will be made time-dependent later
  return(wideTransitionTable)
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

prepareDataForModel <- function (transitionMatrix, wideTransitionTable) {
  preparedData <- msprep(time = c(NA, "diabetes", "mi", "post-mi", "stroke", "post-stroke", "death"), 
                         trans = transitionMatrix,
                         data = wideTransitionTable,
                         status = c(NA, "diabetes.s", "mi.s", "post-mi.s", "stroke.s", "post-stroke.s", "death.s"),
                         keep = c("age", "gender", "deprivationIndex", "diabetesFH", "cvdFH"),
                         id = "patid")
  preparedData <- addTimeDependentAge(preparedData)
  return(preparedData)
}

addTimeDependentAge <- function(msdata) {
  # Calculate age at each Tstart (in years)
  msdata$age_at_Tstart <- msdata$age + (msdata$Tstart / 365.25)

  # Define the age categories based on the updated age
  msdata$age <- cut(msdata$age_at_Tstart, 
                                  breaks = c(-Inf, 24, 34, 44, 54, 64, Inf),
                                  labels = c("18-24", "25-34", "35-44", "45-54", "55-64", ">65"),
                                  right = FALSE)  # Adjust to include the lower bound correctly

  # Check for any NAs and handle them
  if (any(is.na(msdata$age_category_time))) {
    cat("Warning: NAs found in age_category_time after cutting. Consider reviewing the age ranges or input data.\n")
  }

  # Clean up intermediate columns if not needed
  msdata <- msdata[, !names(msdata) %in% c("age_at_Tstart")]

  return(msdata)
}


# Can generalise for different covariates, but let's worry about that later
includeCovariates <- function(preparedMStateData, transitionMatrix) {
  covariates <- c("age", "gender", "deprivationIndex", "diabetesFH", "cvdFH")
  
  preparedMStateData <- expand.covs(preparedMStateData, covariates, append = TRUE, longnames = TRUE)
  return(preparedMStateData)
}

################### FINE-GRAY ############################################

getCumInc <- function(preparedMStateData) {
  cif <- cuminc(ftime = preparedMStateData$time, fstatus = preparedMStateData$status)
  return(cif)
}

##########################################################################

# Below is Cox-related, prototype and initial testing only

fitWeightedCox <- function(MStateData) {
  coxfit <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = MStateData, method = "breslow")
  return(coxfit)
}

fitWeightedCox_WithCovariates <- function(MStateDate) {
                                        msprep(data = ebmt, trans = tmat, 
                                        time = c(NA, "rec", "ae","recae", "rel", "srv"), 
                                        status = c(NA, "rec.s", "ae.s", "recae.s","rel.s", "srv.s"), 
                                        keep = c("match", "proph", "year", "agecl"))}


convertDaysToYears <- function(MStateData) {
  MStateData[c("Tstart", "Tstop", "time")] <- MStateData[c("Tstart", "Tstop", "time")] / 365.25
  return(MStateData)
}

investigateDiabetesFHCovariate <- function(MStateData) {
  # SAutomatically include all expanded diabetesFH covariates
  covariates <- paste(paste0("diabetesFHTRUE.", 1:13), collapse = " + ")

  cox_model <- coxph(Surv(Tstart, Tstop, status) ~ diabetesFHTRUE.1, data = MStateData)

}

plotTest <- function(fittedCox, transitionMatrix) {
  resultsOfFit <- msfit(object = fittedCox, vartype = "greenwood", trans = transitionMatrix)
  probabilityTransform <- probtrans(resultsOfFit, predt = 0, method = "greenwood")
  #plot(resultsOfFit, las = 1, xlab = "Days since study began")
  state_labels <- c("Death", "Post-Stroke", "Stroke", "Post-MI", "MI", "Diabetes", "Healthy")

  statecols <- c("black", "purple", "blue", "maroon", "orange", "yellow", "lightgray")
  ord <- c(7, 6, 5, 4, 3, 2, 1)
  plot(probabilityTransform, ord = ord, xlab = "Years since study began", type = "filled", col = statecols[ord], legend.pos = "topleft", legend = c("", "", "", "", "", "", ""))
  legend("topleft", legend = state_labels, fill = statecols, cex = 0.8)
}

getProbabilitiesAtTime <- function(fittedCox, transitionMatrix, years) {
  library(mstate)
  
  # Fit the model
  resultsOfFit <- msfit(object = fittedCox, vartype = "greenwood", trans = transitionMatrix)
  probabilityTransform <- probtrans(resultsOfFit, predt = 0, method = "greenwood")
  
  # Convert years to days (assuming 365.25 days per year to account for leap years)
  target_time <- years * 365.25
  
  # Find the closest time point in the probabilityTransform
  closest_time_index <- which.min(abs(probabilityTransform[[1]]$time - target_time))
  closest_time <- probabilityTransform[[1]]$time[closest_time_index]
  
  # Get probabilities at the closest time point for each state
  state_probs <- sapply(probabilityTransform, function(x) x$prob[closest_time_index, ])
  
  # Debugging: Print the closest time point and probabilities
  print(paste("Closest time (days):", closest_time))
  print("Probabilities at closest time:")
  print(state_probs)
  
  # Create a table with the state labels and probabilities
  state_labels <- names(probabilityTransform)
  
  # Ensure probabilities vector matches the length of state_labels
  if (length(state_probs) != length(state_labels)) {
    stop("The number of probabilities does not match the number of states.")
  }
  
  prob_table <- data.frame(State = state_labels, Probability = state_probs)
  
  return(prob_table)
}


