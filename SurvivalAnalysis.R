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
    epikey = integer()
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
  
  message("Tranforming observation and hospital tables so that the columns match...")
  # Ensure date columns are of Date type and type columns are character
  observations[, date := as.Date(obsdate)]
  hospitalDiagnoses[, date := as.Date(epistart)]
  observations[, type := as.character(type)]
  hospitalDiagnoses[, type := as.character(type)]
  
  # Add a placeholder for the 'epikey' column in observations
  observations[, epikey := NA_integer_]
  
  # Select and rename columns to align
  observations <- observations[, .(patid, date, type, epikey)]
  hospitalDiagnoses <- hospitalDiagnoses[, .(patid, date, type, epikey)]
  
  message("Combining observations and hospital diagnoses...")
  # Merge observations and hospital diagnoses
  allEvents <- rbindlist(list(observations, hospitalDiagnoses), use.names = TRUE, fill = TRUE)
  
  message("Adding death events...")
  # Add death events
  sampledPatient[, date := as.Date(cprd_ddate)]
  sampledPatient[, type := "death"]
  sampledPatient[, epikey := NA_integer_] # death events don't have an epikey
  
  deathEvents <- sampledPatient[!is.na(date), .(patid, date, type, epikey)]
  
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

relabelFirstSecondCVD <- function(stateTransitionTable) {
  # Step 1: Retain the first and second entries for "stroke" and "mi"
  retainFirstSecond <- function(dt, eventType) {
    dt[type == eventType, eventCount := seq_len(.N), by = patid]
    firstSecondEntries <- dt[type == eventType & eventCount <= 2]
    firstSecondEntries[eventCount == 2, type := paste0("post-", eventType)]
    firstSecondEntries[, eventCount := NULL]
    return(firstSecondEntries)
  }
  
  # Handle "stroke" and "mi" events separately
  strokes <- retainFirstSecond(stateTransitionTable, "stroke")
  mis <- retainFirstSecond(stateTransitionTable, "mi")
  
  # Extract other events
  otherEvents <- stateTransitionTable[!type %in% c("stroke", "mi")]
  
  # Combine the results
  combinedEvents <- rbind(strokes, mis, otherEvents, fill = TRUE)
  
  # Sort the combined table
  setorder(combinedEvents, patid, date)
  
  combinedEvents[, eventCount := NULL]
  
  return(combinedEvents)
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
  
  return(stateTransitionTable)
}

cleanAndFilter <- function(stateTransitionTable) {
  message("Removing follow-on diabetes events...")
  stateTransitionTable <- filterFirstDiabetes(stateTransitionTable)
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
addExtraColumns <- function(dataTable, sampledPatient, sampledIMD2010) {
  # Define the start date of the study
  startDate <- as.Date("1990-01-01")
  studyStartYear <- as.numeric(format(startDate, "%Y"))
  
  # Calculate age at the start of the study
  sampledPatient[, age := studyStartYear - yob]
  
  # Map gender codes to full names
  genderMapping <- c("F" = "Female", "M" = "Male", "I" = "Intersex")
  sampledPatient[, gender := genderMapping[gender]]
  
  # Merge with sampledPatient to get age and gender
  dataTable <- merge(dataTable, sampledPatient[, .(patid, age, gender)], by = "patid", all.x = TRUE)
  
  # Merge with sampledIMD2010 to get deprivationIndex
  dataTable <- merge(dataTable, sampledIMD2010[, .(patid, imd2010_5)], by = "patid", all.x = TRUE)
  
  # Rename the deprivation index column
  setnames(dataTable, "imd2010_5", "deprivationIndex")
  
  return(dataTable)
}

# Gets all of the event data and turns it into wide format
performDataPreparation <- function (sampledPatient) {
  longStateTransitionTable <- generateStateTransitionTable(sampledPatient)
  longStateTransitionTable <- cleanAndFilter(longStateTransitionTable)
  longStateTransitionTable <- addBackPatientsWithNoEvents(longStateTransitionTable)
  wideTransitionTable <- convertToWideFormat(longStateTransitionTable)
  wideTransitionTable <- addExtraColumns(wideTransitionTable, sampledPatient, sampledIMD2010)
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
                         keep = c("age", "gender", "deprivationIndex"))
  return(preparedData)
}

# Can generalise for different covariates, but let's worry about that later
includeCovariates <- function(preparedMStateData, transitionMatrix) {
  covariates <- c("age", "gender", "deprivationIndex")
  
  preparedMStateData <- expand.covs(preparedMStateData, covariates, append = TRUE, longnames = TRUE)
  return(preparedMStateData)
}

fitWeightedCox <- function(MStateData) {
  coxfit <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = MStateData, method = "breslow")
  return(coxfit)
}

convertDaysToYears <- function(MSStateData) {
  MSStateData[, c("Tstart", "Tstop", "time")] <- MSStateData[, c("Tstart", "Tstop", "time")]/365.25
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

plotTest2 <- function(fittedCox, transitionMatrix) {
  resultsOfFit <- msfit(object = fittedCox, vartype = "greenwood", trans = transitionMatrix)
  probabilityTransform <- probtrans(resultsOfFit, predt = 0, method = "greenwood")
  
  state_labels <- c("Death", "Post-Stroke", "Stroke", "Post-MI", "MI", "Diabetes", "Healthy")
  statecols <- c("black", "purple", "blue", "maroon", "orange", "yellow", "lightgray")
  ord <- c(7, 6, 5, 4, 3, 2, 1)
  
  par(mfrow = c(4, 2))  # Create a 4x2 grid for individual plots
  for (i in 1:length(ord)) {
    plot(probabilityTransform$time, probabilityTransform$probs[, ord[i]], type = "l", col = statecols[ord[i]],
         xlab = "Years since study began", ylab = "Probability", main = state_labels[ord[i]], lwd = 2)
  }
  par(mfrow = c(1, 1))  # Reset the plotting grid to default
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


