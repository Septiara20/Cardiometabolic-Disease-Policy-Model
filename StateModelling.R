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

getFilteredDiagEpi <- function(patIds) {
  # Combine ICD codes into one vector for efficient filtering
  combinedICDCodes <- unique(c(icdCodeStroke$icd, icdCodeMyoInf$icd, icdCodeT2Diab$icd))
  return(sampledDiagEpi[patid %in% patIds & ICD %in% combinedICDCodes])
}

getFilteredObservations <- function(patIds) {
  # Combine medical codes into one vector for efficient filtering
  combinedMedcodes <- unique(c(medcodeStroke$medcodeid, medcodeDiabetesT2$medcodeid, medcodeMyoInf$medcodeid))
  # Pre-filter relevant observations and diagnosis, adding a source column
  relevantObservations <- observationDataset %>%
    filter(patid %in% patIds & medcodeid %in% combinedMedcodes) %>%
    collect()
  relevantObservations <- setDT(relevantObservations)
  setkey(relevantObservations, "patid")
  return(relevantObservations)
}

getDiabToCVDPats <- function(survivalTable, firstCMDDates) {
  output <- c()
  possiblePatIds <- survivalTable[diabetes == TRUE & (mi == TRUE | stroke == TRUE), patid]
  diagnoses <- getFilteredDiagEpi(possiblePatIds)
  observations <- getFilteredObservations(possiblePatIds)
  
  nonDiabMedcodes <- c(medcodeStroke$medcodeid, medcodeMyoInf$medcodeid)
  nonDiabICD <- c(icdCodeMyoInf$icd, icdCodeStroke$icd)

  for(possiblePatId in possiblePatIds) {
    patientObservations <- observations[patid == possiblePatId]
    patientDiagnoses <- observations[patid == possiblePatId]
    
    patientFirstCMD <- firstCMDDates[patid == possiblePatId]
    if (patientFirstCMD[, source] == 'GP') {
      observation <- observations[patid == possiblePatId & obsdate == patientFirstCMD[, date]]
      if (!any(observation$medcodeid %in% nonDiabMedcodes)) {
        append(output, possiblePatId)
      }
    } else if (patientFirstCMD[, source] == 'Hospital') {
      diagnosis <- diagnoses[patid == possiblePatId & epistart == patientFirstCMD[, date]]
      if (!any(diagnosis[, ICD] %in% nonDiabICD)) {
        append(output, possiblePatId)
      }
    }
  }
  return(output)
}

getDiabToCVDPats <- function(survivalTable, firstCMDDates) {
  # Filter possible patient IDs
  possiblePatIds <- survivalTable[diabetes == TRUE & (mi == TRUE | stroke == TRUE), patid]
  
  # Get filtered diagnoses and observations
  diagnoses <- getFilteredDiagEpi(possiblePatIds)
  observations <- getFilteredObservations(possiblePatIds)
  
  # Define non-diabetic medcodes and ICD codes
  nonDiabMedcodes <- c(medcodeStroke$medcodeid, medcodeMyoInf$medcodeid)
  nonDiabICD <- c(icdCodeMyoInf$icd, icdCodeStroke$icd)
  
  # Convert to data.table
  setDT(observations)
  setDT(diagnoses)
  setDT(firstCMDDates)
  
  # Convert datetime columns to Date type
  observations[, obsdate := as.IDate(obsdate)]
  diagnoses[, epistart := as.IDate(epistart)]
  firstCMDDates[, date := as.IDate(date)]
  
  # Filter firstCMDDates for relevant patients
  firstCMDDates <- firstCMDDates[patid %in% possiblePatIds]
  
  # Initialize output
  output <- c()
  
  # Create lookup tables for non-diabetic medcodes and ICD codes
  nonDiabMedcodeLookup <- unique(nonDiabMedcodes)
  nonDiabICDLookup <- unique(nonDiabICD)
  
  # Process GP source patients
  firstCMD_GP <- firstCMDDates[source == 'GP']
  if (nrow(firstCMD_GP) > 0) {
    merged_GP <- merge(observations, firstCMD_GP, by.x = c("patid", "obsdate"), by.y = c("patid", "date"))
    filtered_GP <- merged_GP[!(medcodeid %in% nonDiabMedcodeLookup), .(patid)]
    output <- unique(c(output, filtered_GP$patid))
  }
  
  # Process Hospital source patients
  firstCMD_Hospital <- firstCMDDates[source == 'Hospital']
  if (nrow(firstCMD_Hospital) > 0) {
    merged_Hospital <- merge(diagnoses, firstCMD_Hospital, by.x = c("patid", "epistart"), by.y = c("patid", "date"))
    filtered_Hospital <- merged_Hospital[!(ICD %in% nonDiabICDLookup), .(patid)]
    output <- unique(c(output, filtered_Hospital$patid))
  }
  
  return(output)
}

