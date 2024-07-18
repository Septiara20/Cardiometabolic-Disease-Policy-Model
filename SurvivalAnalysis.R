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
  combinedMedcodes <- unique(c(medcodeStroke$medcodeid, medcodeDiabetesT2$medcodeid, medcodeMyoInf$medcodeid))
  relevantObservations <- observationDataset %>%
    filter(patid %in% patIds & medcodeid %in% combinedMedcodes & obsdate >= ymd('1990-01-01')) %>%
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

generateStateTransitionTable <- function() {
  
  message("Getting the qualifying patients...")
  qualifyingPatients <- getQualifyingPatients()
  message("Getting the post-1990 CMD event data set...")
  hospitalDiagnoses <- getHospDiagnosesPost1990(qualifyingPatients)
  observations <- getCMDObservationsPost1990(qualifyingPatients)
  
  # Will optimize in the second pass...
  for(i in seq_len(length(qualifyingPatients)))
  {
    observationsForPatient <- observations[patid == qualifyingPatients[i]]
    hospitalDiagnoses <- hospitalDiagnoses[patid == qualifyingPatients[i]]
  }
  
  
  
}