################################################################################
# Cardiometabolic disease policy model in the UK setting
# Septiara Putri 
# Health Economics and Health Technology Assessment (HEHTA), 
# University of Glasgow, UK [2024]
################################################################################
# Biomarker descriptive statistic extractor
################################################################################
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

getRelevantObservations <- function(cmdDateData) {
  
  # Combine medical codes into one vector for efficient filtering 
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
  
  # Fetch the relevant observations
  relevantObservations <- observationDataset |>
                          filter(patid %in% cmdDateData$patid &
                                 medcodeid %in% combinedBiomarkerMedcodes) |>
                          select(patid, medcodeid, value, numunitid, obsdate) %>%
                          mutate(obsdate = as.Date(obsdate, format = "%Y-%m-%d"), source = "GP") %>%
                          collect()
  
  
  
}