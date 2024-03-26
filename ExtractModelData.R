################################################################################
# Cardiometabolic disease policy model in the UK setting
# Septiara Putri 
# Health Economics and Health Technology Assessment (HEHTA), 
# University of Glasgow, UK [2024]
################################################################################
# Yaldi ya bass wkwkwkwkwkwkwkwkwkwkwkkwkwkwkwkwkwkwkwkwkwkwkwkwkwkwkwkwkwkwkwkw
################################################################################
# Peptiara Pep Pep Tiara, Peptiara Pep Pep
################################################################################
getPostStrokeGPData <- function() {
  postStrokeData <- data.frame(patid = character(), 
                   firstObsDate = as.Date(character()), 
                   numDistinctDiagnoses = integer(), 
                   lastObsDate = as.Date(character()),
                   deathDate = as.Date(character())) 

  relevantObservations <- observationDataset |> 
      filter(patid %in% strokePats$patid & medcodeid %in% medcodeStroke$medcodeid) |>
      select(patid, obsdate, obsid, medcodeid) |>
      collect()
  
  relevantDeath <- sampledDeath |> filter(patid %in% strokePats$patid) |> select(patid, dod)
  
  for (i in 1:nrow(strokePats))
  {
    patientId <- strokePats$patid[i]
    allStrokeRowsForPatient <- relevantObservations |>
      filter(patid == patientId) |>
      arrange(obsdate) |>
      select(patid, obsdate, obsid, medcodeid)
      
    #message(glue("Patient {strokePats$patid[i]} had {nrow(allObsRowsForPatient)} entries in Observation"))
    
    #dateOfFirstObservation <- allObsRowsForPatient$obsdate[1]
    
    #message(glue("The first observation for this patient was {dateOfFirstObservation}"))
    #message(glue("The medcodeid of the first observation was {allObsRowsForPatient$medcodeid[1]} and the obsid was {allObsRowsForPatient$obsid[1]}"))
    
    #allStrokeRowsForPatient <- allObsRowsForPatient |>
    #  filter(medcodeid %in% medcodeStroke$medcodeid) |>
    #  arrange(obsdate, '%Y-%m-%d')
    
    #message(glue("This patient had {nrow(allStrokeRowsForPatient)} stroke-related observation(s)"))
            
    dateOfFirstStrokeObservation <- allStrokeRowsForPatient$obsdate[1]
    
    #message(glue("This patient's first stroke observation was on {dateOfFirstStrokeObservation}"))

    listOfObsDates <- unique(allStrokeRowsForPatient$obsdate)

    numDistinctDiagnoses <- length(listOfObsDates)  
    #message(glue("This patient had {numDistinctDiagnoses} distinct diagnoses"))
    
    deathDate <- relevantDeath |>
      filter(patid == patientId) |>
      select(dod)
    
    #message(glue("This patient died on {deathDate}"))
    lastObsDate <- allStrokeRowsForPatient$obsdate[length(allStrokeRowsForPatient$obsdate)]  
    #message(glue("This patient's last observation was on {lastObsDate}"))
    
    postStrokeData[nrow(postStrokeData) + 1,] <- list(patid = patientId, firstObsDate = dateOfFirstStrokeObservation,
                               numDistinctDiagnoses = numDistinctDiagnoses, lastObsDate = lastObsDate,
                               deathDate = as.Date(deathDate$dod))
    
    if (i %% 100 == 0)
    {
      message(glue("Processed {i} patients. Please wait..."))
    }
  }
  return (postStrokeData)
}

getPostStrokeEpiData <- function() {
    postStrokeEpiData <- data.frame(patid = character(), 
                                 firstEpiDate = as.Date(character()), 
                                 numDistinctDiagnoses = integer(), 
                                 lastEpiDate = as.Date(character()),
                                 deathDate = as.Date(character())) 
    
    relevantEpiData <- sampledDiagEpi |> 
      filter(ICD %in% icdCodeStroke$icd & patid %in% filteredPatients$patid) |>
      select(patid, epistart, epikey, spno, ICD)
    
    strokeEpiPats <- relevantEpiData |> select(patid) |> distinct()
    
    relevantDeath <- sampledDeath |> filter(patid %in% strokeEpiPats$patid) |> select(patid, dod)
    
    for (i in 1:nrow(strokeEpiPats))
    {
      patientId <- strokeEpiPats$patid[i]
      allStrokeRowsForPatient <- relevantEpiData |>
        filter(patid == patientId) |>
        arrange(epistart) |>
        select(patid, epistart, epikey, spno, ICD)
      
      
      dateOfFirstStrokeEpisode <- allStrokeRowsForPatient$epistart[1]
      
      #message(glue("This patient's first stroke observation was on {dateOfFirstStrokeObservation}"))
      
      listOfEpiKeys <- unique(allStrokeRowsForPatient$epikey)
      
      numDistinctDiagnoses <- length(listOfEpiKeys)  
      #message(glue("This patient had {numDistinctDiagnoses} distinct diagnoses"))
      
      deathDate <- relevantDeath |>
        filter(patid == patientId) |>
        select(dod)
      
      #message(glue("This patient died on {deathDate}"))
      lastEpiDate <- allStrokeRowsForPatient$epistart[length(allStrokeRowsForPatient$epistart)]  
      #message(glue("This patient's last observation was on {lastObsDate}"))
      
      postStrokeEpiData[nrow(postStrokeEpiData) + 1,] <- list(patid = patientId, firstEpiDate = dateOfFirstStrokeEpisode,
                                                        numDistinctDiagnoses = numDistinctDiagnoses, lastEpiDate = lastEpiDate,
                                                        deathDate = as.Date(deathDate$dod))
      
      if (i %% 1000 == 0)
      {
        message(glue("Processed {i} patients. Please wait..."))
      }
    }
    
    return (postStrokeEpiData)
}

# Returns a combined dataframe indexing number of hospital visits,
# GP visits and the other patient data
# Note: this is an inner join. Patients who only have
# entries in one of the two sets (hospital, GP) are excluded.
getCombinedStrokeData <- function() {
  postStrokeData <- data.frame(patid = double(), 
                               firstGPObsDate = as.Date(character()), 
                               numGPObservations = integer(), 
                               lastGPObsDate = as.Date(character()),
                               firstHospitalEpisode = as.Date(character()),
                               numHospitalisations = integer(),
                               numHospitalEpisodes = integer(),
                               lastHospitalEpisode = as.Date(character()),
                               deathDate = as.Date(character())) 
  
  relevantObservations <- data.table(observationDataset |> 
    filter(medcodeid %in% medcodeStroke$medcodeid & patid %in% filteredPatients$patid) |>
    select(patid, obsdate, obsid, medcodeid) |>
    collect())
  
  setKey(relevantObservations, patid)
  
  relevantEpiData <- data.table(sampledDiagEpi |> 
    filter(ICD %in% icdCodeStroke$icd & patid %in% filteredPatients$patid) |>
    select(patid, epistart, epikey, spno, ICD))
  
  setKey(relevantEpiData, patid)
  
  strokeEpiPats <- relevantEpiData |> select(patid) |> distinct()
  strokeHospPats <- relevantObservations |> select(patid) |> distinct()
  
  strokePatientIds <- union(strokeEpiPats, strokeHospPats, 'patid')
  
  relevantDeath <- sampledDeath |> filter(patid %in% strokePatientIds$patid) |> select(patid, dod)
  
  for (i in 1:nrow(strokePatientIds))
  {
    patientId <- strokePatientIds$patid[i]
    allGPRowsForPatient <- relevantObservations |>
      filter(patid == patientId) |>
      arrange(obsdate) |>
      select(patid, obsdate, obsid, medcodeid)
    
    allHospRowsForPatient <- relevantEpiData |>
      filter(patid == patientId) |>
      arrange(epistart) |>
      select(patid, epistart, epikey, spno, ICD)
    
    
    firstObsDate <- allGPRowsForPatient$obsdate[1]
    lastObsDate <- allGPRowsForPatient$obsdate[length(allGPRowsForPatient$obsdate)]  

    listOfObsDates <- unique(allGPRowsForPatient$obsdate)
    
    numDistinctObservations <- length(listOfObsDates)  
    
    firstEpiDate <- allHospRowsForPatient$epistart[1]
    lastEpiDate <- allHospRowsForPatient$epistart[length(allHospRowsForPatient$epistart)]  
    
    listOfEpiKeys <- unique(allHospRowsForPatient$epikey)
    listOfSpNo <- unique(allHospRowsForPatient$spno)

    numDistinctEpisodes <- length(listOfEpiKeys)
    numHospitalisations <- length(listOfSpNo)
    
    deathDate <- relevantDeath |>
      filter(patid == patientId) |>
      select(dod)
    
    #message(glue("This patient died on {deathDate}"))
    #message(glue("This patient's last observation was on {lastObsDate}"))
    
    postStrokeData[nrow(postStrokeData) + 1,] <- list(patid = patientId, 
                                                      firstGPObsDate = firstObsDate,
                                                      numGPObservations = numDistinctObservations, 
                                                      lastGPObsDate = lastObsDate,
                                                      firstHospitalEpisode = firstEpiDate, 
                                                      numHospitalisations = numHospitalisations,
                                                      numHospitalEpisodes = numDistinctEpisodes,
                                                      lastHospitalEpisode = lastEpiDate,
                                                      deathDate = as.Date(deathDate$dod))
    
    if (i %% 100 == 0)
    {
      message(glue("Processed {i} patients. Please wait..."))
    }
  }
  return (data.table(postStrokeData))
}

# Returns a combined dataframe indexing number of hospital visits,
# GP visits and the other patient data
# Note: this is an inner join. Patients who only have
# entries in one of the two sets (hospital, GP) are excluded.
getCombinedMyoInfData <- function() {
  postMIData <- data.frame(patid = double(), 
                               firstGPObsDate = as.Date(character()), 
                               numGPObservations = integer(), 
                               lastGPObsDate = as.Date(character()),
                               firstHospitalEpisode = as.Date(character()),
                               numHospitalisations = integer(),
                               numHospitalEpisodes = integer(),
                               lastHospitalEpisode = as.Date(character()),
                               deathDate = as.Date(character())) 
  
  relevantObservations <- data.table(observationDataset |> 
    filter(medcodeid %in% medcodeMyoInf$medcodeid & patid %in% filteredPatients$patid) |>
    select(patid, obsdate, obsid, medcodeid) |>
    collect())
  setkey(relevantObservations, patid)
  relevantEpiData <- data.table(sampledDiagEpi |> 
    filter(ICD %in% icdCodeMyoInf$icd & patid %in% filteredPatients$patid) |>
    select(patid, epistart, epikey, spno, ICD))
  setkey(relevantObservations, patid)
  
  strokeEpiPats <- relevantEpiData |> select(patid) |> distinct()
  strokeHospPats <- relevantObservations |> select(patid) |> distinct()
  
  strokePatientIds <- getIntersection(strokeEpiPats, strokeHospPats, 'patid')
  
  relevantDeath <- sampledDeath |> filter(patid %in% strokePatientIds$patid) |> select(patid, dod)
  
  for (i in 1:nrow(strokePatientIds))
  {
    patientId <- strokePatientIds$patid[i]
    allGPRowsForPatient <- relevantObservations |>
      filter(patid == patientId & !is.na(obsdate)) |>
      arrange(obsdate) |>
      select(patid, obsdate, obsid, medcodeid)
    
    allHospRowsForPatient <- relevantEpiData |>
      filter(patid == patientId) |>
      arrange(epistart) |>
      select(patid, epistart, epikey, spno, ICD)
    
    
    firstObsDate <- allGPRowsForPatient$obsdate[1]
    lastObsDate <- allGPRowsForPatient$obsdate[length(allGPRowsForPatient$obsdate)]  
    
    listOfObsDates <- unique(allGPRowsForPatient$obsdate)
    
    numDistinctObservations <- length(listOfObsDates)  
    
    firstEpiDate <- allHospRowsForPatient$epistart[1]
    lastEpiDate <- allHospRowsForPatient$epistart[length(allHospRowsForPatient$epistart)]  
    
    listOfEpiKeys <- unique(allHospRowsForPatient$epikey)
    listOfSpNo <- unique(allHospRowsForPatient$spno)
    
    numDistinctEpisodes <- length(listOfEpiKeys)
    numHospitalisations <- length(listOfSpNo)
    
    deathDate <- relevantDeath |>
      filter(patid == patientId) |>
      select(dod)
    
    #message(glue("This patient died on {deathDate}"))
    #message(glue("This patient's last observation was on {lastObsDate}"))
    
    postMIData[nrow(postMIData) + 1,] <- list(patid = patientId, 
                                                      firstGPObsDate = firstObsDate,
                                                      numGPObservations = numDistinctObservations, 
                                                      lastGPObsDate = lastObsDate,
                                                      firstHospitalEpisode = firstEpiDate, 
                                                      numHospitalisations = numHospitalisations,
                                                      numHospitalEpisodes = numDistinctEpisodes,
                                                      lastHospitalEpisode = lastEpiDate,
                                                      deathDate = as.Date(deathDate$dod))
    
    if (i %% 100 == 0)
    {
      message(glue("Processed {i} patients. Please wait..."))
    }
  }
  return (data.table(postMIData))
}

# Returns a combined dataframe indexing number of hospital visits,
# GP visits and the other patient data
# Note: this is an inner join. Patients who only have
# entries in one of the two sets (hospital, GP) are excluded.
getCombinedT2DData <- function() {
  postT2Data <- data.frame(patid = double(), 
                           firstGPObsDate = as.Date(character()), 
                           numGPObservations = integer(), 
                           lastGPObsDate = as.Date(character()),
                           firstHospitalEpisode = as.Date(character()),
                           numHospitalisations = integer(),
                           numHospitalEpisodes = integer(),
                           lastHospitalEpisode = as.Date(character()),
                           deathDate = as.Date(character())) 
  
  message("Preparing data...")
  relevantObservations <- data.table(observationDataset |> 
                                       filter(medcodeid %in% medcodeDiabetesT2$medcodeid & patid %in% filteredPatients$patid) |>
                                       select(patid, obsdate, obsid, medcodeid) |>
                                       collect())
  setkey(relevantObservations, patid)
  relevantEpiData <- data.table(sampledDiagEpi |> 
                                  filter(ICD %in% icdCodeT2Diab$icd & patid %in% filteredPatients$patid) |>
                                  select(patid, epistart, epikey, spno, ICD))
  setkey(relevantEpiData, patid)
  
  strokeEpiPats <- relevantEpiData |> select(patid) |> distinct()
  strokeHospPats <- relevantObservations |> select(patid) |> distinct()
  
  strokePatientIds <- getIntersection(strokeEpiPats, strokeHospPats, 'patid')
  
  relevantDeath <- sampledDeath |> filter(patid %in% strokePatientIds$patid) |> select(patid, dod)
  
  message("Data prepared! Please wait...")
  
  for (i in 1:nrow(strokePatientIds))
  {
    patientId <- as.double(strokePatientIds$patid[i])
    allGPRowsForPatient <- relevantObservations[patid == patientId & !is.na(obsdate)] |>
      #filter(patid == patientId & !is.na(obsdate)) |>
      arrange(obsdate) |>
      select(patid, obsdate, obsid, medcodeid)
    
    allHospRowsForPatient <- relevantEpiData[patid == patientId] |>
      #filter(patid == patientId) |>
      arrange(epistart) |>
      select(patid, epistart, epikey, spno, ICD)
    
    
    firstObsDate <- allGPRowsForPatient$obsdate[1]
    lastObsDate <- allGPRowsForPatient$obsdate[length(allGPRowsForPatient$obsdate)]  
    
    listOfObsDates <- unique(allGPRowsForPatient$obsdate)
    
    numDistinctObservations <- length(listOfObsDates)  
    
    firstEpiDate <- allHospRowsForPatient$epistart[1]
    lastEpiDate <- allHospRowsForPatient$epistart[length(allHospRowsForPatient$epistart)]  
    
    listOfEpiKeys <- unique(allHospRowsForPatient$epikey)
    listOfSpNo <- unique(allHospRowsForPatient$spno)
    
    numDistinctEpisodes <- length(listOfEpiKeys)
    numHospitalisations <- length(listOfSpNo)
    
    deathDate <- relevantDeath |>
      filter(patid == patientId) |>
      select(dod)
    
    #message(glue("This patient died on {deathDate}"))
    #message(glue("This patient's last observation was on {lastObsDate}"))
    
    postT2Data[nrow(postT2Data) + 1,] <- list(patid = patientId, 
                                              firstGPObsDate = firstObsDate,
                                              numGPObservations = numDistinctObservations, 
                                              lastGPObsDate = lastObsDate,
                                              firstHospitalEpisode = firstEpiDate, 
                                              numHospitalisations = numHospitalisations,
                                              numHospitalEpisodes = numDistinctEpisodes,
                                              lastHospitalEpisode = lastEpiDate,
                                              deathDate = as.Date(deathDate$dod))
    
    if (i %% 100 == 0)
    {
      message(glue("Processed {i} patients. Please wait..."))
    }
  }
  return (data.table(postT2Data))
}


getCombinedT2DData_Optimized <- function()
{
  # Assuming 'observationDataset', 'sampledDiagEpi', and 'sampledDeath' are already loaded as data.table
  setDT(sampledDeath)
  setDT(sampledDiagEpi)
  setkey(sampledDiagEpi, patid)
  setkey(sampledDeath, patid)
  
  message("Preparing data...")
  
  # Filter observations and episodes only once
  relevantObservations <- data.table(observationDataset |> 
                                       filter(medcodeid %in% medcodeDiabetesT2$medcodeid & patid %in% filteredPatients$patid) |>
                                       select(patid, obsdate, obsid, medcodeid) |>
                                       collect())
  relevantEpiData <-
    sampledDiagEpi[ICD %in% icdCodeT2Diab$icd &
                     patid %in% filteredPatients$patid, .(patid, epistart, epikey, spno, ICD)]
  
  # Compute unique diabetes patients based on the filtered data
  strokePatientIds <-
    unique(c(relevantObservations$patid, relevantEpiData$patid))
  
  relevantDeath <-
    sampledDeath[patid %in% strokePatientIds, .(patid, dod)]
  
  # Pre-allocate the output data frame with the correct size
  numPatients <- length(strokePatientIds)
  postT2Data <- data.table(
    patid = numeric(numPatients),
    firstGPObsDate = as.Date(rep(NA, numPatients)),
    numGPObservations = integer(numPatients),
    lastGPObsDate = as.Date(rep(NA, numPatients)),
    firstHospitalEpisode = as.Date(rep(NA, numPatients)),
    numHospitalisations = integer(numPatients),
    numHospitalEpisodes = integer(numPatients),
    lastHospitalEpisode = as.Date(rep(NA, numPatients)),
    deathDate = as.Date(rep(NA, numPatients))
  )
  
  message("Data prepared! Please wait...")
  
  for (i in seq_len(numPatients)) {
    patientId <- strokePatientIds[i]
    allGPRowsForPatient <-
      relevantObservations[patid == patientId, .SD[order(obsdate)], .SDcols = c("obsdate", "obsid", "medcodeid")]
    allHospRowsForPatient <-
      relevantEpiData[patid == patientId, .SD[order(epistart)], .SDcols = c("epistart", "epikey", "spno", "ICD")]
    
    firstObsDate <- allGPRowsForPatient$obsdate[1]
    lastObsDate <-
      allGPRowsForPatient$obsdate[.N]  # Using .N for the last row
    
    numDistinctObservations <- uniqueN(allGPRowsForPatient$obsdate)
    
    firstEpiDate <- allHospRowsForPatient$epistart[1]
    lastEpiDate <- allHospRowsForPatient$epistart[.N]
    
    numDistinctEpisodes <- uniqueN(allHospRowsForPatient$epikey)
    numHospitalisations <- uniqueN(allHospRowsForPatient$spno)
    
    patientDeathDate <- relevantDeath[patid == patientId, dod]
    
    # Directly assign values to the pre-allocated data.table
    postT2Data[i, `:=`(
      patid = patientId,
      firstGPObsDate = firstObsDate,
      numGPObservations = numDistinctObservations,
      lastGPObsDate = lastObsDate,
      firstHospitalEpisode = firstEpiDate,
      numHospitalisations = numHospitalisations,
      numHospitalEpisodes = numDistinctEpisodes,
      lastHospitalEpisode = lastEpiDate,
      deathDate = as.Date(patientDeathDate)
    )]
    if (i %% 1000 == 0)
    {
      message(glue("Processed {i} patients. Please wait..."))
    }
  }
}
#experiment
previousHistoryMetadata <- function() {
  historyInfo <- data.frame(patid = double(), 
                               numberPrevGPObservations = integer(),
                               numberPrevHospVisits = integer())
  
  message("Loading working data into memory...")
  relevantDiagEpi <- data.table((sampledDiagEpi |> filter(patid %in% combinedStrokeData$patid & !(ICD %in% icdCodeMyoInf$icd) & !(ICD %in% icdCodeStroke$icd) &!(ICD %in% icdCodeT2Diab$icd))))
  
  setkey(relevantDiagEpi, "patid")
  relevantObservation <- data.table(observationDataset |> filter(!(medcodeid %in% medcodeMyoInf$medcodeid) & !(medcodeid %in% medcodeStroke$medcodeid) & !(medcodeid %in% medcodeDiabetesT2$medcodeid) & patid %in% combinedStrokeData$patid) |> select(patid, obsdate) |> collect())
  setkey(relevantObservation, "patid")
  relevantPatientIds <- (data.table(combinedStrokeData))[, patid]
  message("Loading completed! Processing...")

  

  for (i in 1:length(relevantPatientIds))
  {
    #message("Getting patient ID...")
    patientId <- as.double(relevantPatientIds[i])
    #message(glue("patient id is {patientId}"))
    #message("Getting first observation date...")
    firstObsDate <- as.Date(combinedStrokeData[patid == patientId, firstGPObsDate])
    #message("Getting first hospital visit date...")
    firstHospVisit <- as.Date(combinedStrokeData[patid == patientId, firstHospitalEpisode])
    #message("Getting relevant observations...")
    obsRecordsForPatient <- relevantObservation[patid == patientId]
    #message("Getting relevant episodes...")
    epiRecordsForPatient <- relevantDiagEpi[patid == patientId]
    #message("Getting number of previous GP records...")
    numOfPrevGPRecords <- nrow(obsRecordsForPatient[obsdate < firstObsDate & obsdate < firstHospVisit])
    hospitalVisitsWithDups <- epiRecordsForPatient[as.Date(epistart) < firstHospVisit & as.Date(epistart) < firstObsDate, spno]
    #numOfPrevHospVisits <- nrow(unique(epiRecordsForPatient[as.Date(epistart) < firstHospVisit & as.Date(epistart) < firstObsDate, spno]))
    #message("Inserting row...")
    

        #                         numberPrevGPObservations = numOfPrevGPRecords, 
    #                         numberPrevHospVisits = length(unique(hospitalVisitsWithDups))))
          
    historyInfo[nrow(historyInfo) + 1,] <- list(patid = patientId, 
                                                numberPrevGPObservations = numOfPrevGPRecords, 
                                                numberPrevHospVisits = length(unique(hospitalVisitsWithDups)))
  
    if (i %% 5000 == 0)
    {
      message(glue("Processed {i} patients. Please wait..."))
    }
    
  }
  
  return(historyInfo)
}


getDescriptiveStatisticMetadata = function(){
  
  setkey(observationDataset, patid)
  setkey(sampledDiagEpi, patid)
  setkey(sampledDeath, patid)
  
  message("Preparing data...")
  # Pre-allocate the output data frame with the correct size
  cohort <- unique(c(combinedStrokeData[,patid], combinedMIData[,patid], combinedDiabData[,patid]))
  numPatients <- length(cohort)
  postT2Data <- data.table(
    patid = numeric(numPatients),
    age = integer(numPatients),
    hasBMI = boolean(numPatients),
    hasHDL = boolean(numPatients),
    hasLDL = boolean(numPatients),
    hasTriglycerides = boolean(numPatients),
    hasGlucose = boolean(numPatients),
    hasSBP = boolean(numPatients),
    hasHbA1C = boolean(numPatients),
    hasSmokerStatus = boolean(numPatients),
    hasCurrentSmoker = boolean(numPatients),
    hasExSmoker = boolean(numPatients),
    hasAlcoholStatus = boolean(numPatients),
    hasHypertension = boolean(numPatients),
    hasHyperglycemia = boolean(numPatients),
    hasAtrialFibrillation = boolean(numPatients)
  )
  
  relevantObservations <- observationDataset[patid %in% cohort]
  relevantPatientData <- sampledPatient[patid %in% cohort]
  
  for (i in seq_len(numPatients))
  {
    patientYob - relevantPatientData[, yob]
  }
}


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
  combined <- intersect(cmdPatientsInObservation$patid, cmdPatientsInHosp)
  
  message("Calculating inverse...")
  inverse <- fullCohort[!patid %in% combined]$patid
  
  newData <- data.table(numCMDPatients = length(combined), 
                        numNonCMDPatients = length(inverse))
  data <- rbind(data, newData)
  
  return(data)
}


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
###
# HDL (mmol/l)
# LDL (mmol/l)
# Triglycerides (mmol/l)
# Total cholesterol (mmol/l)

# Blood glucose (mmol/l)
# SBP (mmHg)
# HbA1C (mmol/l)
# Smoking status
# Current smoker
# Ex- smoker
# Alcohol status
# Hypertension
# Hyperglycemia
# Atrial fibrillation

###
