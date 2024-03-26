################################################################################
# Cardiometabolic disease policy model in the UK setting
# Septiara Putri 
# Health Economics and Health Technology Assessment (HEHTA), 
# University of Glasgow, UK [2024]
################################################################################
# Cohort development
# Inclusion criteria: adult patients >=18 years old, registered between 2000-2020, linked to HES, IMD and ONS data
# This R script filters the usable cohort from CPRD AURUM data based on inclusion criteria
 ################################################################################


# Patient eligibility criteria
filterPatients <- function() {
  #message("Filtering out patients registered for under a year...")
  # Extract patients who were registered for over two years
  #regOver2Y <- sampledPatient |>
  #  filter(getDateDiff(regenddate, regstartdate) > (365*2)) |>
  #  select(patid)
  
  message("Filtering out patients who were under 18 when first registered...")
  # Extract the patient IDs who were over 18 when first registered
  over18WhenReg <- sampledPatient |> 
    filter(getDateDiff(regstartdate, glue("{yob}-01-01}")) > (365*18)) %>% 
    select(patid)
  
  # Find the patients who were both under 18 and registered for less than 2y
  #excludedByBoth <- sampledPatient |>
  # filter(!patid %in% over18WhenReg$patid & 
  #           !patid %in% regOver2Y$patid)
  
  #combinedFilter <- getIntersection(over18WhenReg, regOver2Y, 'patid')
  
  # Log out the number of patients discarded due to their age
  message(glue("Number of patients discarded due to being <18y when registered: 
           {nrow(patientIds) - nrow(over18WhenReg) }"))
  # Log out the number of patients discarded due to their registration period
  #message(glue("Number of patients registered for under a year:
  #         {nrow(patientIds) - nrow(regOver2Y)}"))
  # Log out the number of patients who were both under 18 and registered for <2y
  #message(glue("Number of patients under 18 registered for under two years:
  #         {nrow(excludedByBoth)}"))
  # Overall patient cohort remaining
  #message(glue("Overall number of patients remaining: {nrow(combinedFilter)}"))
  message(glue("Overall number of patients remaining: {nrow(over18WhenReg)}"))
  
  
  return(over18WhenReg)
}


# Cardiometabolic disease (events of interest)

# Type 2 diabetes
getDiabetesPatients = function(filteredPatients)
{
  allDiabetesPatients <- observationDataset |> 
    filter(medcodeid %in% medcodeDiabetesT2$medcodeid) |>
    select(patid) |>
    distinct() |>
    collect()
  
  message(glue("Number of type 2 diabetes patients in total: 
               {nrow(allDiabetesPatients)}"))
  
  filteredDiabetesPatients <- allDiabetesPatients |> 
    filter(patid %in% filteredPatients$patid)
  
  message(glue("Number of type 2 diabetes patients once filtered: 
               {nrow(filteredDiabetesPatients)}"))

  return(filteredDiabetesPatients)
}

  # Stroke - using GP observations
getStrokePatients <- function(filteredPatients)
{
  allStrokePatients <- observationDataset |> 
    filter(medcodeid %in% medcodeStroke$medcodeid) |>
    select(patid) |>
    distinct() |>
    collect()
  
  message(glue("Number of stroke patients in total: 
               {nrow(allStrokePatients)}"))
  
  filteredStrokePatients <- allStrokePatients |> 
    filter(patid %in% filteredPatients$patid)
  
  message(glue("Number of stroke patients once filtered: 
               {nrow(filteredStrokePatients)}"))
  
  return(filteredStrokePatients)
}

# Stroke - using DiagEpi (hospital diagnoses)
getStrokeDiagEpiPatients <- function(filteredPatients)
{
  allStrokePatients <- sampledDiagEpi |> 
    filter(ICD %in% icdCodeStroke$icd) |>
    select(patid) |>
    distinct()
  
  message(glue("Number of stroke patients fr in total: 
               {nrow(allStrokePatients)}"))
  ß
  filteredStrokePatients <- allStrokePatients |> 
    filter(patid %in% filteredPatients$patid)
  
  message(glue("Number of stroke patients once filtered: 
               {nrow(filteredStrokePatients)}"))
  
  return(filteredStrokePatients)
}

  # Myocardial Infarction
getMIPatients <-function(filteredPatients)
{
  allMIPatients <- observationDataset |> 
    filter(medcodeid %in% medcodeMyocardialInf$medcodeid) |>
    select(patid) |>
    distinct() |>
    collect()
  
  message(glue("Number of myocardial infarction patients in total: 
               {nrow(allMIPatients)}"))
  
  filteredMIPatients <- allMIPatients |> 
    filter(patid %in% filteredPatients$patid)
  
  message(glue("Number of myocardial infarction patients once filtered: 
               {nrow(filteredMIPatients)}"))
  
  return(filteredMIPatients)
}

  # Death

  
# Potential covariates
  # Age and sex will be tabulated in descriptive statistic

  # Smoking status
getSmokingStatPatients <- function(filteredPatients)
{
  allSmokingStatPatients <- observationDataset |> 
    filter(medcodeid %in% medcodeSmokingStatus$medcodeid) |>
    select(patid) |>
    distinct() |>
    collect()
  
  message(glue("Number of smoking status patients in total: 
               {nrow(allSmokingstatPatients)}"))
  
  filteredSmokingstatPatients <- allSmokingstatPatients |> 
    filter(patid %in% filteredPatients$patid)
  
  message(glue("Number of smoking status patients once filtered: 
               {nrow(filteredSmokingstatPatients)}"))
  
  return(filteredSmokingStatPatients)
}
  
  # Alcohol status
getAlcoholPatients <-   function(filteredPatients)
{
  allAlcoholPatients <- observationDataset |> 
    filter(medcodeid %in% medcodeAlcohol$medcodeid) |>
    select(patid) |>
    distinct() |>
    collect()
  
  message(glue("Number of alcohol use of patients in total: 
               {nrow(allAlcoholPatients)}"))
  
  filteredAlcoholPatients <- allAlcoholPatients |> 
    filter(patid %in% filteredPatients$patid)
  
  message(glue("Number of alcohol use of patients once filtered: 
               {nrow(filteredAlcoholPatients)}"))
  
  return(filteredAlcoholPatients)
}


  # Ethnicity
getEthnicityPatients <- function(filteredPatients)
{
  allEthnicityPatients <- observationDataset |> 
    filter(medcodeid %in% medcodeSmokingStatus$medcodeid) |>
    select(patid) |>
    distinct() |>
    collect()
  
  message(glue("Number of ethnicity group patients in total: 
               {nrow(allSmokingstatPatients)}"))
  
  filteredEthnicityPatients <- allEthnicityPatients |> 
    filter(patid %in% filteredPatients$patid)
  
  message(glue("Number of ethnicity group patients once filtered: 
               {nrow(filteredSmokingstatPatients)}"))
  
  return(filteredEthnicityPatients)
}


  # Atrial Fibrillation
getAFPatients <- function(filteredPatients)
{
  allAFPatients <- observationDataset |> 
    filter(medcodeid %in% medcodeAtrialFibrillation$medcodeid) |>
    select(patid) |>
    distinct() |>
    collect()
  
  message(glue("Number of atrial fibrillation patients in total: 
               {nrow(allAFPatients)}"))
  
  filteredAFPatients <- allAFPatients |> 
    filter(patid %in% filteredPatients$patid)
  
  message(glue("Number of atrial fibrillation patients once filtered: 
               {nrow(filteredAFPatients)}"))
  
  return(filteredAFPatients)
}

  # Hypertension
getHypertensionPatients <- function(filteredPatients)
{
  allHypertensionPatients <- observationDataset |> 
    filter(medcodeid %in% medcodeHypertension$medcodeid) |>
    select(patid) |>
    distinct() |>
    collect()
  
  message(glue("Number of hypertension patients in total: 
               {nrow(allHypertensionPatients)}"))
  
  filteredHypertensionPatients <- allHypertensionPatients |> 
    filter(patid %in% filteredPatients$patid)
  
  message(glue("Number of hypertension patients once filtered: 
               {nrow(filteredHypertensionPatients)}"))
  
  return(filteredHypertensionPatients)
}

  #Hyperlipidaemia
getHyperlipidaemiaPatients <- function(filteredPatients)
{
  allHyperlipidaemiaPatients <- observationDataset |> 
    filter(medcodeid %in% medcodeHyperlipidaemia$medcodeid) |>
    select(patid) |>
    distinct() |>
    collect()
  
  message(glue("Number of hyperlipidaemia patients in total: 
               {nrow(allHyperlipidaemiaPatients)}"))
  
  filteredHyperlipidaemiaPatients <- allHyperlipidaemiaPatients |> 
    filter(patid %in% filteredPatients$patid)
  
  message(glue("Number of hyperlipidaemia patients once filtered: 
               {nrow(filteredHyperlipidaemiaPatients)}"))
  
  return(filteredHyperlipidaemiaPatients)
}

# EHR Biomarker
  #Fasting blood glucose
getBloodGlucosePatients <- function(filteredPatients)
{
  allBloodGlucosePatients <- observationDataset |> 
    filter(medcodeid %in% medcodeBloodGlucose$medcodeid) |>
    select(patid) |>
    distinct() |>
    collect()
  
  message(glue("Number of fasting blood glucose measurement in patients in total: 
               {nrow(allBloodGlucosePatients)}"))
  
  filteredBloodGlucosePatients <- allBloodGlucosePatients |> 
    filter(patid %in% filteredPatients$patid)
  
  message(glue("Number of fasting blood glucose measurement in patients once filtered: 
               {nrow(filteredBloodGlucosePatients)}"))
  
  return(filteredBloodGlucosePatients)
}
  
  #Diastolic Blood Pressure (DBP)
getDBPPatients <- function(filteredPatients)
{
  allDBPPatients <- observationDataset |> 
    filter(medcodeid %in% medcodeDBP$medcodeid) |>
    select(patid) |>
    distinct() |>
    collect()
  
  message(glue("Number of DBP measurement in patients in total: 
               {nrow(allDBPPatients)}"))
  
  filteredDBPPatients <- allDBPPatients |> 
    filter(patid %in% filteredPatients$patid)
  
  message(glue("Number of DBP measurement in patients once filtered: 
               {nrow(filteredDBPPatients)}"))
  
  return(filteredDBPPatients)
}

  #Systolic Blood Pressure (SBP)
getSBPPatients = function(filteredPatients)
{
  allSBPPatients <- observationDataset |> 
    filter(medcodeid %in% medcodeSBP$medcodeid) |>
    select(patid) |>
    distinct() |>
    collect()
  
  message(glue("Number of SBP measurement in patients in total: 
               {nrow(allBloodPressurePatients)}"))
  
  filteredSBPPatients <- allSBPPatients |> 
    filter(patid %in% filteredPatients$patid)
  
  message(glue("Number of SBP measurement in patients once filtered: 
               {nrow(filteredSBPPatients)}"))
  
  return(filteredSBPPatients)
}

  #Body Mass Index
getBMIPatients <- function(filteredPatients)
{
  allBMIPatients <- observationDataset |> 
    filter(medcodeid %in% medcodeBMI$medcodeid) |>
    select(patid) |>
    distinct() |>
    collect()
  
  message(glue("Number of BMI measurement in patients in total: 
               {nrow(allBMIPatients)}"))
  
  filteredBMIPatients <- allBMIPatients |> 
    filter(patid %in% filteredPatients$patid)
  
  message(glue("Number of BMI measurement in patients once filtered: 
               {nrow(filteredBMIPatients)}"))
  
  return(filteredBMIPatients)
}

  #Hba1c
getHba1cPatients <- function(filteredPatients)
{
  allHba1cPatients <- observationDataset |> 
    filter(medcodeid %in% medcodeHba1c$medcodeid) |>
    select(patid) |>
    distinct() |>
    collect()
  
  message(glue("Number of Hba1c measurement in patients in total: 
               {nrow(allHba1cPatients)}"))
  
  filteredHba1cPatients <- allHba1cPatients |> 
    filter(patid %in% filteredPatients$patid)
  
  message(glue("Number of Hba1c measurement in patients once filtered: 
               {nrow(filteredHba1cPatients)}"))
  
  return(filteredHba1cPatients)
}

  #HDL
getHDLPatients = function(filteredPatients)
{
  allHDLPatients <- observationDataset |> 
    filter(medcodeid %in% medcodeHDL$medcodeid) |>
    select(patid) |>
    distinct() |>
    collect()
  
  message(glue("Number of HDL measurement in patients in total: 
               {nrow(allHDLPatients)}"))
  
  filteredHDLPatients <- allHDLPatients |> 
    filter(patid %in% filteredPatients$patid)
  
  message(glue("Number of HDL measurement in patients once filtered: 
               {nrow(filteredHDLPatients)}"))
  
  return(filteredHDLPatients)
}

  #LDL
getLDLPatients <- function(filteredPatients)
{
  allLDLPatients <- observationDataset |> 
    filter(medcodeid %in% medcodeLDL$medcodeid) |>
    select(patid) |>
    distinct() |>
    collect()
  
  message(glue("Number of LDL measurement in patients in total: 
               {nrow(allLDLcPatients)}"))
  
  filteredLDLPatients <- allLDLPatients |> 
    filter(patid %in% filteredPatients$patid)
  
  message(glue("Number of LDLmeasurement in patients once filtered: 
               {nrow(filteredLDLPatients)}"))
  
  return(filteredLDLPatients)
}

  #Total Cholesterol
getTotalcholesterolPatients <- function(filteredPatients)
{
  allTotalcholesterolPatients <- observationDataset |> 
    filter(medcodeid %in% medcodeTotalcholesterol$medcodeid) |>
    select(patid) |>
    distinct() |>
    collect()
  
  message(glue("Number of Totalcholesterol measurement in patients in total: 
               {nrow(allTotalcholesterolcPatients)}"))
  
  filteredTotalcholesterolPatients <- allTotalcholesterolPatients |> 
    filter(patid %in% filteredPatients$patid)
  
  message(glue("Number of Totalcholesterol in patients once filtered: 
               {nrow(filteredTotalcholesterolPatients)}"))
  
  return(filteredTotalcholesterolPatients)
}


# Tryglyceride
getTriglyceridePatients <- function(filteredPatients)
{
  allTriglyceridePatients <- observationDataset |> 
    filter(medcodeid %in% medcodeTriglyceride$medcodeid) |>
    select(patid) |>
    distinct() |>
    collect()
  
  message(glue("Number of Triglyceride measurement in patients in total: 
               {nrow(allTriglyceridecPatients)}"))
  
  filteredTriglyceridePatients <- allTriglyceridePatients |> 
    filter(patid %in% filteredPatients$patid)
  
  message(glue("Number of Triglyceride measurement in patients once filtered: 
               {nrow(filteredTriglyceridePatients)}"))
  
  return(filteredTriglyceridePatients)
}



# All descriptive statistic

# Prevalence and Incidence (age, sex standarized)