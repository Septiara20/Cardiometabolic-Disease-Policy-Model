###############################################################################
# This R script loads patient data for those patients with HES-complete data.
# Medcodes are also loaded into R dataframes.
# The relevant CSV files must be provided.

# Note: CSV files are reproducible by using the sampled patient ID index.
###############################################################################
message("Loading data...")
# Load smaller tables directly into dataframes
message("Loading smaller tables...")
sampledPatient <- read.csv('Data/sampled-patient.csv') 
sampledIMD2010 <- read.csv('Data/sampled-patient-imd2010.csv')
sampledDeath <- read.csv('Data/sampled-aurum-death.csv')
sampledPrimaryDiag <- read.csv('Data/sampled-primary-diag-hosp.csv')
sampledDiagEpi <- read.csv('Data/sampled-diag-epi.csv')
sampledDiagHosp <- read.csv('Data/sampled-diag-hosp.csv')
sampledEpiHes <- read.csv('Data/sampled-epi-hes.csv')

setDT(sampledDeath)
setDT(sampledIMD2010)
setDT(sampledDiagEpi)
setDT(sampledPrimaryDiag)
setDT(sampledEpiHes)
setDT(sampledPatient)
setorder(sampledPatient, patid)
setDT(sampledDiagHosp)

setkey(sampledDeath, patid)
setkey(sampledIMD2010, patid)
setkey(sampledDiagEpi, patid)
setkey(sampledPrimaryDiag, patid)
setkey(sampledEpiHes, patid)
setkey(sampledPatient, patid)

message("Loading observation dataset...")
# Load observation data in Parquet format
observationDataset <- open_dataset('Data/Observation/')

message("Loading medcodes...")
# Load diagnosis medcodes
medcodeDiabetesT2 <- read.csv('Data/Medcodes/medcode_type2diabetes.csv')
medcodeMyoInf <- read.csv('Data/Medcodes/medcode_myocardialinfarction.csv')
medcodeStroke <- read.csv('Data/Medcodes/medcode_stroke.csv')

# Load ICD codelists
icdCodeStroke <- read.csv('Data/ICD/icd10_stroke.csv')
icdCodeMyoInf <- read.csv('Data/ICD/icd10_mi.csv')
icdCodeT2Diab <- read.csv('Data/ICD/icd10_t2dm.csv')

# Load biomarker medcodes
medcodeAtrialFib <- read.csv('Data/Medcodes/medcode_af.csv')
medcodeAlcoholStatus <- read.csv('Data/Medcodes/medcode_alcohol.csv')
medcodeBMI <- read.csv('Data/Medcodes/medcode_bmi.csv')
medcodeDBP <- read.csv('Data/Medcodes/medcode_dbp.csv')
medcodeFastingGlucose <- read.csv('Data/Medcodes/medcode_fastingglucose.csv')
medcodeHBA1C <- read.csv('Data/Medcodes/medcode_hba1c.csv')
medcodeHDL <- read.csv('Data/Medcodes/medcode_hdl.csv')

# The two below are important indicators for CMD conditions
medcodeHypertension <- read.csv('Data/Medcodes/medcode_hypertension.csv')
medcodeHyperlipidaemia <- read.csv('Data/Medcodes/medcode_hyperlipidaemia.csv')

medcodeLDL <- read.csv('Data/Medcodes/medcode_ldl.csv')
medcodeSBP <- read.csv('Data/Medcodes/medcode_sbp.csv')
medcodeSmokingStatus <- read.csv('Data/Medcodes/medcode_smoking.csv')
medcodeTotalChol <- read.csv('Data/Medcodes/medcode_totalcholesterol.csv')
medcodeTriglycerides <- read.csv('Data/Medcodes/medcode_triglycerides.csv')

medcodeCVD_FH <- read.csv('Data/Medcodes/medcode_fh_cvd.csv')
medcodeDiab_FH <- read.csv('Data/Medcodes/medcode_fh_diabetes.csv')
setDT(medcodeCVD_FH)
setDT(medcodeDiab_FH)
setDT(medcodeDiabetesT2)
setDT(medcodeMyoInf)
setDT(medcodeStroke)
setDT(medcodeAtrialFib)
setDT(medcodeAlcoholStatus)
setDT(medcodeBMI)
setDT(medcodeDBP)
setDT(medcodeFastingGlucose)
setDT(medcodeHBA1C)
setDT(medcodeHDL)
setDT(medcodeHypertension)
setDT(medcodeHyperlipidaemia)
setDT(medcodeLDL)
setDT(medcodeSBP)
setDT(medcodeSmokingStatus)
setDT(medcodeTotalChol)
setDT(medcodeTriglycerides)
# Setting medcodeid as integer64 for each medcode dataset
medcodeDiabetesT2 <- medcodeDiabetesT2[, medcodeid := as.integer64(medcodeid)]
medcodeMyoInf <- medcodeMyoInf[, medcodeid := as.integer64(medcodeid)]
medcodeStroke <- medcodeStroke[, medcodeid := as.integer64(medcodeid)]
medcodeAtrialFib <- medcodeAtrialFib[, medcodeid := as.integer64(medcodeid)]
medcodeAlcoholStatus <- medcodeAlcoholStatus[, medcodeid := as.integer64(medcodeid)]
medcodeBMI <- medcodeBMI[, medcodeid := as.integer64(medcodeid)]
medcodeDBP <- medcodeDBP[, medcodeid := as.integer64(medcodeid)]
medcodeFastingGlucose <- medcodeFastingGlucose[, medcodeid := as.integer64(medcodeid)]
medcodeHBA1C <- medcodeHBA1C[, medcodeid := as.integer64(medcodeid)]
medcodeHDL <- medcodeHDL[, medcodeid := as.integer64(medcodeid)]
medcodeHypertension <- medcodeHypertension[, medcodeid := as.integer64(medcodeid)]
medcodeHyperlipidaemia <- medcodeHyperlipidaemia[, medcodeid := as.integer64(medcodeid)]
medcodeLDL <- medcodeLDL[, medcodeid := as.integer64(medcodeid)]
medcodeSBP <- medcodeSBP[, medcodeid := as.integer64(medcodeid)]
medcodeSmokingStatus <- medcodeSmokingStatus[, medcodeid := as.integer64(medcodeid)]
medcodeTotalChol <- medcodeTotalChol[, medcodeid := as.integer64(medcodeid)]
medcodeTriglycerides <- medcodeTriglycerides[, medcodeid := as.integer64(medcodeid)]

medcodeDiab_FH <- medcodeDiab_FH[, medcodeid := as.integer64(medcodeid)]
medcodeCVD_FH <- medcodeCVD_FH[, medcodeid := as.integer64(medcodeid)]
# Ensure keys are set correctly for fast operations
setkey(medcodeDiabetesT2, medcodeid)
setkey(medcodeMyoInf, medcodeid)
setkey(medcodeStroke, medcodeid)
setkey(medcodeAtrialFib, medcodeid)
setkey(medcodeAlcoholStatus, medcodeid)
setkey(medcodeBMI, medcodeid)
setkey(medcodeDBP, medcodeid)
setkey(medcodeFastingGlucose, medcodeid)
setkey(medcodeHBA1C, medcodeid)
setkey(medcodeHDL, medcodeid)
setkey(medcodeHypertension, medcodeid)
setkey(medcodeHyperlipidaemia, medcodeid)
setkey(medcodeLDL, medcodeid)
setkey(medcodeSBP, medcodeid)
setkey(medcodeSmokingStatus, medcodeid)
setkey(medcodeTotalChol, medcodeid)
setkey(medcodeTriglycerides, medcodeid)
