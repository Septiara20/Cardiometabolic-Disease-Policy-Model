###############################################################################
# This R script loads patient data for those patients registered on the AURUM 
# system between 2000 and 2020. Medcodes are also loaded into R dataframes.
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

message("Loading observation dataset...")
# Load observation data in Parquet format
observationDataset <- open_dataset('Data/Observation/')

message("Loading medcodes...")
# Load diagnosis medcodes
medcodeDiabetesT2 <- read.csv('Data/Medcodes/medcode_type2diabetes.csv')
medcodeMyoInf <- read.csv('Data/Medcodes/medcode_mi.csv')
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
medcodeSmokingStatus <- read.csv('Data/Medcodes/medcode_smokingstatus.csv')
medcodeTotalChol <- read.csv('Data/Medcodes/medcode_totalclstrl.csv')
medcodeTriglycerides <- read.csv('Data/Medcodes/medcode_triglycerides.csv')