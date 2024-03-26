###############################################################################
# Septiara Putri, University of Glasgow
# RunPipeline.R
# This is the entry point for the application and is responsible for initiating
# a cold-start run of the end-to-end data processing pipeline.
###############################################################################

source("LoadDependencies.R")
source("LoadData.R")
source("HelperFunctions.R")
source("GetCohort.R")

runPipeline <- function(installPackages = FALSE) {
  loadModules()
  loadDependencies(installPackages)
  
  filteredPatients <- filteredPatients()
  
  strokePats <- getStrokePatients(filteredPatients)
  miPats <- getMIPatients(filteredPatients)
  diabetesPats <- getDiabetesPatients(filteredPatients)
  hyperlipidaemiaPats <- getHyperlipidaemiaPatients(filteredPatients)
  hypertensionPats <- getHypertensionPatients(filteredPatients)
}