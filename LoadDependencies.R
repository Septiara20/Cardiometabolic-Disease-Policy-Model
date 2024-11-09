###############################################################################
# This R script loads all dependencies required to use the modelling project
# Created for convenience
###############################################################################
loadDependencies <- function() {
  message("Loading dependencies...")
  library("arrow", 
          "dplyr", 
          "glue",
          "data.table")
  message("Dependencies loaded!")
}

