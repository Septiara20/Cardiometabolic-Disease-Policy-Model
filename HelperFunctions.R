###############################################################################
# Helper functions used throughout DataProcessing.R
###############################################################################

# Accepts two dates in string format (e.g. '2020-01-20')
# Optional format specifier defaulting to the format used in AURUM data
# Returns an integer value with the difference between the dates in days
# Example usage: getDateDiff('2006-01-01', '2005-01-01') returns 365
getDateDiff <- function(date1Str, date2Str, format ='%Y-%m-%d') {
  date1 <- strptime(date1Str, format = format)
  date2 <- strptime(date2Str, format = format)
  diffInDays <- difftime(date1, date2, units = "days")
  return(abs(diffInDays))
}

# Accepts two comparable dataframes, where a column is of a matching type
# Returns all values in listToSearch that are also in listToExclude
# Example usage:
# Assume dataframes A and B both have a column named testVal
# Example usage: getIntersection(A, B, 'testVal')
getIntersection <- function(listToSearch, listToIntersect, columnName) {
  values <- listToSearch |> filter(!!as.symbol(columnName) %in% 
                                   listToIntersect[[columnName]])
  return(values)
}