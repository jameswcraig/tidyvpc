#' Example observed data from vpc package.
#'
#' An observed dataset from a hypothetical PK model.
#'
#' @format A data frame with 600 rows and 7 variables:
#' \describe{
#'   \item{ID}{Subect identifier}
#'   \item{TIME}{Time}
#'   \item{DV}{Concentration of drug}
#'   \item{AMT}{Amount of dosage initially administered at DV = 0, TIME = 0}
#'   \item{DOSE}{Dosage amount}
#'   \item{AMT}{Amount of drug administered}
#'   \item{MDV}{Dummy indiciating missing dependent variable value}
#'   \item{ISM}{Dummy variable indicating subject's gender (ISM = 0, ISM = 1)}
#' }
#' @source \code{\link[vpc]{simple_data$obs}} 
"exampleobs"

#' Example simulated data from vpc package.
#'
#' A simulated dataset from a hypothetical PK model with 100 replicates.
#'
#' @format A data frame with 60000 rows and 10 variables:
#' \describe{
#'   \item{ID}{Subect identifier}
#'   \item{TIME}{Time}
#'   \item{DV}{Concentration of drug}
#'   \item{IPRED}{Individual prediction variable}
#'   \item{PRED}{Population prediction variable}
#'   \item{AMT}{Amount of dosage initially administered at DV = 0, TIME = 0}
#'   \item{DOSE}{Dosage amount}
#'   \item{AMT}{Amount of drug administered}
#'   \item{ISM}{Dummy variable indicating subject's gender (ISM = 0, ISM = 1)}
#'   \item{MDV}{Dummy indiciating missing dependent variable value}
#' }
#' @source \code{\link[vpc]{simple_data$sim}} 
"examplesim"

