#' anticoagulation
#' 
#' @description
#' Matrix anticoagulation
#' 
#' @name anticoagulation
#' @aliases anticoagulation
#' 
#' @docType data
#' 
#' @format A matrix with 5000 simulated odds ratios (rows) for 18 studies
#' (columns).
#' 
#' @source
#' Add reference to source (if available).
#'
#' @keywords datasets
#' 
#' @seealso \code{\link{inc}}
#'
#' @examples
#' data(anticoagulation)
#' dis <- inc(log(anticoagulation), br = 0.3, t1 = 20, sm = "OR")
#' dis

NULL
