#' Simulated effect sizes from Akl et al. (2017)
#' 
#' @description
#' Matrix with simulated effect sizes from Akl et al. (2017).
#' 
#' @name anticoagulation
#' @aliases anticoagulation
#' 
#' @docType data
#' 
#' @format A matrix with 5000 simulated odds ratios (rows) for 18 studies
#' (columns) following a Bayesian meta-analysis.
#' 
#' @details
#' The \code{anticoagulation} data is from Akl et al. (2017) and displays
#' results presented as mean differences. Each column corresponds
#' to a different primary study.
#' 
#' @source
#' Akl, E. A., Kahale, L. A., Hakoum, M. B., Matar, C. F., Sperati, F.,
#' Barba, M., et al. (2017).
#' \dQuote{Parenteral anticoagulation in ambulatory patients with cancer.}
#' \emph{Cochrane Database of Systematic Reviews},
#' \bold{9}: CD006652.
#'
#' @keywords datasets
#' 
#' @seealso \code{\link{inc}}
#'
#' @examples
#' data(anticoagulation)
#' dis <- inc(log(anticoagulation), br = 0.3, dt1 = 20, sm = "OR")
#' dis

NULL
