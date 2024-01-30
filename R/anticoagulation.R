#' Samples of effect sizes (odds ratio)
#' 
#' @description
#' Matrix with samples of effect sizes from a meta-analysis by Akl et
#' al. (2017).
#' 
#' @name anticoagulation
#' @aliases anticoagulation
#' 
#' @docType data
#' 
#' @format A matrix with 5000 sampled odds ratios (rows) for 18 studies
#' (columns) following a Bayesian meta-analysis.
#' 
#' @details
#' The \code{anticoagulation} data is from Akl et al. (2017) and
#' displays results presented as odds ratios. Each column corresponds
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
#' 
#' # Since the results are already presented as odds ratios, we need
#' # to indicate that effects have not been transformed to log odds
#' # ratios yet (\code{transf = FALSE}).
#' dis1 <- inc(anticoagulation, br = 0.504, dt1 = 16, sm = "OR",
#'   transf = FALSE)
#' dis1
#' 
#' # Alternatively, we may simply apply the \code{inc()} function to
#' # log odds ratios.
#' data(anticoagulation)
#' dis2 <- inc(log(anticoagulation), br = 0.504, dt1 = 16, sm = "OR")
#' dis2

NULL
