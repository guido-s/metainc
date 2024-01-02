#' Simulated effect sizes from Krishnamoorthy et al. (2020)
#' 
#' @description
#' Matrix with simulated effect sizes from Krishnamoorthy et al. (2020).
#' 
#' @name montelukast
#' @aliases montelukast
#' 
#' @docType data
#' 
#' @format A matrix with 5000 simulated mean differences (rows) for 9 studies
#' (columns) following a Bayesian meta-analysis.
#' 
#' @details
#' The \code{montelukast} data is from Krishnamoorthy et al. (2020) (Figure 7)
#' and displays results presented as mean differences. Each column corresponds
#' to a different primary study.
#' 
#' @source
#' Krishnamoorthy, M., Mohd Noor, N., Mat Lazim, N., & Abdullah, B. (2020).
#' \dQuote{Efficacy of montelukast in allergic rhinitis treatment:
#' a systematic review and meta-analysis.}
#' \emph{Drugs},
#' \bold{80}, 1381--1851.
#'
#' @keywords datasets
#' 
#' @seealso \code{\link{inc}}
#'
#' @examples
#' data(montelukast)
#' inc(montelukast, t1 = 0.2, t2 = 0.4, t3 = 0.6, sm = "md")

NULL
