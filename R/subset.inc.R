#' Decision Inconsistency index and Across-Studies Inconsistency index in
#' subsets of studies
#'
#' @description
#' Computation of the Decision Inconsistency index and the Across-Studies
#' Inconsistency index for specific subsets of studies (allowing for
#' subgroup analysis).
#'
#' @param x An object of the \code{inc} class.
#' @param data A data frame with the same number of primary studies as those
#'   included in meta-analysis and containing the variables based on which
#'   subset analysis is to be performed.
#' @param subset A logical vector to select studies for subset analysis.
#' @param \dots Additional arguments (ignored).
#'
#' @details
#' This function computes the Decision Inconsistency index and
#' the Across-Studies Inconsistency index for a subset (subgroup) of studies.
#' 
#' @return
#' An object of class \code{\link{inc}}.
#' 
#' @author Bernardo Sousa-Pinto \email{bernardo@@med.up.pt},
#'   Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @examples
#' data(anticoagulation)
#' inc_anticoagulation <-
#'   inc(anticoagulation, dt1 = 16, dt2 = 31, dt3 = 60,
#'     br = 0.5, sm = "OR", transf = FALSE)
#' inc_anticoagulation
#' 
#' # Example with subset analysis restricted to studies with a low risk of bias:
#' data(anticoagulation_df)
#' subset1_anticoagulation <-
#'   subset(inc_anticoagulation, data = anticoagulation_df,
#'     RoB == "Low")
#' subset1_anticoagulation
#'  
#' # Example with subset analysis excluding studies with a high risk of bias
#' # (resulting in the same subset of studies as no study has RoB "Moderate"):
#' subset2_anticoagulation <-
#'    subset(inc_anticoagulation, anticoagulation_df, RoB != "High")
#' subset2_anticoagulation
#'
#' @method subset inc
#' @export

subset.inc <- function(x, data, subset, ...) {
  
  catch <- function(argname, matchcall, data, encl)
    eval(matchcall[[match(argname, names(matchcall))]], data, enclos = encl)
  
  nulldata <- is.null(data)
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  
  if (nulldata)
    data <- sfsp
  
  subset <- catch("subset", mc, data, sfsp)
  
  res <- inc(x$x[, subset, drop = FALSE],
             dt1 = x$dt1, dt2 = x$dt2, dt3 = x$dt3,
             sm = x$sm, br = x$br, scale = x$scale,
             transf = x$transf, transf.dt = x$transf.dt)
  res
}
