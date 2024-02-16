#' Calculate risk ratios from odds ratios and baseline risk
#'
#' @description
#' Calculate risk ratios from odds ratios and baseline risk.
#'
#' @param or Odds ratio(s).
#' @param br Baseline risk (also known as assumed comparator risk),
#'   i.e., risk that the outcome of interest would occur in the comparator 
#'   intervention. It must be expressed as a value between 0 and 1.
#'
#' @details
#' This function converts odds ratios (OR) into risk ratios (RR) using the
#' formula available in Schunemann et al. (2019), Chapter 15:
#' 
#' RR = OR / (1 - br x (1 - OR)),
#' 
#' with br corresponding to the baseline risk (the assumed comparator risk;
#' i.e., the risk that the outcome of interest would occur in the comparison
#' intervention).
#'
#' @return
#' A vector or matrix with risk ratios.
#' 
#' @author Bernardo Sousa-Pinto \email{bernardo@@med.up.pt},
#'   Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @references
#' Schunemann HJ, Vist GE, Higgins JPT, et al. (2019).
#' \dQuote{Interpreting results and drawing conclusions.}
#' \emph{Cochrane Handbook for Systematic Reviews of Interventions},
#' 403--431.
#' 
#' @examples
#' data(anticoagulation)
#' RRs <- or2rr(anticoagulation, br = 0.5)
#' head(RRs)
#'
#' @export or2rr

or2rr <- function(or, br) {
  chknumeric(br, min = 0, max = 1, length = 1)
  or / (1 - br + br * or)
}
