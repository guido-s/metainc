#' metainc: Brief overview of methods and general hints
#' 
#' R package \bold{metainc} allows for an assessment of inconsistency in
#' meta-analysis by calculating the Decision Inconsistency index (DI) and
#' the Across-Studies Inconsistency (ASI) index. These indices quantify
#' inconsistency taking into account outcome-level decision thresholds.
#'
#' @name metainc-package
#' 
#' @docType package
#'
#' @details
#' R package \bold{metainc} provides the following methods:
#' \itemize{
#' \item Generate or extract samples of effect sizes from primary studies
#'   included in a meta-analysis (\code{\link{getsamples}})
#' \item Calculate the Decision Inconsistency and Across-Studies Inconsistency
#'   index (\code{\link{inc}})
#' }
#' 
#' In addition, the following methods for sensitivity and subset analysis are
#' available:
#' \itemize{
#' \item Sensitivity analysis based on the baseline risk (\code{\link{sens_br}})
#' \item Sensitivity analysis based on decision thresholds
#'   (\code{\link{sens_dt}})
#' \item Leave-one-out analysis (\code{\link{sens_inf}})
#' \item Subset analysis (\code{\link{subset.inc}})
#' }
#' 
#' Type \code{help(package = "metainc")} for a listing of R functions
#' available in \bold{metainc}.
#' 
#' Type \code{citation("metainc")} on how to cite \bold{metainc} in
#' publications.
#' 
#' The development version of \bold{metainc} is available on GitHub
#' \url{https://github.com/guido-s/metainc}.
#' 
#' @author Bernardo Sousa-Pinto \email{bernardo@@med.up.pt},
#'   Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @references
#' Cohen J. (1998).
#' \dQuote{Statistical Power Analysis in the Behavioral Sciences},
#' 2nd edition ed. Hillsdale (NJ): Lawrence Erlbaum Associates, Inc.
#' 
#' Schunemann HJ, Higgins JPT, Vist GE, et al. (2019).
#' \dQuote{Completing ‘Summary of findings’ tables and grading the certainty of
#' the evidence.}
#' \emph{Cochrane Handbook for Systematic Reviews of Interventions},
#' 375--402.
#' 
#' Schunemann HJ, Vist GE, Higgins JPT, et al. (2019).
#' \dQuote{Interpreting results and drawing conclusions.}
#' \emph{Cochrane Handbook for Systematic Reviews of Interventions},
#' 403--431.
#' 
#' Skoetz N, Goldkuhle M, van Dalen EC, et al. (2020).
#' \dQuote{GRADE guidelines 27: how to calculate absolute effects for
#' time-to-event outcomes in summary of findings tables and Evidence Profiles.}
#' \emph{Journal of Clinical Epidemiology}, \bold{118},
#' 124--131.
#'
#' @keywords package
#'
#' @importFrom meta blup backtransf gs metagen transf
#' @importFrom stats reshape rnorm update
#' @importFrom utils head
#' @importFrom confintr ci_cramersv
#' @importFrom ggplot2 aes geom_raster ggplot labs scale_fill_gradient
#' @importFrom graphics lines mtext polygon title 


NULL
