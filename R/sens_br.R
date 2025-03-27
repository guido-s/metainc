#' Sensitivity analysis (based on the baseline risk) for the Decision
#' Inconsistency and Across-Studies Inconsistency index
#'
#' @description
#' Sensitivity analysis on the Decision Inconsistency index and the
#' Across-Studies Inconsistency index based on a range of baseline
#' risks. It is applicable only to meta-analyses with a binary outcome
#' (effect size measures expressed as risk ratios, odds ratios or
#' hazard ratios).
#'
#' @param x An R object created with \code{\link{getsamples}} or a matrix
#'   containing sampled effect sizes of primary studies. Note, log
#'   transformed effect sizes must be provided (e.g., log odds ratios
#'   instead of odds ratios).
#' @param br1 Smallest baseline risk considered.
#' @param br2 Largest baseline risk considered.
#' @param dt1 A single numeric defining the decision threshold to
#'   distinguish (i) meaningful from trivial effects, if arguments
#'   \code{dt2} and \code{dt3} are not provided, or (ii) small from
#'   trivial effects if arguments \code{dt2} and \code{dt3} are
#'   provided.
#' @param dt2 A single numeric defining the decision threshold to
#'   distinguish moderate from small effects provided.
#' @param dt3 A single numeric defining the decision threshold to
#'   distinguish large from moderate effects.
#' @param sm A character string indicating the summary measure used in
#'   primary studies (either \code{sm = "OR"}, \code{sm = "RR"} or
#'   \code{sm = "HR"}).
#' @param by Increment of the sequence from \code{br1} to \code{br2}.
#' @param scale The number of people per which absolute decision
#'   thresholds are provided (default: 1000, i.e., absolute decision
#'   threshold values are defined per 1000 people).
#' @param ylim1 The y limits (min, max) of the plot showing the
#'   Decision Inconsistency index.
#' @param ylim2 The y limits (min, max) of the plot showing the
#'   Across-Studies Inconsistency index.
#' @param ylab1 A label for the y-axis (Decision Inconsistency index).
#' @param ylab2 A label for the y-axis (Across-Studies Inconsistency
#'   index).
#' @param \dots Additional graphical arguments (ignored).
#'
#' @details
#' Computes the Decision Inconsistency index (DI) and the
#' Across-Studies Inconsistency index (ASI) across a range of baseline risks.
#' It can only be applied for meta-analyses with binary outcome data (effect
#' size measures expressed as (log) risk ratios, odds ratios or hazard ratios),
#' with the DI and the ASI being calculated based on absolute effects.
#' As a result, the decision threshold values (\code{dt1}, \code{dt2},
#' \code{dt3}) must be provided as absolute effects. By default, it is assumed
#' that threshold values are provided as numbers of events per 1000
#' persons (\code{scale = 1000}).
#'
#' @return
#' A data frame containing
#' \item{br}{Baseline risk}
#' \item{ASI}{Decision Inconsistency index at baseline risk}
#' \item{DI}{Across-Studies Inconsistency index at baseline risk}
#' 
#' @author Bernardo Sousa-Pinto \email{bernardo@@med.up.pt},
#'   Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @references
#' Schunemann HJ, Higgins JPT, Vist GE, et al. (2019).
#' \dQuote{Completing ‘Summary of findings’ tables and grading the certainty of
#' the evidence.}
#' \emph{Cochrane Handbook for Systematic Reviews of Interventions},
#' 375--402.
#' 
#' Skoetz N, Goldkuhle M, van Dalen EC, et al. (2020).
#' \dQuote{GRADE guidelines 27: how to calculate absolute effects for
#' time-to-event outcomes in summary of findings tables and Evidence Profiles.}
#' \emph{Journal of Clinical Epidemiology}, \bold{118},
#' 124--131.
#' 
#' Sousa-Pinto B, Neumann I, Vieira RJ, et al. (2025).
#' \dQuote{Quantitative assessment of inconsistency in meta-analysis using
#' decision thresholds with two new indices.}
#' \emph{Journal of Clinical Epidemiology}, \bold{181},
#' 111725.
#' 
#' @examples
#' \donttest{
#' data(anticoagulation)
#' dis <- sens_br(log(anticoagulation),  br1 = 0.3, br2 = 0.7, dt1 = 20,
#'   sm = "OR", by = 0.1)
#' dis
#' plot(dis, ylim1 = c(0, 100), ylim2 = c(0, 50))
#' }
#'
#' @export sens_br

sens_br <- function(x, br1, br2, dt1, dt2 = NULL, dt3 = NULL,
                    sm, by = 0.01, scale = 1000) {
  
  #
  # (1) Check arguments
  #
  
  if (inherits(x, "samples_metainc")) {
    samdat <- x$data
    sm <- x$sm
  }
  else {
    samdat <- x
    #
    if (missing(sm))
      stop("Argument 'sm' must be provided.")
    #
    sm <- setchar(sm,
                  c("OR", "RR", "HR", "MD", "SMD", "RD",
                    "GEN_ratio", "GEN_diff"),
                  stop.at.error = FALSE, return.NULL = FALSE,
                  nchar.equal = TRUE)
  }
  #
  chknumeric(dt1, length = 1)
  chknumeric(dt2, length = 1)
  chknumeric(dt3, length = 1)
  #
  chknumeric(br1, min = 0, length = 1)
  chknumeric(br2, min = 0, length = 1)
  #
  chknumeric(by, length = 1)
  #
  seq_br <- seq(br1, br2, by)
  #
  chknumeric(scale, min = 0, zero = TRUE, length = 1)
  #
  if (!(sm %in% c("OR", "RR", "HR")))
    stop("This function is only available for summary measures ",
         "\"OR\", \"RR\" and \"HR\".",
         call. = FALSE)
  
  #
  # (2) Calculate DIS and AI indices
  #
  asi <- di <- vector("numeric", length(seq_br))
  classif <- data.frame("higher" = c(), "lower" = c(), "trivial" = c())
  #
  for (i in seq_along(seq_br)) {
    inc.i <- inc(samdat, dt1 = dt1, dt2 = dt2, dt3 = dt3,
                 sm = sm, br = seq_br[i], scale = scale)
    #
    asi[i] <- inc.i$ASI
    di[i] <- inc.i$DI
  }

  res <- data.frame(br = seq_br, ASI = asi, DI = di)
  class(res) <- c("sens_br", class(res))
  #
  attr(res, "br1") <- br1
  attr(res, "br2") <- br2
  attr(res, "dt1") <- dt1
  attr(res, "dt2") <- dt2
  attr(res, "dt3") <- dt3
  attr(res, "sm") <- sm
  attr(res, "by") <- by
  attr(res, "scale") <- scale
  #
  res
}





#' @rdname sens_br
#' @method plot sens_br
#' @export

plot.sens_br <- function(x, ylim1 = c(0, 100), ylim2 = c(0, 100),
                         ylab1 = "DI index (%)",
                         ylab2 = "ASI index (%)",
                         ...) {
  
  plot(x$br, x$DI, type = "l",
       ylim = ylim1, 
       xlab = "Baseline risk",
       ylab = ylab1)
  #
  title("Decision Inconsistency index (DI)")
  
  plot(x$br, x$ASI, type = "l",
       ylim = ylim2,
       xlab = "Baseline risk",
       ylab = ylab2)
  #
  title("Across-Studies Inconsistency index (ASI)")
  if (FALSE) {
    plot(x$br, 1 - classif[, 1], type = "l", ylim = c(0, 1),
         xlab = "Baseline risk", ylab = "Proportion of samples")
    #
    polygon(c(x$br, rev(x$br)), c(x$br * 0, rev(classif[, 2])),
            col = "#6BD7AF", lty = 0)
    polygon(c(x$br, rev(x$br)), c(classif[, 2], rev(1 - classif[, 1])),
            col = "gray85", lty = 0)
    polygon(c(x$br, rev(x$br)), c(1 - classif[, 1], x$br / x$br), col = "red",
            lty = 0)
    lines(x$br, 1 - classif[, 1], ylim=c(0, 1))
    lines(x$br, classif[, 2], ylim = c(0, 1))
    mtext("red: higher; grey: trivial; green: lower", side = 3, adj = 1)
    title(paste("Proportion of samples higher, lower and",
                "within the decision threshold values"))
  }
  #
  invisible(NULL)
}
