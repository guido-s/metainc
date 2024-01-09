#' Sensitivity analysis (based on decision thresholds) for the Decision Inconsistency index and Across-Studies Inconsistency index
#'
#' @description
#' Sensitivity analysis on the Decision Inconsistency index and the
#' Across-Studies Inconsistency index based on a range of decision thresholds.
#'
#' @param x An R object of class \code{sims} or a matrix containing the 
#'   simulated effect sizes of primary studies. Note, transformed effect
#'   sizes must be provided, for example, log odds ratios instead of odds
#'   ratios.
#' @param br Baseline risk (only considered for odds, risk or hazard ratio).
#' @param min1 A single numeric defining the lowest value for the lower
#'   decision threshold. Must be lower than \code{max1}, \code{min2} and
#'    \code{max2}.
#' @param max1 A single numeric defining the highest value for the lower
#'   decision threshold. Must be higher than \code{min1}, but lower than
#'   \code{min2} and \code{max2}.
#' @param min2 A single numeric defining the lowest value for the higher
#'   decision threshold. Must be higher than \code{min1} and \code{max1},
#'   but lower than \code{max2}.
#' @param max2 A single numeric defining the highest value for the higher
#'   decision threshold. Must be higher than \code{min1}, \code{max1} and
#'   \code{min2}.
#' @param sm A character string indicating the summary measure used in
#'   primary studies (either \code{sm = "OR"}, \code{sm = "RR"} or 
#'   \code{sm = "HR"}).
#' @param by Increment of the sequences from \code{min1} to \code{max1} and
#'   \code{min2} to \code{max2}.
#' @param scale The number of people per which absolute decision thresholds are
#'  provided (default: 1000, i.e., absolute decision threshold values are
#'  defined per 1000 people).
#' @param limits1 Limits for the colour range in the heatplot showing the
#'   Decision Inconsistency index.
#' @param limits2 Limits for the colour range in the heatplot showing the
#'   Across-Studies Inconsistency index.
#' @param \dots Additional graphical arguments (ignored).
#'
#' @details
#' The \code{\link{inc}} function computes the Decision Inconsistency index
#' (DI) and the Across-Studies Inconsistency index (ASI) for a single
#' baseline risk. This function allows for performing sensitivity analysis
#' according to the decision thresholds.
#' 
#'  Two possibilities are considered:
#'  \itemize{
#'    \item The DI and the ASI are calculated based on absolute effects for meta-analyses of odds ratios, risk ratios or hazard ratios: This requires the setting of a baseline risk (i.e., \code{br} must be defined), as well as the definition of the effect size measure computed for the primary studies (i.e., \code{sm} must be defined either as “odds ratio” [\code{"or"}], “risk ratio” [\code{"rr"}] or “hazard ratio” [\code{"hr"}]). The decision threshold values (\code{min1}, \code{max1}, \code{min2}, \code{max2}) should be provided as absolute effects (i.e., minimum number of additional or diminished events that should occur for the effect to be considered appreciable/important instead of trivial). By default, it is assumed that these threshold values are provided per 1000 people (argument \code{scale}).
#'    \item The DI and the ASI are calculated by directly comparing the effect size measures of primary studies with decision thresholds (\code{min1}, \code{max1}, \code{min2}, \code{max2}). This is the most appropriate approach except when assessing absolute effects for meta-analysis of odds ratios, risk ratios or hazard ratios. The \code{br}, \code{sm} and \code{den} arguments should not be defined.
#'  }
#'
#' @return
#' A dataframe containing
#' \item{dt}{Decision threshold}
#' \item{ASI}{Decision Inconsistency index at baseline risk}
#' \item{DI}{Across-Studies Inconsistency index at Baseline risk}
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
#' @examples
#' data(anticoagulation)
#' dis <- sens_dt(log(anticoagulation),
#'   br = 0.5, min1 = -30, max1 = -10, min2 = 10, max2 = 30, sm = "or")
#' #
#' head(dis)
#' summary(dis$DI)
#' summary(dis$ASI)
#' 
#' heatplot(dis)
#'
#' @export sens_dt

sens_dt <- function(x, br = NULL, min1, max1, min2, max2, sm, by = 1,
                    scale = 1000) {
  
  #
  # (1) Check arguments
  #
  
  if (inherits(x, "sims")) {
    simdat <- x$data
    sm <- x$sm
  }
  else {
    simdat <- x
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
  if (!is.null(br))
    chknumeric(br, min = 0, length = 1)
  #
  chknumeric(min1, length = 1)
  chknumeric(max1, length = 1)
  chknumeric(min2, length = 1)
  chknumeric(max2, length = 1)
  #
  chknumeric(by, length = 1)
  #
  if (min1 >= max1 |
      max1 >= min2 |
      min2 >= max2)
    stop("Invalid decision threshold values. Please make sure that ",
         "decision thresholds are in increasing order: ",
         "min1 < max1 < min2 < max2'.",
         call. = FALSE)
  #
  seq1 <- seq(min1, max1, by)
  seq2 <- seq(min2, max2, by)
  #
  seq1 <- rep(seq1, length(seq1))
  seq2 <- rep(seq2, length(seq2))
  seq2 <- seq2[order(seq2)]
  #
  N <- length(seq1)
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
  asi <- di <- vector("numeric", N)
  classif <- data.frame("higher" = c(), "lower" = c(), "trivial" = c())
  #
  for (i in seq_len(N)) {
    inc.i <- inc(simdat, dt1 = seq1[i], dt2 = seq2[i], sm = sm, br = br,
                 scale = scale)
    #
    asi[i] <- inc.i$ASI
    di[i] <- inc.i$DI
  }
  
  res <- data.frame(dt1 = seq1, dt2 = seq2, ASI = asi, DI = di)
  class(res) <- c("sens_dt", class(res))
  #
  attr(res, "br") <- br
  attr(res, "min1") <- min1
  attr(res, "max1") <- max1
  attr(res, "min2") <- min2
  attr(res, "max2") <- max2
  attr(res, "sm") <- sm
  attr(res, "by") <- by
  attr(res, "scale") <- scale
  #
  res
}





#' @rdname sens_dt
#' @export

heatplot <- function(x, limits1 = NULL, limits2 = NULL, ...) {
  #
  seq1 <- x$dt1
  seq2 <- x$dt2
  dsi <- x$DI
  asi <- x$ASI
  #
  if (is.null(limits1))
    limits1 <- range(dsi, na.rm = TRUE)
  #
  if (is.null(limits2))
    limits2 <- range(asi, na.rm = TRUE)
  #
  p1 <- ggplot(mapping = aes(x = seq1, y = seq2)) +
    geom_raster(aes(fill = dsi), interpolate = TRUE) +
    scale_fill_gradient(high = "red", low = "green", limits = limits1,
                        name = "DI") +
    labs(x = "Lower threshold of appreciable effect",
         y ="Upper threshold of appreciable effect",
         title = "Decision Inconsistency index")
  #
  p2 <- ggplot(mapping = aes(x = seq1, y = seq2)) +
    geom_raster(aes(fill = asi), interpolate = TRUE) +
    scale_fill_gradient(high = "red", low = "green", limits = limits2,
                        name = "V") + 
    labs(x = "Lower threshold of appreciable effect",
         y = "Upper threshold of appreciable effect",
         title = "Across-studies Inconsistency index")
  #
  #p3 <- ggplot(mapping = aes(x = seq1, y = seq2)) + geom_raster(aes(fill = prop_higher),interpolate=T)+scale_fill_gradient(low="white", high="black",limits=c(0,1),name="Proportion") + labs(x = "Lower threshold of appreciable effect", y="Upper threshold of appreciable effect",title="Proportion of simulations higher than \n the highest decision threshold value")
  #p4 <- ggplot(mapping = aes(x = seq1, y = seq2)) + geom_raster(aes(fill = prop_lower),interpolate=T)+scale_fill_gradient(low="white", high="black",limits=c(0,1),name="Proportion") + labs(x = "Lower threshold of appreciable effect", y="Upper threshold of appreciable effect",title="Proportion of simulations lower than \n the lowest decision threshold value")
  #p5 <- ggplot(mapping = aes(x = seq1, y = seq2)) + geom_raster(aes(fill = prop_trivial),interpolate=T)+scale_fill_gradient(low="white", high="black",limits=c(0,1),name="Proportion") + labs(x = "Lower threshold of appreciable effect", y="Upper threshold of appreciable effect",title="Proportion of simulations with trivial effect")
  #
  plot(p1)
  plot(p2)
  #plot(x$p3)
  #plot(x$p4)
  #plot(x$p5)
  #
  invisible(NULL)  
}