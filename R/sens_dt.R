#' Sensitivity analysis (based on decision thresholds) for the Decision
#' Inconsistency and Across-Studies Inconsistency index
#'
#' @description
#' Sensitivity analysis on the Decision Inconsistency index and the
#' Across-Studies Inconsistency index based on a range of decision thresholds
#' distinguishing meaningful versus trivial effects. It is applicable only to 
#' meta-analyses with binary outcome data (effect size measures expressed as 
#' risk ratios, odds ratios or hazard ratios).
#'
#' @param x An R object created with \code{\link{getsamples}} or a matrix
#'   containing the samples of the effect sizes of primary
#'   studies. Note that log-transformed effect sizes must be provided
#'   (e.g., log odds ratios instead of odds ratios).
#' @param br Baseline risk.
#' @param min1 A single numeric defining the lowest value for the
#'   lower decision threshold. Must be lower than \code{max1},
#'   \code{min2} and \code{max2}.
#' @param max1 A single numeric defining the highest value for the
#'   lower decision threshold. Must be higher than \code{min1}, but
#'   lower than \code{min2} and \code{max2}.
#' @param min2 A single numeric defining the lowest value for the
#'   higher decision threshold. Must be higher than \code{min1} and
#'   \code{max1}, but lower than \code{max2}.
#' @param max2 A single numeric defining the highest value for the
#'   higher decision threshold. Must be higher than \code{min1},
#'   \code{max1} and \code{min2}.
#' @param sm A character string indicating the summary measure used in
#'   primary studies (either \code{sm = "OR"}, \code{sm = "RR"} or
#'   \code{sm = "HR"}).
#' @param by Increment of the sequences from \code{min1} to
#'   \code{max1} and \code{min2} to \code{max2}.
#' @param scale The number of people per which absolute decision
#'   thresholds are provided (default: 1000, i.e., absolute decision
#'   threshold values are defined per 1000 people).
#' @param limits1 Limits for the colour range in the heatplot showing
#'   the Decision Inconsistency index.
#' @param limits2 Limits for the colour range in the heatplot showing
#'   the Across-Studies Inconsistency index.
#' @param \dots Additional graphical arguments (ignored).
#'
#' @details
#' Computes the Decision Inconsistency index
#' (DI) and the Across-Studies Inconsistency index (ASI) across a range of 
#' decision thresholds distinguishing meaningful vs trivial effects. This
#' function can only be applied to dichotomous outcomes expressed as (log-) odds
#' ratio, risk ratio and hazard ratio.
#' Graphical representations can be obtained using the \code{heatplot} function.
#'
#' @return
#' A data frame containing
#' \item{dt1}{Lower decision threshold}
#' \item{dt2}{Higher decision threshold}
#' \item{ASI}{Decision Inconsistency index for each combination of decision thresholds}
#' \item{DI}{Across-Studies Inconsistency index for each combination of decision thresholds}
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
#' dis <- sens_dt(log(anticoagulation),
#'   br = 0.5, min1 = -20, max1 = -10, min2 = 10, max2 = 20, sm = "or")
#' #
#' head(dis)
#' summary(dis$DI)
#' summary(dis$ASI)
#' 
#' heatplot(dis)
#' }
#'
#' @export sens_dt

sens_dt <- function(x, br = NULL, min1, max1, min2, max2, sm, by = 1,
                    scale = 1000) {
  
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
  # (2) Calculate DIs and AIs indices
  #
  asi <- di <- vector("numeric", N)
  classif <- data.frame("higher" = c(), "lower" = c(), "trivial" = c())
  #
  for (i in seq_len(N)) {
    inc.i <- inc(samdat, dt1 = seq1[i], dt2 = seq2[i], sm = sm, br = br,
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
  #p3 <- ggplot(mapping = aes(x = seq1, y = seq2)) + geom_raster(aes(fill = prop_higher),interpolate=T)+scale_fill_gradient(low="white", high="black",limits=c(0,1),name="Proportion") + labs(x = "Lower threshold of appreciable effect", y="Upper threshold of appreciable effect",title="Proportion of samples higher than \n the highest decision threshold value")
  #p4 <- ggplot(mapping = aes(x = seq1, y = seq2)) + geom_raster(aes(fill = prop_lower),interpolate=T)+scale_fill_gradient(low="white", high="black",limits=c(0,1),name="Proportion") + labs(x = "Lower threshold of appreciable effect", y="Upper threshold of appreciable effect",title="Proportion of samples lower than \n the lowest decision threshold value")
  #p5 <- ggplot(mapping = aes(x = seq1, y = seq2)) + geom_raster(aes(fill = prop_trivial),interpolate=T)+scale_fill_gradient(low="white", high="black",limits=c(0,1),name="Proportion") + labs(x = "Lower threshold of appreciable effect", y="Upper threshold of appreciable effect",title="Proportion of samples with trivial effect")
  #
  plot(p1)
  plot(p2)
  #plot(x$p3)
  #plot(x$p4)
  #plot(x$p5)
  #
  invisible(NULL)  
}
