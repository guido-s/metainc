#' Decision Inconsistency and Across-Studies Inconsistency index
#'
#' @description
#' Calculates the Decision Inconsistency (DI) and Across-Studies Inconsistency 
#' (ASI) indices.
#'
#' @param x An R object created with \code{\link{getsamples}} or a matrix
#'   containing the sampled effect sizes of primary studies.
#' @param dt1 A single numeric defining the decision threshold to
#'   distinguish (i) meaningful from trivial effects, if arguments
#'   \code{dt2} and \code{dt3} are not provided, (ii) negative / harmful
#'   from trivial effects, if only argument \code{dt2} is also
#'   provided, or (iii) small from trivial effects if arguments
#'   \code{dt2} and \code{dt3} are provided.
#' @param dt2 A single numeric defining the decision threshold to
#'   distinguish (i) positive / beneficial from trivial effects if
#'   argument \code{dt3} is not provided, or (ii) moderate from small
#'   effects if argument \code{dt3} is provided.
#' @param dt3 A single numeric defining the decision threshold to
#'   distinguish large from moderate effects.
#' @param sm A character string indicating the summary measure used in
#'   primary studies (see Details).
#' @param br Baseline risk (only considered for odds, risk or hazard
#'   ratio).
#' @param utility Utility value.
#' @param scale The number of people per which absolute decision
#'   thresholds are provided (default: 1000, i.e., absolute decision
#'   threshold values are defined per 1000 persons). Only considered if
#'   \code{br} is not missing.
#' @param transf A logical indicating whether the values of an effect
#'   size matrix (argument \code{x}) are to be transformed. By default
#'   \code{transf = TRUE}, it is assumed that the matrix
#'   contains, e.g., log odds ratios instead of odds ratios.
#' @param transf.dt A logical indicating whether relative decision
#'   thresholds are transformed or on the original scale. If
#'   \code{transf.dt = FALSE} (default), relative decision thresholds
#'   are expected to be on the natural scale (e.g., odds ratios
#'   instead of log odds ratios for \code{sm = "OR"}). Note, the GRADE working
#'   group recommends to use absolute instead of
#'   relative decision thresholds.
#' @param digits Minimal number of significant digits to print
#'   percentages, see \code{print.default}.
#' @param object R object of class \code{inc}.
#' @param \dots Additional arguments (ignored)
#'
#' @details
#' Calculates the Decision Inconsistency index (DI) and
#' the Across-Studies Inconsistency index (ASI) for a meta-analysis. The
#' following possibilities are considered depending on the type of
#' effect size measures:
#' \itemize{
#' \item Effect size measure corresponding to a ratio (\code{sm =
#'   "OR"}, \code{"RR"} or \code{"HR"}) with the DI and the ASI being
#'   calculated based on absolute effects: This requires the specification
#'   of a baseline risk (i.e., \code{br} must be defined). The
#'   decision threshold values (\code{dt1}, \code{dt2} and \code{dt3}
#'   must be provided as absolute effects (i.e., number of additional
#'   or diminished events per N people. By default, it is assumed that
#'   these threshold values are provided per 1000 people. However,
#'   this can be changed using the \code{scale} argument).
#' \item Effect size measure corresponding to a ratio (\code{sm} =
#'   \code{"OR"}, \code{"RR"}, \code{"HR"} or \code{"GEN_ratio"}) with
#'   the DI and the ASI being calculated based on relative effect size
#'   measures: The sampled effect sizes of primary studies are
#'   directly compared with decision thresholds (\code{dt1},
#'   \code{dt2}, \code{dt3}) also expressed as relative effect sizes.
#'   This is the adopted approach when no information is
#'   provided on the baseline risk (\code{br}).
#' \item Effect size measure corresponding to a difference (\code{sm}
#'   = \code{"MD"}, \code{"SMD"}, \code{"RD"} or \code{"GEN_diff"}):
#'   The sampled effect sizes of primary studies are directly compared with
#'   decision thresholds (\code{dt1}, \code{dt2}, \code{dt3}) also expressed
#'   as differences.
#' }
#'   
#' Of note, when dealing with relative effect size measures,
#' judgements based on absolute effects tend to be considered more
#' important for decision making. The formulae for calculating
#' absolute effects based on relative effect size measures are those
#' used by the GRADE approach (see references below).
#' 
#' Ideally, arguments \code{dt1}, \code{dt2} and \code{dt3} should be
#' provided. If only one decision threshold is available, it is either
#' possible to provide (i) only \code{dt1}, or (ii) both \code{dt1}
#' and \code{dt2} (if the threshold distinguishing clinically relevant
#' benefits vs trivial effects is different from that distinguishing
#' clinically relevant harms vs trivial effects).
#' 
#' Argument \code{sm} must be \code{"OR"} (odds ratio), \code{"RR"}
#' (risk ratio), \code{"HR"} (hazard ratio), \code{"MD"} (mean
#' difference), \code{"SMD"} (standardised mean difference),
#' \code{"RD"} (risk difference), \code{"GEN_diff"} (generic
#' difference), or \code{"GEN_ratio"} (generic ratio).
#' 
#' The baseline risk (\code{br}) must be a numeric value between 0 and
#' 1. It can be provided when \code{sm = "OR"}, \code{"RR"} or
#' '\code{"HR"}. The baseline risk is also known as assumed comparator
#' risk (i.e., the risk that the outcome of interest occurs in the
#' comparison intervention).
#'
#' @return
#' An object of class \code{inc}, for which some standard methods are
#' available, see \code{\link{metainc-package}}. Some of the
#' components include:
#' \item{DI}{A percentage corresponding to the Decision Inconsistency
#'   index. The higher / closer to 100\% the value, the higher the
#'   inconsistency.}
#' \item{ASI}{A percentage corresponding to the Across-Studies
#'   Inconsistency index. The higher / closer to 100\% the value, the
#'   higher the across-studies inconsistency.}
#' \item{class_distribution}{A data frame containing the proportion of
#'   samples indicating (if three decision thresholds had been
#'   provided):
#'   \itemize{
#'     \item Large positive effects (effect sizes higher than
#'       \code{dt3}): "large (higher)" row;
#'     \item Moderate positive effects (efect sizes between \code{dt2}
#'       and \code{dt3}): "moderate (higher)" row;
#'     \item Small positive effects (effect sizes between \code{dt1}
#'       and \code{dt2}): "small (higher)" row;
#'     \item Non meaningful effects (effect sizes between \code{-dt1}
#'       and \code{dt1}): "not meaningful" row;
#'     \item Small negative effects (effect sizes between \code{-dt1}
#'       and \code{-dt2}): "small (lower)" row;
#'     \item Moderate negative effects (effect sizes between
#'       \code{-dt2} and \code{-dt3}): "moderate (lower)" row;
#'     \item Large negative effects (effect sizes lower than
#'       \code{-dt3}): "large (lower)" row.
#'     }
#' }
#' \item{prop_over_null}{A numeric value indicating the proportion of
#'   samples with a value higher than the value representing no
#'   difference between the groups.}
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
#' #' Sousa-Pinto B, Neumann I, Vieira RJ, et al. (2025).
#' \dQuote{Quantitative assessment of inconsistency in meta-analysis using
#' decision thresholds with two new indices.}
#' \emph{Journal of Clinical Epidemiology}, \bold{181},
#' 111725.
#' 
#' @examples
#' 
#' # Example with effect sizes measures expressed as ratios and with
#' # calculation of the Decision Inconsistency index and the Across-Studies
#' # Inconsistency index based on absolute effects:
#' 
#' data(anticoagulation)
#' inc_anticoagulation <-
#'   inc(anticoagulation, dt1 = 16, dt2 = 31, dt3 = 60, br = 0.5, sm = "OR",
#'       transf = FALSE)
#' inc_anticoagulation
#' 
#' # Same result
#' inc_anticoagulation <-
#'   inc(log(anticoagulation), dt1 = 16, dt2 = 31, dt3 = 60,
#'     br = 0.5, sm = "OR")
#' inc_anticoagulation
#' 
#' # Example with calculation of the Decision Inconsistency index and the 
#' # Across-Studies Inconsistency index based on effect size measures expressed
#' # as mean differences:
#' 
#' data(montelukast)
#' inc_montelukast <- inc(montelukast, dt1 = 0.2, dt2 = 0.4, dt3 = 0.6, sm = "md")
#' inc_montelukast
#'
#' @export inc

inc <- function(x, dt1, dt2 = NULL, dt3 = NULL, sm, br = NULL,
                utility = NULL, scale = 1000,
                transf = TRUE, transf.dt = FALSE) {
  
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
    #
    chklogical(transf)
    #
    if (!transf)
      samdat <- transf(samdat, sm)
  }
  #
  chknumeric(scale, min = 0, zero = TRUE, length = 1)
  #
  avail.br <- !is.null(br)
  if (avail.br && (br <= 0 | br > 1))
    stop("Baseline risk (argument 'br') must be a value between 0 and 1.")
  #
  avail.util <- !is.null(utility)
  if (avail.util && (utility <= 0 | utility > 1))
    stop("Utility value (argument 'utility') must be a value between 0 and 1.")
  #
  avail.dt1 <- !missing(dt1) && !is.null(dt1)
  avail.dt2 <- !missing(dt2) && !is.null(dt2)
  avail.dt3 <- !missing(dt3) && !is.null(dt3)
  #
  if (avail.util & avail.br & sm %in% c("OR", "RR", "HR") &
      !avail.dt1 & !avail.dt2 & !avail.dt3) {
    dt1 <- scale * (0.073 - (0.061 * (1 - utility)))
    dt2 <- scale * (0.18 - (0.146 * (1 - utility)))
    dt3 <- scale * (0.338 - (0.271 * (1 - utility)))
    #
    avail.dt1 <- TRUE
    avail.dt2 <- TRUE
    avail.dt3 <- TRUE
  }
  #
  only.dt1 <- avail.dt1 & !avail.dt2 & !avail.dt3
  only.dt12 <- avail.dt1 & avail.dt2 & !avail.dt3
  avail.dt123 <- avail.dt1 & avail.dt2 & avail.dt3
  #
  if (!only.dt1 & !only.dt12 & !avail.dt123)
    stop("You must provide an input for either (i) argument 'dt1' or ",
         "(ii) argument 'dt1' and 'dt2' or ",
         "(iii) arguments 'dt1', 'dt2' and 'dt3'.",
         call. = FALSE)
  #
  if (only.dt12) {
    chknumeric(dt1, max = 0, zero = TRUE, length = 1)
    chknumeric(dt2, min = 0, zero = TRUE, length = 1)
  }
  else
    chknumeric(dt1, min = 0, zero = TRUE, length = 1)
  #
  if (avail.dt123) {
    chknumeric(dt2, min = 0, zero = TRUE)
    chknumeric(dt3, min = 0, zero = TRUE)
    #
    if (dt2 <= dt1 | dt3 <= dt1 | dt3 <= dt2)
      stop("Values provided for arguments 'dt1', 'dt2' and 'dt3' ",
           "must be in increasing order: dt1 < dt2 < dt3")
  }
  #
  chklogical(transf.dt)
  
  
  #
  # (2) Transform / back transform effect sizes and decision thresholds
  #
  
  if (!transf.dt && !(sm %in% c("OR", "RR", "HR") & avail.br)) {
      dt1 <- transf(dt1, sm)
      #
      if (only.dt12 | avail.dt123)
        dt2 <- transf(dt2, sm)
      if (avail.dt123)
        dt3 <- transf(dt3, sm)
  }
  #
  if ((is_relative_effect(sm) & !avail.br) || sm == "GEN_ratio")
    warning(
      paste("Please note that judgements based on absolute effects",
            "tend to be considered more important for decision-making"),
            call. = FALSE)
  else if (sm %in% c("OR", "RR", "HR") & avail.br) {
    samdat <- backtransf(samdat, sm)
    #
    if (sm == "OR")
      samdat <- -scale * (br - ((br * samdat) / (1 - br + (br * samdat))))
    else if (sm == "RR")
      samdat <- scale * (br * samdat - br)
    else
      samdat <- -scale * (br - (1 - (1 - br)^samdat))
  }
  #
  if (only.dt1) {
    n.cat <- 3
    labs <- c("lower", "trivial", "higher")
    if (dt1 == 0)
      dt1 <- 1 + 1e-12
    cuts <- sort(c(-dt1, dt1))
  }
  else if (only.dt12) {
    n.cat <- 3
    labs <- c("lower", "trivial", "higher")
    cuts <- sort(c(dt1, dt2))
  }
  else {
    n.cat <- 7
    labs <- c("large (lower)", "moderate (lower)", "small (lower)",
              "not meaningful",
              "small (higher)", "moderate (higher)", "large (higher)")
    #
    if (dt1 == 0)
      dt1 <- 1 + 1e-12
    cuts <- sort(c(-dt3, -dt2, -dt1, dt1, dt2, dt3))
  }
  
  #
  # (3) Calculate DI and ASI
  #
  
  samdat <- as.data.frame(samdat)
  N <- ncol(samdat)
  n.obs <- prod(dim(samdat))
  #
  propnull <- sum(samdat > 0) / n.obs
  #
  dtcut <- function(x, cutpoints, labels)
    cut(x, c(-Inf, cutpoints, Inf), labels = labels)
  #
  samdat.dt <-
    as.data.frame(lapply(samdat, dtcut, cutpoints = cuts, labels = labs))
  colnames(samdat.dt) <- colnames(samdat)
  #
  tab.dt <- do.call("rbind", lapply(samdat.dt, table))
  #
  sum0 <- function(x)
    sum(x) == 0
  #
  drp <- apply(tab.dt, 2, sum0)
  tab.dt <- tab.dt[, !drp, drop = FALSE]
  nc <- ncol(tab.dt)
  #
  if (nc == 1) {
    ds <- 0
  }
  else {
    suppressWarnings(v <- ci_cramersv(stats::chisq.test(tab.dt)))
    #
    if (!avail.dt123) {
      if (nc == 2 & N > 2)
        ds <- 100 * sqrt(as.numeric(v$estimate)^2 / 2)
      else
        ds <- 100 * as.numeric(v$estimate)
    }
    else {
      if (nc < 7 & N >= 7)
        ds <- 100 * sqrt(as.numeric(v$estimate)^2 / (1 / ((nc - 1) / 6)))
      else if (nc < 7 & N < 7 & nc < N)
        ds <- 100 * sqrt(as.numeric(v$estimate)^2 / (1 / ((nc - 1) /(N - 1))))
      else
        ds <- 100 * as.numeric(v$estimate)
    }
  }
  #
  dtVec <- factor(as.vector(unlist(samdat.dt)), levels = labs)
  n.obs <- length(dtVec)
  #
  tab.dtVec <- table(dtVec)
  n.dt <- as.vector(tab.dtVec)
  #
  cl <- data.frame(n = n.dt,
                   proportion = n.dt / n.obs,
                   d = abs((1 / n.cat) - n.dt / n.obs))
  #
  rownames(cl) <- names(tab.dtVec)
  #
  di <- 100 * (1 - 0.5 * sum(cl$d) / ((n.cat - 1) / n.cat))
  
  #
  # (4) Return results
  #
  
  res <- list(ASI = ds, DI = di, class_distribution = cl,
              prop_over_null = propnull,
              stud_class = samdat.dt,
              dt1 = dt1,
              dt2 = if (only.dt12 | avail.dt123) dt2 else NA,
              dt3 = if (avail.dt123) dt3 else NA,
              br = if (avail.br) br else NA,
              utility = if (avail.util) utility else NA,
              scale = scale, sm = sm,
              transf = transf, transf.dt = transf.dt,
              x = x,
              call = match.call())
  class(res) <- "inc"
  #
  res
}





#' @rdname inc
#' @keywords print
#' @method print inc
#' @export 

print.inc <- function (x, digits = 1, ...) {
  chknumeric(digits, min = 0, length = 1)
  #
  cat(paste0("Decision Inconsistency index (DI):  ",
             formatN(round(x$DI, digits), digits), "%\n"))
  #
  cat(paste0("Across-Studies Inconsistency (ASI): ",
             formatN(round(x$ASI, digits), digits), "%\n\n"))
  #
  cat("Proportion of samples by classification in relation to ",
      "the decision threshold", if (!is.na(x$dt2)) "s", ":\n",
      sep = "")
  props <- x$class_distribution
  props$proportion <- paste0(formatN(round(100 * props$proportion, digits),
                                     digits), "%")
  props$n <- props$d <- NULL
  prmatrix(props, quote = FALSE, right = TRUE)
  #
  cat(paste0("\nProportion of samples larger than null effect: ",
             formatN(round(100 * x$prop_over_null, digits), digits),
             "%", "\n"))
  #
  cat("\nSettings:\n")
  #
  if (!is.na(x$br))
    cat("- baseline risk = ", x$br, "\n", sep = "")
  #
  if (!is.na(x$utility))
    cat("- utility = ", x$utility, "\n", sep = "")
  #
  cat("- decision threshold", if (!is.na(x$dt2)) "s", ": ", sep = "")
  cat("dt1 =", x$dt1)
  if (!is.na(x$dt2))
    cat(", dt2 =", x$dt2)
  if (!is.na(x$dt3))
    cat(", dt3 =", x$dt3)
  #
  if (!is.na(x$br))
    cat(if (!is.na(x$dt2)) "\n  (" else " ",
        "events per ", x$scale, " observations",
        if (!is.na(x$dt2)) ")", "\n",
        sep = "")
  #
  invisible(NULL)
}





#' @rdname inc
#' @method summary inc
#' @export 

summary.inc <- function(object, ...) {
  output <- object
  class(output) <- "summary.inc"
  output
}





#' @rdname inc
#' @keywords print
#' @method print summary.inc
#' @export 

print.summary.inc <- function (x, digits = 1, ...) {
  chknumeric(digits, min = 0, length = 1)
  #
  class(x) <- "inc"
  #
  print(x, digits = digits, ...)
  #
  invisible(NULL)
}
