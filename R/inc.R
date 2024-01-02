#' Decision Inconsistency and Across-Studies Inconsistency index
#'
#' @description
#' 
#' Calculate Decision Inconsistency (DI) and Across-Studies Inconsistency (ASI)
#' index
#'
#' @param x An R object of class \code{sims} or a matrix containing the 
#'   simulated effect sizes of primary studies.
#' @param t1 A single numeric defining the decision threshold to distinguish
#'   (i) meaningful from trivial effects, if arguments \code{t2} and \code{t3}
#'   are not provided, (ii) negative from trivial effects, if only argument
#'   \code{t2} is also provided or (iii) small from trivial effects if
#'   arguments \code{t2} and \code{t3} are provided.
#' @param t2 A single numeric defining the decision threshold to distinguish
#'   (i) positive from trivial effects if argument \code{t3} is not provided
#'   or (ii) moderate from small effects, if argument \code{t3} is provided.
#' @param t3 A single numeric defining the decision threshold to distinguish
#'   large from moderate effects.
#' @param sm A character string indicating the summary measure used in
#'   primary studies (see Details).
#' @param br Baseline risk (only considered for odds, risk or hazard ratio).
#' @param scale The number of people per which absolute decision thresholds are
#'  provided (default: 1000, i.e., absolute decision threshold values are
#'  defined per 1000 people). Only considered if \code{br} is not missing.
#' @param transf A logical indicating whether the values of an effect size
#'   matrix (argument \code{x}) are transformed. By default
#'   (\code{transf = TRUE}), it is assumed that the matrix contains log odds
#'   ratios instead of odds ratios, for example.
#' @param transf.dt A logical indicating whether relative decision thresholds are
#'   transformed or on the original scale. If \code{transf.dt = FALSE} (default),
#'   relative decision thresholds are expected to be odds ratios instead of
#'   log odds ratios for \code{sm = "OR"}, for example.
#' @param digits Minimal number of significant digits to print percentages, see
#'   \code{print.default}.
#' @param object R object of class \code{inc}.
#' @param \dots Additional arguments (ignored)
#'
#' @details
#' 
#' This function calculates the Decision Inconsistency index (DI) and the
#' Across-Studies Inconsistency index (ASI) for meta-analyses. The following
#' possibilities are considered depending on the type of effect size measures:
#'   \itemize{
#'     \item Effect size measure corresponding to a ratio (\code{sm = "OR"}, \code{"RR"} or \code{"HR"}) with the DI and the ASI being calculated based on absolute effects: This requires the setting of a baseline risk (i.e., \code{br} must be defined). The decision threshold values (\code{t1}, \code{t2} and \code{t3} must be provided as absolute effects (i.e., number of additional or diminished events per N people. By default, it is assumed that these threshold values are provided per 1000 people. However, this can be changed using the \code{scale} argument).
#'     \item Effect size measure corresponding to a ratio (\code{sm} = \code{"OR"}, \code{"RR"}, \code{"HR"} or \code{"GEN_ratio"}) with the DI and the ASI being calculated based on relative effect size measures: The simulation results of relative effect size measures of primary studies are directly compared with decision thresholds (\code{t1}, \code{t2}, \code{t3}) also expressed as relative effect size measures. This is the adopted approach when no information is provided on the baseline risk (\code{br}).
#'     \item Effect size measure corresponding to a difference (\code{sm = "MD"}, \code{"SMD"}, \code{"RD"} or \code{"GEN_diff"}): The simulation results of the effect size measures of primary studies are directly compared with decision thresholds (\code{t1}, \code{t2}, \code{t3}) also expressed as differences.
#'   }
#'   
#' Of note, when dealing with relative effect size measures, judgements based
#' on absolute effects tend to be considered more important for
#' decision making. The formulae for calculating absolute effects based on
#' relative effect size measures are those used in the GRADE approach
#' (see References below).
#' 
#' Either only argument \code{t1} or arguments  \code{t1}, \code{t2} and
#' \code{t3} must be provided.
#' 
#' Argument \code{sm} must be \code{"OR"} (odds ratio),
#' \code{"RR"} (risk ratio), \code{"HR"} (hazard ratio),
#' \code{"MD"} (mean difference), \code{"SMD"} (standardised mean difference),
#' \code{"RD"} (risk difference), \code{"GEN_diff"} (generic difference), or
#' \code{"GEN_ratio"} (generic ratio).
#' 
#' The baseline risk (optional argument \code{br} which must be a numeric value
#' between 0 and 1) can be provided when \code{sm = "OR"}, \code{"RR"} or
#' \code{"HR"}. The baseline risk is also known as assumed comparator risk,
#' i.e., the risk that the outcome of interest occurs in the comparison
#' intervention.
#'
#' @return
#' An object of class \code{inc}, for which some standard methods are available, see \code{\link{metainc-package}}. Some of the components include:
#' \item{DI}{A percentage corresponding to the Decision Inconsistency index. The higher/closer to 100\% the value, the higher the inconsistency.}
#' \item{ASI}{A percentage corresponding to the Across-Studies Inconsistency index. The higher/closer to 100\% the value, the higher the across-studies inconsistency.}
#' \item{class_distribution}{A data frame containing the proportion of simulations indicating:
#'     \itemize{
#'       \item Large positive effects (effect sizes higher than \code{tl}): “large (higher)” row;
#'       \item Moderate positive effects (efect sizes between \code{tm} and \code{tl}): “moderate (higher)” row;
#'       \item Small positive effects (effect sizes between \code{ts} and \code{tm}): “small (higher)” row;
#'       \item Non meaningful effects (effect sizes between \code{-ts} and \code{ts}): “not meaningful” row;
#'       \item Small negative effects (effect sizes between \code{-ts} and \code{-tm}): “small (lower)” row;
#'       \item Moderate negative effects (effect sizes between \code{-tm} and \code{-tl}): “moderate (lower)” row;
#'       \item Large negative effects (effect sizes lower than \code{-tl}): “large (lower)” row.
#'     }
#' }
#' \item{prop_over_null}{A numeric value indicating the proportion of simulations with a value higher than the value representing no difference between the groups.}
#' 
#' @author Bernardo Sousa-Pinto \email{bernardo@@med.up.pt},
#'   Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @references
#' 
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
#' @examples
#' # Example with effect sizes measures expressed as ratios and with
#' # calculation of the Decision Inconsistency index and the Across-Studies
#' # Inconsistency index being calculated based on absolute effects:
#' 
#' data(anticoagulation)
#' inc_anticoagulation <-
#'   inc(anticoagulation, t1 = 20, t2 = 30, t3 = 40, br = 0.5, sm = "or",
#'       transf = FALSE)
#' inc_anticoagulation
#' 
#' \dontrun{
#' # Same result
#' inc_anticoagulation <-
#'   inc(log(anticoagulation), t1 = 20, t2 = 30, t3 = 40, br = 0.5, sm = "or")
#' inc_anticoagulation
#' 
#' # Example with calculation of the Decision Inconsistency index and the 
#' # Across-Studies Inconsistency index based on effect size measures expressed
#' # as mean differences:
#' 
#' data(montelukast)
#' inc_montelukast <- inc(montelukast, t1 = 0.2, t2 = 0.4, t3 = 0.6, sm = "md")
#' inc_montelukast
#' }
#'
#' @export inc
#' 
#' @importFrom meta gs transf backtransf
#' @importFrom stats reshape rnorm update
#' @importFrom utils head
#' @importFrom confintr ci_cramersv

inc <- function(x, t1, t2, t3, sm, br = NULL, scale = 1000,
                transf = TRUE, transf.dt = FALSE) {
  
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
    #
    chklogical(transf)
    #
    if (!transf)
      simdat <- transf(simdat, sm)
  }
  #
  only.t1 <- !missing(t1) & missing(t2) & missing(t3)
  only.t12 <- !(missing(t1) | missing(t2)) & missing(t3)
  avail.t123 <- !(missing(t1) | missing(t2) | missing(t3))
  #
  if (!only.t1 & !only.t12 & !avail.t123)
    stop("You must provide an input for either (i) argument 't1' or ",
         "(ii) argument 't1' and 't2' or ",
         "(iii) arguments 't1', 't2' and 't3'.",
         call. = FALSE)
  #
  if (only.t12) {
    chknumeric(t1, max = 0, zero = TRUE, length = 1)
    chknumeric(t2, min = 0, zero = TRUE, length = 1)
  }
  else
    chknumeric(t1, min = 0, zero = TRUE, length = 1)
  #
  if (avail.t123) {
    chknumeric(t2, min = 0, zero = TRUE)
    chknumeric(t3, min = 0, zero = TRUE)
    #
    if (t2 <= t1 | t3 <= t1 | t3 <= t2)
      stop("Values provided for arguments 't1', 't2' and 't3' ",
           "must be in increasing order: t1 < t2 < t3")
  }
  #
  avail.br <- !is.null(br)
  if (avail.br && (br <= 0 | br > 1))
    stop("Baseline risk (argument 'br') must be a value between 0 and 1.")
  #
  chknumeric(scale, min = 0, zero = TRUE, length = 1)
  #
  chklogical(transf.dt)
  
  
  #
  # (2) Transform / back transform effect sizes and decision thresholds
  #
  
  if (!transf.dt && !(sm %in% c("OR", "RR", "HR") & avail.br)) {
      t1 <- transf(t1, sm)
      #
      if (only.t12 | avail.t123)
        t2 <- transf(t2, sm)
      if (avail.t123)
        t3 <- transf(t3, sm)
  }
  #
  if ((is_relative_effect(sm) & !avail.br) || sm == "GEN_ratio")
    warning(
      paste("Please note that judgements based on absolute effects",
            "tend to be considered more important for decision-making"),
            call. = FALSE)
  else if (sm %in% c("OR", "RR", "HR") & avail.br) {
    simdat <- backtransf(simdat, sm)
    #
    if (sm == "OR")
      simdat <- -scale * (br - ((br * simdat) / (1 - br + (br * simdat))))
    else if (sm == "RR")
      simdat <- scale * (br * simdat - br)
    else
      simdat <- -scale * (br - (1 - (1 - br)^simdat))
  }
  #
  if (only.t1) {
    n.cat <- 3
    labs <- c("lower", "trivial", "higher")
    if (t1 == 0)
      t1 <- 1 + 1e-12
    cuts <- sort(c(-t1, t1))
  }
  else if (only.t12) {
    n.cat <- 3
    labs <- c("lower", "trivial", "higher")
    cuts <- sort(c(t1, t2))
  }
  else {
    n.cat <- 7
    labs <- c("large (lower)", "moderate (lower)", "small (lower)",
              "not meaningful",
              "small (higher)", "moderate (higher)", "large (higher)")
    #
    if (t1 == 0)
      t1 <- 1 + 1e-12
    cuts <- sort(c(-t3, -t2, -t1, t1, t2, t3))
  }
  
  #
  # (3) Calculate DI and ASI
  #
  
  simdat <- as.data.frame(simdat)
  N <- ncol(simdat)
  n.obs <- prod(dim(simdat))
  #
  propnull <- sum(simdat > 0) / n.obs
  #
  dtcut <- function(x, cutpoints, labels)
    cut(x, c(-Inf, cutpoints, Inf), labels = labels)
  #
  simdat.dt <-
    as.data.frame(lapply(simdat, dtcut, cutpoints = cuts, labels = labs))
  colnames(simdat.dt) <- colnames(simdat)
  #
  tab.dt <- do.call("rbind", lapply(simdat.dt, table))
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
    if (!avail.t123) {
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
  dtVec <- factor(as.vector(unlist(simdat.dt)), levels = labs)
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
              stud_class = simdat.dt,
              t1 = t1,
              t2 = if (only.t12 | avail.t123) t2 else NA,
              t3 = if (avail.t123) t3 else NA,
              br = if (avail.br) br else NA,
              scale = scale, sm = sm,
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
  cat(paste("Proportion of simulations by classification in relation to",
            "the decision threshold:\n"))
  props <- x$class_distribution
  props$proportion <- paste0(formatN(round(100 * props$proportion, digits),
                                     digits), "%")
  props$n <- props$d <- NULL
  prmatrix(props, quote = FALSE, right = TRUE)
  #
  cat(paste0("\nProportion of simulations larger than null effect: ",
             formatN(round(100 * x$prop_over_null, digits), digits),
             "%", "\n"))
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
  cat(paste0("Decision Inconsistency index (DI):  ",
             formatN(round(x$DI, digits), digits), "%\n"))
  #
  cat(paste0("Across-Studies Inconsistency (ASI): ",
             formatN(round(x$ASI, digits), digits), "%\n\n"))
  #
  cat(paste("Proportion of simulations by classification in relation to",
            "the decision threshold:\n"))
  props <- x$class_distribution
  props$proportion <- paste0(formatN(round(100 * props$proportion, digits),
                                     digits), "%")
  props$n <- props$d <- NULL
  prmatrix(props, quote = FALSE, right = TRUE)
  #
  cat(paste0("\nProportion of simulations larger than null effect: ",
             formatN(round(100 * x$prop_over_null, digits), digits),
             "%", "\n"))
  #
  invisible(NULL)
}