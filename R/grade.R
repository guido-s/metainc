#' Title ...
#'
#' @description
#' One (or few) sentence(s) description.
#'
#' @param yi ...
#' @param vi ...
#' @param rob ...
#' @param data ...
#' @param br ...
#' @param sm ...
#' @param scale ...
#' @param threshold.pval.Q ...
#' @param threshold.I2 ...
#' @param nsim ...
#' @param keep ...
#' @param \dots Additional arguments (ignored).
#'
#' @details
#' Details ...
#' 
#' @return
#' The \code{grade} function returns a list with
#' \item{ASI}{...}
#' \item{DI}{...}
#' \item{inc.Q}{...}
#' \item{inc.I2}{...}
#' \item{pval.Q}{...}
#' \item{threshold.pval.Q}{...}
#' \item{I2 = I2}{...}
#' \item{threshold.I2}{...}
#' \item{yi}{...}
#' \item{vi}{...}
#' \item{rob}{...}
#' \item{sm}{...}
#' \item{ES}{...}
#' \item{intervention.events}{...}
#' \item{baseline.events}{...}
#' \item{br}{...}
#' \item{scale}{...}
#' 
#' @author Bernardo Sousa-Pinto \email{bernardo@@med.up.pt},
#'   Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @examples
#' data(anticoagulation)
#' inc_anticoagulation <-
#'   inc(log(anticoagulation), dt1 = 20, dt2 = 30, dt3 = 40,
#'     br = 0.5, sm = "or")
#' sens <- sens_inf(inc_anticoagulation)
#' sens
#'
#' @export grade

grade <- function(
    yi, vi, rob = NULL, data = NULL,
    br = NA, sm = "GEN", scale = 1000,
    threshold.pval.Q = 0.1, threshold.I2 = 0.60,
    nsim = 10000,
    keep = TRUE,
    ...
    ) {
  
  # Catch 'yi', 'vi', and 'rob' from data
  nulldata <- is.null(data)
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  #
  if (nulldata)
    data <- sfsp
  #
  yi <- catch("yi", mc, data, sfsp)
  chknull(yi)
  k <- length(yi)
  #
  vi <- catch("vi", mc, data, sfsp)
  chknull(vi)
  chklength(vi, k, "Length of 'yi' and 'vi' differs.")
  rob <- catch("rob", mc, data, sfsp)
  if (!is.null(rob))
    chklength(rob, k, "Length of 'yi' and 'rob' differs.")
  
  # Run meta-analysis
  m1 <- metagen(yi, vi, sm = sm)
  #
  if (!is.null(rob)) {
    m1_low <- update(m1, subset = tolower(rob) == "low")
    m1_no_high <- update(m1,
                         subset = tolower(rob) !=
                           c("high" | "very high" | "critical" | "serious"))
  }
  #
  ES <- as.vector(m1$TE.random)
  
  # 1a) Point estimates vary widely across studies
  
  # 2a) Confidence intervals show minimal or no overlap
  
  # 1) Inconsistency indices
  
  sims <- getsims(m1, nsim = nsim)
  res_inc <- inc(sims, sm = sm, ...)
  
  # 2) Statistical test for heterogeneity shows a low p-value
  pval.Q <- m1$pval.Q

  # 3) I2 is large
  I2 <- m1$I2

  # 4) Calculate absolute risk from baseline risk
  if (tolower(sm) %in% c("or", "rr", "hr", "rd")) {
    if (tolower(sm) %in% c("or", "rr", "hr"))
      ES <- exp(ES)
    #
    if (tolower(sm) == "or")
      intervention.events <-
        scale * ES * br /
        (1 - br + ES * br)
    else if (tolower(sm) == "rr")
      intervention.events <- scale * br * ES
    else if (tolower(sm) == "hr")
      intervention.events <-
        scale - exp(log(1 - br) * ES) * scale
    else
      intervention.events <- scale * (br - ES)
  }
  else
    intervention.events <- NA
  
  res <- list(ASI = res_inc$ASI,
              DI = res_inc$DI,
              inc.Q = pval.Q < threshold.pval.Q,
              inc.I2 = I2 > threshold.I2,
              #
              pval.Q = pval.Q, threshold.pval.Q = threshold.pval.Q,
              I2 = I2, threshold.I2 = threshold.I2,
              #
              yi = yi, vi = vi, rob = rob, sm = sm,
              #
              ES = ES,
              #
              intervention.events = intervention.events,
              baseline.events = scale * br,
              br = br,
              scale = scale)
  #
  names(res)[names(res) == "ES"] <- sm
  #
  if (!keep) {
    res$yi <- NULL
    res$vi <- NULL
    res$rob <- NULL
  }
  #
  res
}
