#' Extract simulation results with effect size measures of primary studies
#'
#' @description
#' Obtains a matrix containing simulation results of effect size measures of
#' the primary studies.
#'
#' @param obj Meta-analysis object (see Details).
#' @param param Parameters for effect sizes (see Details).
#' @param package Character string with name of R package used to create
#'   meta-analysis object \code{obj}.
#' @param nsim Number of simulation runs if meta-analysis had been performed
#'   with \code{\link[metafor]{rma}} from R package \bold{metafor}.
#' @param sm A character string indicating the summary measure used in
#'   primary studies.
#' @param transf A logical indicating whether effect sizes are transformed or
#'   on the original scale. If \code{transf = TRUE} (default), effect estimates
#'   are expected to be log odds ratios instead of odds ratios for
#'   \code{sm = "OR"}, for example.
#'
#' @details
#' This function can be used to create a matrix with simulated effect sizes of
#' primary studies (number of rows equal to the number of simulation runs,
#' number of columns equal to the number of studies in the meta-analysis.
#' 
#' Main input to this function is a meta-analysis object created with
#' \code{bugs} from R package \bold{R2OpenBUGS}, \code{jags} from \bold{R2jags},
#' \code{jags.samples} from \bold{rjags}, \code{brm} from \bold{brms}, or
#' \code{rma} from \bold{metafor}.
#' 
#' Parameter corresponding to the effect size measures of the primary studies
#' (R2OpenBUGS, R2jags, rjags packages) or to the name of the column that
#' contains the unique identification of each primary study in the
#' data frame containing the results of primary studies used for meta-analysis
#' (brms package). This argument is not defined if meta-analysis has been
#' performed using \bold{metafor} package.
#' 
#' Information of the R package used to conduct the meta-analysis can be
#' provided in argument \code{package}. This information must not be provided
#' for an object created with \bold{metafor} or \bold{brms}. The value of
#' \code{package} should be "R2OpenBUGS", "R2jags", or "rjags".
#' 
#' The number of simulation runs (argument \code{nsim}) can be specified for
#' meta-analysis objects created with \bold{metafor}. Default is 10000 runs.
#' 
#'
#' @return
#' A matrix containing simulated effect sizes for primary study.
#'
#' @author Bernardo Sousa-Pinto \email{bernardo@@med.up.pt}
#' 
#' @examples
#' \dontrun{
#' # Example with an object obtained using the R2OpenBUGS package ("es",
#' # whose simulations on the effect sizes of primary studies are those
#' # of the "delta" parameter). If the object had been obtained using a
#' # different package, all remaining arguments would be the same and only
#' # the "package" argument would have a different input.
#' 
#' # In this example, argument \code{transf = TRUE} (default), as simulation
#' # results had been obtained as ln(odds ratio).
#' 
#' load(url("https://raw.github.com/BernardoSousaPinto/metainc_extra_files/main/es.Rdata"))
#' 
#' sims <- getsims(es, param = "delta", package = "R2OpenBUGS", sm = "OR")
#' sims
#' }
#' 
#' # Example with a dataset containing variables indicating the effect sizes
#' # for each primary study (yi) and the respective variances (vi), and for
#' # which frequentist meta-analysis using the meta package is applied.
#' 
#' data("anticoagulation_df")
#' m1 <- meta::metagen(yi, sqrt(vi), sm = "OR", data = anticoagulation_df,
#'   studlab = LETTERS[1:18])
#' set.seed(1090) # Make simulated effect sizes reproducible
#' sims1 <- getsims(m1)
#' sims1
#' 
#' \dontrun{
#' # Same simulated effect sizes using R package metafor (must be installed)
#' 
#' m1f <- metafor::rma(anticoagulation_df, measure = "OR", slab = LETTERS[1:18])
#' set.seed(1090) # Make simulated effect sizes reproducible
#' sims2 <- getsims(m1f)
#' sims2
#' 
#' all.equal(sims1, sims2) # Only difference: package name
#' }
#'
#' @export getsims

getsims <- function(obj, param = NULL, package = NULL, nsim = 10000,
                    sm, transf = TRUE) {
  
  #
  # (1) Check arguments
  #
  if (inherits(obj, "meta"))
    package <- "meta"
  else if (inherits(obj, "rma"))
    package <- "metafor"
  else if (inherits(obj, "brmsfit"))
    package <- "brms"
  else if (is.null(package))
    stop("Please indicate the name of the package for meta-analysis")
  else if (!(tolower(package) %in% c("r2openbugs", "r2jags", "rjags")))
    stop("Argument 'package' must be \"R2OpenBUGS\", \"R2jags\", or ",
         "\"rjags\" as meta-analysis object was not created with R package ",
         "meta, metafor, or brms")
  
  #
  # (2) Extract or generate simulated data
  #
  if (package %in% c("meta", "metafor")) {
    m0 <- blup(obj)
    #
    if (package == "meta") {
      studlab <- obj$studlab
      N <- length(studlab)
      sm <- obj$sm
      #
      yi <- m0$blup
      sei <- m0$se.blup
    }
    else {
      studlab <- obj$slab
      N <- length(studlab)
      sm <- obj$measure
      #
      yi <- m0$pred
      sei <- m0$se
    }
    #
    db <- data.frame("V1" = seq_len(nsim))
    #
    for (i in seq_len(N)) {
      a <- rnorm(nsim, yi[i], sei[i])
      db <- cbind(db, a)
    }
    #
    sims <- db[, 1 + seq_len(N)]
    sims <- as.matrix(sims)
    rownames(sims) <- seq_len(nsim)
    colnames(sims) <- studlab
  }
  else {
    if (missing(sm))
      sm <- ""
    #
    if (package == "R2OpenBUGS") {
      sims <- as.data.frame(obj$sims.list[param])
      sims <- as.matrix(sims)
    }
    else if (package == "R2jags") {
      sims <- as.data.frame(obj$BUGSoutput$sims.list[param])
      sims <- as.matrix(sims)
    }
    else if (package == "rjags") {
      N <- as.numeric(nrow(obj[[param]]))
      Np <- N + 1
      sims2 <- data.frame(study = vector("numeric", 0),
                          sims0 = vector("numerci", 0))
      for (i in 1:N) {
        sims0 <- obj[[param]][i, , ]
        NN <- length(sims0)
        sims1 <- data.frame("n_sim" = 1:NN, study = i, sims0)
        sims2 <- rbind(sims2, sims1)
      }
      #
      sims3 <- reshape(data = sims2, direction = "wide",
                       idvar = "n_sim", timevar = "study")
      sims3 <- sims3[, 2:Np]
      sims <- as.matrix(sims3)
    }
    else if (package == "brms") {
      sims_df <- as.data.frame(obj$fit)
      sims_l <- sims_df[, grepl(paste0("r_", param), names(sims_df))]
      sims <- as.matrix(sims_l)
    }
    else
      stop("Please note that the R package for Bayesian analysis ",
           "should either be 'R2OpenBUGS', 'R2jags', 'rjags' or 'brms'")
  }
  #
  if (!transf)
    sims <- transf(sims, sm)
  
  #
  # (3) Return results
  #
  res <- list(data = sims, param = param, package = package, nsim = nsim,
              sm = sm, transf = transf)
  #
  class(res) <- "sims"
  #
  res
}


#' @method print sims
#' @export

print.sims <- function(x, n = 5, ...) {
  chkclass(x, "sims")
  #
  print(head(x$data, n), ...)
  #
  invisible(NULL)
}
