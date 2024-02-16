#' Generate or extract samples of effect sizes from primary studies
#'
#' @description
#' Obtains a matrix containing sampled effect sizes from primary studies
#' of a meta-analysis.
#'
#' @param obj Meta-analysis object (see Details).
#' @param param Parameters for effect sizes (see Details).
#' @param package Character string with name of R package used to
#'   create the meta-analysis object \code{obj}.
#' @param n.samples Number of samples per study if meta-analysis had been
#'   performed with R package \bold{meta} or
#'   \code{\link[metafor]{rma}} from R package \bold{metafor}.
#' @param sm A character string indicating the summary measure used in
#'   primary studies.
#' @param transf A logical indicating whether effect sizes are
#'   transformed or on the original scale. If \code{transf = TRUE}
#'   (default), effect estimates are expected to be on the log scale
#'   (for example, log odds ratios instead of odds ratios for \code{sm
#'   = "OR"}).
#'
#' @details
#' This function can be used to create a matrix with sampled effect sizes of
#' primary studies (number of rows equal to the number of samples,
#' number of columns equal to the number of studies in the meta-analysis).
#' 
#' Main input to this function is argument \code{obj} containing a
#' meta-analysis object created with \code{bugs} from R package
#' \bold{R2OpenBUGS}, \code{jags} from \bold{R2jags},
#' \code{jags.samples} from \bold{rjags}, \code{brm} from \bold{brms},
#' \code{rma} from \bold{metafor}, or \code{metabin}, \code{metacont}
#' or \code{metagen} from \bold{meta}.
#' 
#' Argument \code{param} corresponds to the effect size measures of
#' the primary studies (for R packages \bold{R2OpenBUGS}, \bold{R2jags},
#' and \bold{rjags}) or to the name of the column that contains the unique
#' identification of each primary study in the data frame containing the
#' results of primary studies used for meta-analysis (\bold{brms} package).
#' This argument is not required if the meta-analysis has been performed using
#' R package \bold{meta} or \bold{metafor}.
#' 
#' Information of the R package used to conduct the meta-analysis can
#' be provided in argument \code{package}. The value of \code{package}
#' can be "R2OpenBUGS", "R2jags", "rjags", "meta", "metafor" or
#' "brms". This information is not required for an object created with
#' \bold{meta}, \bold{metafor} or \bold{brms}.
#' 
#' Argument \code{n.samples} contains the number of samples per study for
#' meta-analyses created with \bold{meta} or \bold{metafor}. Default
#' is 10000 samples.
#'
#' Summary measure used in primary studies (argument \code{sm}) can be
#' either "OR" (odds ratio), "RR" (risk ratio), "HR" (hazard ratio),
#' "RD" (risk difference), "MD" (mean difference), "SMD" (standardised
#' mean difference), "GEN_ratio" (generic ratio) or "GEN_diff"
#' (generic difference). Does not need to be provided if the meta-analysis
#' was performed with \bold{meta} or \bold{metafor}.
#' 
#' @return
#' A matrix containing sampled effect sizes (rows) for primary study (columns).
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
#' # In this example, argument \code{transf = TRUE} (default), as sampled
#' # effect sizes are log odds ratios.
#' 
#' load(url("https://raw.github.com/BernardoSousaPinto/metainc_extra_files/main/es.Rdata"))
#' 
#' sample <- getsamples(es, param = "delta", package = "R2OpenBUGS", sm = "OR")
#' sample
#' }
#' 
#' # Example using a dataset providing effect sizes for primary studies (yi)
#' # and respective variances (vi). A frequentist meta-analysis using the meta
#' # package is conducted.
#' 
#' data("anticoagulation_df")
#' m1 <- meta::metagen(yi, sqrt(vi), sm = "OR", data = anticoagulation_df,
#'   studlab = LETTERS[1:18])
#' set.seed(1090) # Make sampled effect sizes reproducible
#' sample1 <- getsamples(m1)
#' sample1
#' 
#' \dontrun{
#' # Same samples of effect sizes using R package metafor (must be installed)
#' 
#' m2 <- metafor::rma(anticoagulation_df, measure = "OR", slab = LETTERS[1:18])
#' set.seed(1090) # Make sampled effect sizes reproducible
#' sample2 <- getsamples(m2)
#' sample2
#' 
#' all.equal(sample1, sample2) # Only difference: package name
#' }
#'
#' @export getsamples

getsamples <- function(obj, param = NULL, package = NULL, n.samples = 10000,
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
  # (2) Extract or generate sampled data
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
    db <- data.frame("V1" = seq_len(n.samples))
    #
    for (i in seq_len(N)) {
      a <- rnorm(n.samples, yi[i], sei[i])
      db <- cbind(db, a)
    }
    #
    sams <- as.matrix(db[, 1 + seq_len(N)])
    rownames(sams) <- seq_len(n.samples)
    colnames(sams) <- studlab
  }
  else {
    if (missing(sm))
      stop("Argument 'sm' must be provided for R objects created with ",
           "R package ", package, ".")
    #
    if (package == "R2OpenBUGS") {
      sams <- as.data.frame(obj$sims.list[param])
      sams <- as.matrix(sams)
    }
    else if (package == "R2jags") {
      sams <- as.data.frame(obj$BUGSoutput$sims.list[param])
      sams <- as.matrix(sams)
    }
    else if (package == "rjags") {
      N <- as.numeric(nrow(obj[[param]]))
      Np <- N + 1
      sams2 <- data.frame(study = vector("numeric", 0),
                          sams0 = vector("numeric", 0))
      for (i in 1:N) {
        sams0 <- obj[[param]][i, , ]
        NN <- length(sams0)
        sams1 <- data.frame("n_samples" = 1:NN, study = i, sams0)
        sams2 <- rbind(sams2, sams1)
      }
      #
      sams3 <- reshape(data = sams2, direction = "wide",
                       idvar = "n_samples", timevar = "study")
      sams3 <- sams3[, 2:Np]
      sams <- as.matrix(sams3)
    }
    else if (package == "brms") {
      sams_df <- as.data.frame(obj$fit)
      sams_l <- sams_df[, grepl(paste0("r_", param), names(sams_df))]
      sams <- as.matrix(sams_l)
    }
    else
      stop("Please note that the R package for Bayesian analysis ",
           "should either be 'R2OpenBUGS', 'R2jags', 'rjags' or 'brms'")
  }
  #
  if (!transf)
    sams <- transf(sams, sm)
  
  #
  # (3) Return results
  #
  res <- list(data = sams, param = param, package = package,
              n.samples = n.samples,
              sm = sm, transf = transf)
  #
  class(res) <- "samples_metainc"
  #
  res
}


#' @method print samples_metainc
#' @export

print.samples_metainc <- function(x, n = 5, ...) {
  chkclass(x, "samples_metainc")
  #
  print(head(x$data, n), ...)
  #
  invisible(NULL)
}
