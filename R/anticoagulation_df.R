#' Data set with data from Akl et al. (2017)
#' 
#' @description
#' Data set with data from Akl et al. (2017)
#' 
#' @name anticoagulation_df
#' @aliases anticoagulation_df
#' 
#' @docType data
#' 
#' @details
#' 
#' The \code{anticoagulation_df} dataframe displays, for each primary study of
#' the meta-analysis from Akl et al. (2017):
#' (i) the number of events and the total number of participants for each
#' group,
#' (ii) its effect size measure [ln OR] (\code{yi} variable),
#' (iii) the respective variance (\code{vi} variable),
#' and (iv) the risk of bias (\code{RoB} variable).
#' 
#' @source
#' Akl, E. A., Kahale, L. A., Hakoum, M. B., Matar, C. F., Sperati, F.,
#' Barba, M., et al. (2017).
#' \dQuote{Parenteral anticoagulation in ambulatory patients with cancer.}
#' \emph{Cochrane Database of Systematic Reviews},
#' \bold{9}: CD006652.
#'
#' @keywords datasets
#' 
#' @seealso \code{\link{inc}}
#'
#' @examples
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

NULL
