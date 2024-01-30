#' Meta-analysis on the association between parenteral anticoagulation
#' and mortality in patients with cancer
#' 
#' @description
#' Data set from meta-analysis on the association between parenteral
#' anticoagulation and mortality in patients with cancer (Akl et al.,
#' 2017). For each primary study, results are expressed as log odds
#' ratio.
#' 
#' @name anticoagulation_df
#' @aliases anticoagulation_df
#' 
#' @docType data
#' 
#' @details
#' 
#' The \code{anticoagulation_df} data frame displays, for each primary study of
#' the meta-analysis from Akl et al. (2017):
#' \enumerate{
#' 
#' \item number of events (variables \code{rA} and \code{rB}) and
#'   total number of participants (\code{nA} and \code{nB}) for each
#'   group,
#' \item effect sizes [ln OR] (\code{yi}),
#' \item respective variance (\code{vi}),
#' \item risk of bias assessment (\code{RoB}).
#' }
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
#' data(anticoagulation_df)
#' m1 <- meta::metagen(yi, sqrt(vi), sm = "OR", data = anticoagulation_df,
#'   studlab = LETTERS[1:18])
#' set.seed(1090) # Make sampled effect sizes reproducible
#' sample1 <- getsamples(m1)
#' dis1 <- inc(sample1, dt1 = 16, dt2 = 31, dt3 = 60, br = 0.504, sm = "OR")
#' 
#' \dontrun{
#' # Same sampled effect sizes using R package metafor (must be installed)
#' m2 <- metafor::rma(anticoagulation_df, measure = "OR", slab = LETTERS[1:18])
#' set.seed(1090) # Make sampled effect sizes reproducible
#' sample2 <- getsamples(m2)
#' dis2 <- inc(sample2, dt1 = 16, dt2 = 31, dt3 = 60, br = 0.504, sm = "OR")
#' 
#' all.equal(sample1, sample2) # Only difference: package name
#' }

NULL
