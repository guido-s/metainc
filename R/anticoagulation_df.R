#' anticoagulation_df
#' 
#' @description
#' Meta-analysis dataset anticoagulation_df
#' 
#' @name anticoagulation_df
#' @aliases anticoagulation_df
#' 
#' @docType data
#' 
#' @format A dataset with 17 studies and 7 variables.
#' 
#' @source
#' Add reference to source (if available).
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
