\name{getsims}
\alias{getsims}

\title{
Obtention of the simulation results of the effect size measures of the primary studies
}

\description{
Obtains a matrix containing the simulation results of the effect size measures of the primary studies.
}

\usage{
getsims(obj, param, package=c("R2OpenBUGS","R2jags","rjags","brms","metafor","meta"),
 nsims=10000, t_exp = FALSE)
}

\arguments{
  \item{obj}{Object obtained after performing meta-analysis using the functions \code{bugs()} [R2OpenBUGS package], \code{jags()} [R2jags package], \code{jags.samples()} [rjags package], \code{brm()} [brms package], \code{rma()} [metafor package] or those of [meta package] (\code{metabin}, \code{metacont} or \code{metagen}).}
  \item{param}{Parameter corresponding to the effect size measures of the primary studies (R2OpenBUGS, R2jags, rjags packages) or to the name of the column that contains the unique identification of each primary study in the dataframe containing the results of primary studies used for meta-analysis (brms package). This argument is not defined if meta-analysis has been performed using metafor or meta packages.}
  \item{package}{Indication of the package that was used to implement meta-analysis in R. It can either be defined as \code{"R2OpenBUGS"}, \code{"R2jags"}, \code{"rjags"}, \code{"brms"}, \code{"metafor"} or \code{"meta"}.}
  \item{nsims}{Number of simulations to be performed in case meta-analysis had been performed with with the \code{rma()} function of the metafor package or with the functions of the meta package. Default is 10000.}
  \item{t_exp}{A logical value indicating whether to perform an exponential transformation on the effect size measures of the primary studies. If \code{t_exp = TRUE}, an exponential transformation will be applied (potentially useful if effect size measures of the primary studies have been computed in the natural logarithm scale). If \code{t_exp = FALSE}, no such transformation will be applied.}
}
\details{
The \code{getsims} function takes as input an object resulting from the implementation of a meta-analytical model in R2OpenBUGS, R2jags, rjags, brms, metafor or meta. This function allows to obtain a matrix containing the simulation results of the effect size measures of the primary studies. This matrix can then be used as an input for the function \code{\link{inc}} for estimation of the Decision Inconsistency index and the Across-Studies Inconsistency index.
}
\value{A matrix containing the simulation results of the effect size measures for each primary study.
}

\examples{
## Example with an object obtained using the R2OpenBUGS package ("es",
## whose simulations on the effect sizes of primary studies are those
## of the "delta" parameter). If the object had been obtained using a
## different package, all remaining arguments would be the same and only
## the "package" argument would have a different input.
## In this example, t_exp=TRUE, as simulation results had been obtained as
## ln(odds ratio).

## In order to run this example, please download the "es.Rdata" file from
## \link{https://github.com/BernardoSousaPinto/metainc_extra_files}

## es_obj <- es
## sims <- getsims(obj=es_obj, param="delta", package="R2OpenBUGS", t_exp=TRUE)


## Example with a dataset containing variables indicating the effect sizes
## for each primary study (yi) and the respective variances (vi), and for
## which frequentist meta-analysis using the metafor package is applied.
## In this example, t_exp=TRUE, as meta-analytical results had been obtained as
## ln(odds ratio).

## data("anticoagulation_df")
## m1 <- rma(anticoagulation_df)
## sims <- getsims(obj=m1, package="metafor", t_exp=TRUE)
}
