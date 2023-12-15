\name{inc_br_sens}
\alias{inc_br_sens}

\title{
Sensitivity analysis (based on the baseline risk) for the Decision Inconsistency index and Across-Studies Inconsistency index
}

\description{
Sensitivity analysis on the Decision Inconsistency index and the Across-Studies Inconsistency index based on a range of baseline risks. It is applicable only to meta-analyses with binary outcome data (effect size measures expressed as risk ratios, odds ratios or hazard ratios).
}

\usage{
inc_br_sens <- function(sims,br1,br2,t,sm=c("or","rr","hr"),den=1000,y_min=0,y_max=100)
}

\arguments{
  \item{sims}{A matrix containing the simulations of the effect size measures of primary studies (expressed as risk ratios, odds ratios or hazard ratios).}
  \item{br1}{Numeric. Defines the lowest value to be considered for the baseline risk (also known as the assumed comparator risk) – i.e., the lowest plausible risk that the outcome of interest would occur in the comparison intervention. Must be a numeric value between 0 and 1.}
  \item{br2}{Numeric. Defines the highest value to be considered for the baseline risk (also known as the assumed comparator risk) – i.e., the highest plausible risk that the outcome of interest would occur in the comparison intervention. Must be a numeric value between 0 and 1.}
  \item{t}{Numeric. Defines the decision threshold distinguishing meaningful from trivial effects.}
  \item{sm}{A character string indicating the effect size measure computed for the primary studies. It can be stated as \code{"or"} (odds ratio),\code{"rr"} (risk ratio) or \code{"hr"} (hazard ratio).}
  \item{den}{Numeric. The number of people per which decision thresholds are provided . Default is 1000 (i.e., decision threshold values are presented per 1000 people).}
  \item{y_min}{Numeric. Defines the lowest value for the y axis for the Decision Inconsistency index and the across-studies inconsistency plots. Default is 0.}
  \item{y_max}{Numeric. Defines the highest value for the y axis for the Decision Inconsistency index and the across-studies inconsistency plots. Default is 100.}
}
\details{
The \code{inc_bin_sens} function computes the Decision Inconsistency index (DI) and the Across-Studies Inconsistency index (ASI) across a range of defined baseline risks. This function, therefore, allows for performing sensitivity analysis according to the baseline risk. It can only be applied for meta-analyses with binary outcome data (effect size measures expressed as risk ratios, odds ratios or hazard ratios), with the DI and the ASI being calculated based on absolute effects. As a result, the decision threshold value (\code{t}) must be provided as absolute effects (i.e., the minimum number of additional or diminished events that should occur for the effect to be considered appreciable/important instead of trivial). By default, it is assumed that this threshold value is provided per 1000 people. However, this can be changed using the \code{den} argument.
}
\value{The \code{inc_br_sens} function generates:
\itemize{
\item Plot displaying the Decision Inconsistency index according to the baseline risk;
\item Plot displaying the Across-Studies Inconsistency index according to the baseline risk;
\item	Plot displaying the proportion of simulations (according to the baseline risk) (i) higher than the highest decision threshold value (red), (ii) lower than the lowest decision threshold value (green), and (iii) with values within the range of the decision threshold values (grey);
\item \code{inc_sens}{ A dataframe displaying the Decision Inconsistency index and the Across-Studies Inconsistency index according to the baseline risk.}
}
}

\references{
Schunemann HJ, Higgins JPT, Vist GE, et al. (2019).
\dQuote{Completing ‘Summary of findings’ tables and grading the certainty of the evidence.} \emph{Cochrane handbook for systematic reviews of interventions}, 375--402.

Skoetz N, Goldkuhle M, van Dalen EC, et al. (2020). \dQuote{GRADE guidelines 27: how to calculate absolute effects for time-to-event outcomes in summary of findings tables and Evidence Profiles.} \emph{Journal of Clinical Epidemiology}, \bold{118}, 124--131.
}

\examples{
data(anticoagulation)
dis_anticoagulation <- inc_br_sens(sims=anticoagulation, br1=0.3, br2=0.7, t=20, sm="or")
dis_anticoagulation
}
