\name{inc_t_sens}
\alias{inc_t_sens}

\title{
Sensitivity analysis (based on decision thresholds) for the Decision Inconsistency index and Across-Studies Inconsistency index
}

\description{
Sensitivity analysis on the Decision Inconsistency index and the Across-Studies Inconsistency index based on a range of decision thresholds.
}

\usage{
inc_t_sens(sims, br, t1_min, t1_max, t2_min, t2_max, sm=c("or","rr","hr"), den=1000,
 unit=1, null_effect=0)
}

\arguments{
  \item{sims}{A matrix containing the simulations of the effect size measures of primary studies.}
  \item{br}{Numeric. An optional value indicating the baseline risk (also known as the assumed comparator risk) – i.e., the risk that the outcome of interest would occur in the comparison intervention. It shall be expressed between 0 and 1.}
  \item{t1_min}{Numeric. Defines the lowest value for the lower decision threshold. Must be lower than \code{t1_max}, \code{t2_min}, \code{t2_max} and \code{null_effect}.}
  \item{t1_max}{Numeric. Defines the highest value for the lower decision threshold. Must be higher than \code{t1_min}, but lower than \code{t2_min}, \code{t2_max} and \code{null_effect}.}
  \item{t2_min}{Numeric. Defines the lowest value for the higher decision threshold. Must be higher than \code{t1_min}, \code{t1_max} and \code{null_effect}, but lower than \code{t2_max}.}
  \item{t2_max}{Numeric. Defines the highest value for the higher decision threshold. Must  be higher than \code{t1_min}, \code{t1_max}, \code{t2_min} and \code{null_effect}.}
  \item{sm}{A character string indicating the effect size measure computed for the primary studies. It can be stated as \code{"or"} (odds ratio),\code{"rr"} (risk ratio) or \code{"hr"} (hazard ratio). Only required if \code{br} is not missing.}
  \item{den}{Numeric. The number of people per which decision thresholds are provided. Default is 1000 (i.e., decision threshold values are presented per 1000 people). Only required if \code{br} is not missing.}
  \item{unit}{Numeric. Increment of the sequence between t1_min and t1_max and between t2_min and t2_max, based on which sensitivity analysis should be provided. Default is 1 (i.e., a sequence between threshold values is generated with each of the values having one unit of difference from the value immediately higher or lower). Only required if \code{br} is not missing.}
  \item{null_effect}{Numeric. Defines the value representing no difference between the groups in relation to the effect size measures. Default value is 0 (usual null effect when dealing with absolute effects; when directly dealing with relative effects, the null effect is typically 1).}
}
\details{
The \code{inc_t_sens} function computes the Decision Inconsistency index (DI) and the Across-Studies Inconsistency (ASI) index across a range of defined decision thresholds. This function, therefore, allows for performing sensitivity analysis according to the decision thresholds. Two possibilities are considered:
\itemize{
\item The DI and the ASI are calculated based on absolute effects for meta-analyses of odds ratios, risk ratios or hazard ratios: This requires the setting of a baseline risk (i.e., \code{br} must be defined), as well as the definition of the effect size measure computed for the primary studies (i.e., \code{sm} must be defined either as “odds ratio” [\code{"or"}], “risk ratio” [\code{"rr"}] or “hazard ratio” [\code{"hr"}]). The decision threshold values (\code{t1_min}, \code{t1_max}, \code{t2_min}, \code{t2_max}) should be provided as absolute effects (i.e., minimum number of additional or diminished events that should occur for the effect to be considered appreciable/important instead of trivial). By default, it is assumed that these threshold values are provided per 1000 people. However, this can be changed using the \code{den} argument.
\item The DI and the ASI are calculated by directly comparing the effect size measures of primary studies with decision thresholds (\code{t1_min}, \code{t1_max}, \code{t2_min}, \code{t2_max}). This is the most appropriate approach except when assessing absolute effects for meta-analysis of odds ratios, risk ratios or hazard ratios. The \code{br}, \code{sm} and \code{den} arguments should not be defined.
}
}
\value{The \code{inc_t_sens} function generates:
\itemize{
\item A heatmap plot displaying the Decision Inconsistency index according to the different combinations of decision thresholds;
\item A heatmap plot displaying the Across-studies Inconsistency index according to the different combinations of decision thresholds;
\item	Heatmap plots displaying the proportion of simulations (according to the different combinations of decision thresholds) (i) higher than the highest decision threshold value, (ii) lower than the lowest decision threshold value, and (iii) with values within the range of the decision threshold values;
\item \code{inc_sens}{ A dataframe displaying the Decision Inconsistency index and the Across-Studies Inconsistency index according to different combinations of decision thresholds.}
}
}

\references{
Schunemann HJ, Higgins JPT, Vist GE, et al. (2019).
\dQuote{Completing ‘Summary of findings’ tables and grading the certainty of the evidence.} \emph{Cochrane handbook for systematic reviews of interventions}, 375--402.

Skoetz N, Goldkuhle M, van Dalen EC, et al. (2020). \dQuote{GRADE guidelines 27: how to calculate absolute effects for time-to-event outcomes in summary of findings tables and Evidence Profiles.} \emph{Journal of Clinical Epidemiology}, \bold{118}, 124--131.
}

\examples{
data(anticoagulation)
dis_anticoagulation <- inc_t_sens(sims=anticoagulation, br=0.5, t1_min=-30, t1_max=-10,
 t2_min=10, t2_max=30, sm="or")
dis_anticoagulation
}
