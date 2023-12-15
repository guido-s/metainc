\name{inc}
\alias{inc}

\title{
  Decision Inconsistency index and Across-Studies Inconsistency index
}

\description{
  Calculates the Decision Inconsistency index and the Across-Studies Inconsistency index. 
}

\usage{
  inc(sims,t1,t2,t3,t,sm=c("or","rr","hr","md","smd","rd","gen_dif","gen_ratio"),br,den=1000)
}

\arguments{
  \item{sims}{A matrix containing the simulations of the effect size measures of primary studies.}
  \item{t1}{Numeric. Defines the decision threshold differentiating small from trivial effects.}
  \item{t2}{Numeric. Defines the decision threshold differentiating moderate from small effects.}
  \item{t3}{Numeric. Defines the decision threshold differentiating large from moderate effects.}
  \item{t}{Numeric. Defines the decision threshold distinguishing meaningful from trivial effects. Only needs to be defined if t1, t2 and t3 are not defined.}
  \item{sm}{A character string indicating the effect size measure computed for the primary studies. It can be stated as \code{"or"} (odds ratio),\code{"rr"} (risk ratio), \code{"hr"} (hazard ratio), \code{"md"} (mean difference), \code{"smd"} (standardised mean difference), \code{"rd"} (risk difference), \code{"gen_dif"} (generic difference) or \code{"gen_ratio"} (generic ratio).}
  \item{br}{Numeric. An optional value that can be provided when sm = "or", "rr" or "hr", and indicating the baseline risk (also known as assumed comparator risk) – i.e., the risk that the outcome of interest would occur in the comparison intervention. Must be a numeric value between 0 and 1.}
  \item{den}{Numeric. The number of people per which decision thresholds are provided . Default is 1000 (i.e., decision threshold values are presented per 1000 people). Only needs to be defined if \code{br} is not missing.}
}
\details{
  The \code{inc} function calculates the Decision Inconsistency index (DI) and the Across-Studies Inconsistency index (ASI) for meta-analyses. The following possibilities are considered depending on the type of effect size measures:
    \itemize{
      \item Effect size measure corresponding to a ratio (\code{sm} = \code{"or"}, \code{"rr"} or \code{"hr"}) with the DI and the ASI being calculated based on absolute effects: This requires the setting of a baseline risk (i.e., \code{br} must be defined). The decision threshold values (\code{t1}, \code{t2} and \code{t3} or \code{t}) should be provided as absolute effects (i.e., number of additional or diminished events per N people. By default, it is assumed that these threshold values are provided per 1000 people. However, this can be changed using the \code{den} argument).
      \item Effect size measure corresponding to a ratio (\code{sm} = \code{"or"}, \code{"rr"}, \code{"hr"} or \code{"gen_ratio"}) with the DI and the ASI being calculated based on relative effect size measures: The simulation results of relative effect size measures of primary studies are directly compared with decision thresholds (\code{t1}, \code{t2} and \code{t3} or \code{t}) also expressed as relative effect size measures. This is the adopted approach when no information is provided on the baseline risk (\code{br}).
      \item Effect size measure corresponding to a difference (\code{sm} = \code{"md"}, \code{"smd"}, \code{"rd"} or \code{"gen_dif"}): The simulation results of the effect size measures of primary studies are directly compared with decision thresholds (\code{t1}, \code{t2} and \code{t3} or \code{t}) also expressed as differences.
    }
  Of note, when dealing with relative effect size measures, judgements based on absolute effects tend to be considered more important for decision-making. The formulae for calculating absolute effects based on relative effect size measures are those used in the GRADE approach (see References below).
}
\value{An object of class \code{inc}, for which some standard methods are available, see \code{\link{inc-class}}. Some of the components include:
    \item{DI}{A percentage corresponding to the Decision Inconsistency index. The higher/closer to 100\% the value, the higher the inconsistency.}
  \item{ASI}{A percentage corresponding to the Across-Studies Inconsistency index. The higher/closer to 100\% the value, the higher the across-studies inconsistency.}
  \item{class_distribution}{A dataframe displaying the proportion of simulations indicating:
      \itemize{
        \item Large positive effects (efect sizes higher than \code{tl}): “large (higher)” row;
        \item Moderate positive effects (efect sizes between \code{tm} and \code{tl}): “moderate (higher)” row;
        \item Small positive effects (efect sizes between \code{ts} and \code{tm}): “small (higher)” row;
        \item Non meaningful effects (efect sizes between \code{-ts} and \code{ts}): “not meaningful” row;
        \item Small negative effects (efect sizes between \code{-ts} and \code{-tm}): “small (lower)” row;
        \item Moderate negative effects (efect sizes between \code{-tm} and \code{-tl}): “moderate (lower)” row;
        \item Large negative effects (efect sizes lower than \code{-tl}): “large (lower)” row.
      }}
  \item{prop_over_null}{A numeric value indicating the proportion of simulations with a value higher than than the value representing no difference between the groups.}
}

\references{
  Cohen J. (1998). \dQuote{Statistical Power Analysis in the Behavioral Sciences.}, 2nd edition ed. Hillsdale (NJ): Lawrence Erlbaum Associates, Inc.; 1988.
  
  Schunemann HJ, Higgins JPT, Vist GE, et al. (2019).
  \dQuote{Completing ‘Summary of findings’ tables and grading the certainty of the evidence.} \emph{Cochrane handbook for systematic reviews of interventions}, 375--402.
  
  Schunemann HJ, Vist GE, Higgins JPT, et al. (2019).
\dQuote{Interpreting results and drawing conclusions.} \emph{Cochrane handbook for systematic reviews of interventions}, 403--431.
  
  Skoetz N, Goldkuhle M, van Dalen EC, et al. (2020). \dQuote{GRADE guidelines 27: how to calculate absolute effects for time-to-event outcomes in summary of findings tables and Evidence Profiles.} \emph{Journal of Clinical Epidemiology}, \bold{118}, 124--131.
}

\examples{
  ## Example with effect sizes measures expressed as ratios and with calculation 
  ## of the Decision Inconsistency index and the Across-Studies
  ## Inconsistency index being calculated based on absolute effects:
  
  data(anticoagulation)
  inc_anticoagulation <- inc(sims=anticoagulation, t1=20, t2=30, t3=40, br=0.5, sm="or")
  inc_anticoagulation
  
  
  ## Example with effect sizes measures expressed as ratios and with calculation 
  ## of the Decision Inconsistency index and the Across-Studies
  ## Inconsistency index being calculated directly based on relative
  ## effect size measures (and with only one decision threshold being defined):
  
  data(anticoagulation)
  inc_anticoagulation <- inc(sims=anticoagulation, t=1.05, sm="or")
  inc_anticoagulation
  
  
  ## Example with calculation of the Decision Inconsistency index and the 
  ## Across-Studies Inconsistency index based on effect size measures expressed
  ## as differences:

data(montelukast)
inc_montelukast <- inc_md(sims=montelukast, t1=0.2, t2=0.4, t3=0.6)
inc_montelukast
}
