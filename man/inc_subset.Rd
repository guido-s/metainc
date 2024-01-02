\name{inc_subset}
\alias{inc_subset}

\title{
  Decision Inconsistency index and Across-Studies Inconsistency index in subsets of studies}

\description{
  Computation of the Decision Inconsistency index and the Across-Studies Inconsistency index for specific subsets of studies (allowing for subgroup analysis).
}

\usage{
  inc_subset(object,data,variable,subset,exc=FALSE)
}

\arguments{
  \item{object}{An object of the \code{inc} class.}
  \item{data}{A dataframe with the same number of primary studies as those included in meta-analysis and containing the variables based on which subset analysis is to be performed}
  \item{variable}{Variable containing the category based on which subset analysis is to be performed}
  \item{subset}{Category of the \code{variable} based on which subset analysis is to be performed}
  \item{exc}{A logical value indicating whether the indicated \code{subset} corresponds to those primary studies set to be excluded or included. If \code{exc = FALSE} (default), the analysis will be performed considering only those studies of the indicated subset. If \code{exc = TRUE}, the analysis will be performed excluding those studies of the indicated subset.}]
}
\details{
  The \code{inc_subset} function computes the Decision Inconsistency index and the Across-Studies Inconsistency index for subsets (subgroups) of studies. This function takes, as its single argument, an object of the class \code{inc}.
}

\value{An object of class \code{inc_subset}, for which some standard methods are available. Some of the components include:
    \item{DI}{A percentage corresponding to the Decision Inconsistency index for the considered subset. The higher/closer to 100\% the value, the higher the inconsistency.}
  \item{ASI}{A percentage corresponding to the Across-Studies Inconsistency index for the considered subset. The higher/closer to 100\% the value, the higher the across-studies inconsistency.}
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
}

\examples{

  ## Example with subset analysis restricted to studies with a low risk of bias:

  data(anticoagulation)
  inc_anticoagulation <-
    inc(log(anticoagulation), t1=20, t2=30, t3=40, br=0.5, sm="or")
  
  
#  data(anticoagulation_df)
#  subset_anticoagulation <-
#    inc_subset(object = inc_anticoagulation, data = anticoagulation_df,
#    variable = "RoB", subset = "Low")
#  subset_anticoagulation
#  
#  
#  ## Example with subset analysis excluding studies with a high risk of bias:
#
#  inc_anticoagulation <-
#    inc(log(anticoagulation), t1=20, t2=30, t3=40, br=0.5, sm="or")
#  
#  data(anticoagulation_df)
#  subset_anticoagulation <-
#    inc_subset(object = inc_anticoagulation, data = anticoagulation_df,
#    variable = "RoB", subset = "High", exc=TRUE)
#  subset_anticoagulation
}
