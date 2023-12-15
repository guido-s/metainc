\name{inc_inf}
\alias{inc_inf}

\title{
Decision Inconsistency index for individual studies and leave-one out sensitivity analysis for inconsistency measures}

\description{
Computation of the Decision Inconsistency index for each individual primary study and leave-one-out sensitivity analysis for the Decision Inconsistency index and the Across-Studies Inconsistency index.
}

\usage{
inc_inf(object)
}

\arguments{
  \item{object}{An object of the \code{inc} class.}
}
\details{
The \code{inc_inf} function (i) computes the Decision Inconsistency index (DI) for each individual primary study and (ii) performs leave-one-out sensitivity analysis for the DI and for the Across-Studies Inconsistency index. This function takes, as its single argument, an object of the class \code{inc}.
}
\value{The \code{inc_inf} function generates:
\item{di_stud}{A dataframe displaying the Decision Inconsistency Index for each individual primary study.}
\item{di_asi_sens}{A dataframe displaying the results of the leave-one-out sensitivity analysis for the Decision Inconsistency index and the Across-Studies Inconsistency index.}
}

\examples{
data(anticoagulation)
inc_anticoagulation <- inc(sims=anticoagulation, t1=20, t2=30, t3=40, br=0.5, sm="or")
sens_anticoagulation <- inc_inf(inc_anticoagulation)
sens_anticoagulation
}
