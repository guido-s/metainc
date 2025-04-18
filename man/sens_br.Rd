% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sens_br.R
\name{sens_br}
\alias{sens_br}
\alias{plot.sens_br}
\title{Sensitivity analysis (based on the baseline risk) for the Decision
Inconsistency and Across-Studies Inconsistency index}
\usage{
sens_br(x, br1, br2, dt1, dt2 = NULL, dt3 = NULL, sm, by = 0.01, scale = 1000)

\method{plot}{sens_br}(
  x,
  ylim1 = c(0, 100),
  ylim2 = c(0, 100),
  ylab1 = "DI index (\%)",
  ylab2 = "ASI index (\%)",
  ...
)
}
\arguments{
\item{x}{An R object created with \code{\link{getsamples}} or a matrix
containing sampled effect sizes of primary studies. Note, log
transformed effect sizes must be provided (e.g., log odds ratios
instead of odds ratios).}

\item{br1}{Smallest baseline risk considered.}

\item{br2}{Largest baseline risk considered.}

\item{dt1}{A single numeric defining the decision threshold to
distinguish (i) meaningful from trivial effects, if arguments
\code{dt2} and \code{dt3} are not provided, or (ii) small from
trivial effects if arguments \code{dt2} and \code{dt3} are
provided.}

\item{dt2}{A single numeric defining the decision threshold to
distinguish moderate from small effects provided.}

\item{dt3}{A single numeric defining the decision threshold to
distinguish large from moderate effects.}

\item{sm}{A character string indicating the summary measure used in
primary studies (either \code{sm = "OR"}, \code{sm = "RR"} or
\code{sm = "HR"}).}

\item{by}{Increment of the sequence from \code{br1} to \code{br2}.}

\item{scale}{The number of people per which absolute decision
thresholds are provided (default: 1000, i.e., absolute decision
threshold values are defined per 1000 people).}

\item{ylim1}{The y limits (min, max) of the plot showing the
Decision Inconsistency index.}

\item{ylim2}{The y limits (min, max) of the plot showing the
Across-Studies Inconsistency index.}

\item{ylab1}{A label for the y-axis (Decision Inconsistency index).}

\item{ylab2}{A label for the y-axis (Across-Studies Inconsistency
index).}

\item{\dots}{Additional graphical arguments (ignored).}
}
\value{
A data frame containing
\item{br}{Baseline risk}
\item{ASI}{Decision Inconsistency index at baseline risk}
\item{DI}{Across-Studies Inconsistency index at baseline risk}
}
\description{
Sensitivity analysis on the Decision Inconsistency index and the
Across-Studies Inconsistency index based on a range of baseline
risks. It is applicable only to meta-analyses with a binary outcome
(effect size measures expressed as risk ratios, odds ratios or
hazard ratios).
}
\details{
Computes the Decision Inconsistency index (DI) and the
Across-Studies Inconsistency index (ASI) across a range of baseline risks.
It can only be applied for meta-analyses with binary outcome data (effect
size measures expressed as (log) risk ratios, odds ratios or hazard ratios),
with the DI and the ASI being calculated based on absolute effects.
As a result, the decision threshold values (\code{dt1}, \code{dt2},
\code{dt3}) must be provided as absolute effects. By default, it is assumed
that threshold values are provided as numbers of events per 1000
persons (\code{scale = 1000}).
}
\examples{
\donttest{
data(anticoagulation)
dis <- sens_br(log(anticoagulation),  br1 = 0.3, br2 = 0.7, dt1 = 20,
  sm = "OR", by = 0.1)
dis
plot(dis, ylim1 = c(0, 100), ylim2 = c(0, 50))
}

}
\references{
Schunemann HJ, Higgins JPT, Vist GE, et al. (2019).
\dQuote{Completing ‘Summary of findings’ tables and grading the certainty of
the evidence.}
\emph{Cochrane Handbook for Systematic Reviews of Interventions},
375--402.

Skoetz N, Goldkuhle M, van Dalen EC, et al. (2020).
\dQuote{GRADE guidelines 27: how to calculate absolute effects for
time-to-event outcomes in summary of findings tables and Evidence Profiles.}
\emph{Journal of Clinical Epidemiology}, \bold{118},
124--131.

Sousa-Pinto B, Neumann I, Vieira RJ, et al. (2025).
\dQuote{Quantitative assessment of inconsistency in meta-analysis using
decision thresholds with two new indices.}
\emph{Journal of Clinical Epidemiology}, \bold{181},
111725.
}
\author{
Bernardo Sousa-Pinto \email{bernardo@med.up.pt},
  Guido Schwarzer \email{guido.schwarzer@uniklinik-freiburg.de}
}
