\name{metainc-data}
\docType{data}
\alias{metainc-data}
\alias{anticoagulation}
\alias{montelukast}
\alias{anticoagulation_df}

\title{Inconsistency assessment data}

\description{Two matrices with the simulation of effect sizes of primary studies (obtained following a Bayesian meta-analysis) (anticoagulation and montelukast) and one dataframe with the effect size measures and respective variances for each primary study (anticoagulation_df)}


\details{
In each matrix, each column corresponds to a different primary study. The \code{anticoagulation} data is from Akl et al. (2017) and displays results presented as odds ratio (Analysis 1.1.). The \code{montelukast} data is from  Krishnamoorthy et al. (2020) (Figure 7) and displays results presented as mean differences.
The \code{anticoagulation_df} dataframe displays, for each primary study of the meta-analysis from Akl et al. (2017): (i) the number of events and the total number of participants for each group, (ii) its effect size measure [ln OR] (\code{yi} variable), (iii) the respective variance (\code{vi} variable), and (iv) the risk of bias (\code{RoB} variable).

To have access to each of the matrices or dataframe, they should be loaded as follows:

\code{data("anticoagulation")}
\code{data("montelukast")}
\code{data("anticoagulation_df")}

}

\source{
Akl, E. A., Kahale, L. A., Hakoum, M. B., Matar, C. F., Sperati, F., Barba, M., et al. (2017).\dQuote{Parenteral anticoagulation in ambulatory patients with cancer.} \emph{Cochrane Database of Systematic Reviews},
\bold{9}: CD006652.

Krishnamoorthy, M., Mohd Noor, N., Mat Lazim, N., & Abdullah, B. (2020). \dQuote{Efficacy of montelukast in allergic rhinitis treatment: a systematic review and meta-analysis.} \emph{Drugs}, \bold{80}, 1381--1851.
}
