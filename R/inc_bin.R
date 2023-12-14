#' Title line
#'
#' @description
#' 
#' Brief description (more details could be provided under Details).
#'
#' @param sims R object ...
#' @param baseline.risk ...
#' @param t ...
#' @param t1 ...
#' @param t2 ...
#' @param t3 ...
#' @param sm ... "or", "rr", or "hr"  ...
#' @param scale ...
#' @param x R object of class \code{inc_bin} or  \code{summary.inc_bin}.
#' @param object R object of class \code{inc_bin}.
#' @param \dots Additional arguments (ignored)
#'
#' @details
#' ...
#'
#' @return
#' Description of returned R object (data frame?)
#'
#' @author Bernardo Sousa-Pinto \email{bernardo@@med.up.pt}
#' 
#' #' @examples
#' # add example
#'
#' @export inc_bin
#' 
#' @importFrom stats reshape rnorm
#' @importFrom confintr ci_cramersv

inc_bin <- function(sims, baseline.risk, t, t1, t2, t3,
                    sm, scale = 1000) {
  
  if (missing(baseline.risk)) {
    warning("Please note that judgements based on absolute effects tend to be considered more important for decision-making", call. = F)

    if (missing(t1)|missing(t2)|missing(t3)) {
      if (missing(t)) {
        stop("Please either provide (i) a positive value for argument t or (ii) a positive value for t1, t2 and t3")
      }
      else {
        t <- abs(t)

        sim_stud <- ifelse(sims>t,"higher",ifelse(sims< -t,"lower","trivial"))
        sim_studd <- as.data.frame(sim_stud)
        N <- ncol(sim_studd)

        sim_studa <- reshape(data=sim_studd,direction="long",varying=colnames(sim_studd)[1:N],v.names="effect")
        sim_studa <- sim_studa[,1:2]
        sim_studa1 <- table(sim_studa)
        nc <- dim(sim_studa1)[2]
        if (nc==1) {
          ds <- 0
        }else {
          if (nc==2 & N>2) {
            suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
            v <- ci_cramersv(X2)
            ds <- (sqrt((as.numeric(v$estimate)*as.numeric(v$estimate))/2))*100
          }else {
            suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
            v <- ci_cramersv(X2)
            ds <- (as.numeric(v$estimate))*100
          }}


        simsm <- as.matrix(sim_stud)
        cl00 <- c(simsm)
        cl0 <- data.frame("classification"=cl00)
        cl01 <- with(cl0, tapply(classification, classification, FUN = length))
        cl01 <- data.frame(classification = names(cl01), n = unname(cl01))
        cl_n <- data.frame(classification=c("higher","trivial","lower"),n=c(0,0,0))
        cll <- rbind(cl01,cl_n)
        cl1 <- with(cll, tapply(cll[,2], cll[,1], sum))
        cl1 <- data.frame(classification = names(cl1), n = unname(cl1))
        cl1$proportion <- round((cl1$n)/sum(cl1$n),4)
        cl1$d <- abs((1/3) - cl1$proportion)
        cl2 <- cl1[, c("classification", "proportion")]
        propnull <- round(length(sims[sims>1])/sum(cl1$n),3)
        di <- (1-((0.5*sum(cl1$d))/(2/3)))*100

        }}else {

 if (t3 <= t2|t3 <= t1|t2 <= t1) {
          stop("Please note that you should have t3 > t2 > t1")
  }else {
    t3 <- abs(t3)
    t2 <- abs(t2)
    t1 <- abs(t1)
    sim_stud <- ifelse(sims>=t3,"large (higher)",ifelse(sims>=t2,"moderate (higher)",ifelse(sims>=t1,"small (higher)",ifelse(sims> -t1,"not meaningful",ifelse(sims> -t2,"small (lower)",ifelse(sims> -t3,"moderate (lower)","large (lower)"))))))
    sim_studd <- as.data.frame(sim_stud)
    N <- ncol(sim_studd)

    sim_studa <- reshape(data=sim_studd,direction="long",varying=colnames(sim_studd)[1:N],v.names="effect")
    sim_studa <- sim_studa[,1:2]
    sim_studa1 <- table(sim_studa)
    nc <- dim(sim_studa1)[2]
    if (nc==1) {
      ds <- 0
    }else {
      if (nc<7 & N>=7) {
        suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
        v <- ci_cramersv(X2)
        ds <- (sqrt(((as.numeric(v$estimate))*(as.numeric(v$estimate)))/(1/((nc-1)/6))))*100
      }else {
        if (nc<7 & N<7 & nc<N) {
          suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
          v <- ci_cramersv(X2)
          ds <- (sqrt(((as.numeric(v$estimate))*(as.numeric(v$estimate)))/(1/((nc-1)/(N-1)))))*100
        }else {
          suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
          v <- ci_cramersv(X2)
          ds <- (as.numeric(v$estimate))*100
        }}}

    simsm <- as.matrix(sim_stud)
    cl00 <- c(simsm)
    cl0 <- data.frame("classification"=cl00)
    cl01 <- with(cl0, tapply(classification, classification, FUN = length))
    cl01 <- data.frame(classification = names(cl01), n = unname(cl01))
    cl_n <- data.frame(classification=c("large (higher)","moderate (higher)","small (higher)","not meaningful","small (lower)","moderate (lower)","large (lower)"),n=c(0,0,0,0,0,0,0))
    cll <- rbind(cl01,cl_n)
    cl1 <- with(cll, tapply(cll[,2], cll[,1], sum))
    cl1 <- data.frame(classification = names(cl1), n = unname(cl1))
    cl1$proportion <- round((cl1$n)/sum(cl1$n),4)
    cl1$d <- abs((1/7) - cl1$proportion)
    cl2 <- cl1[, c("classification", "proportion")]
    propnull <- round(length(sims[sims>1])/sum(cl1$n),3)
    di <- (1-((0.5*sum(cl1$d))/(6/7)))*100
  }}
  }else {
    if (baseline.risk<=0 | baseline.risk>1) {
      stop("For the baseline risk (baseline.risk), please provide a value between 0 and 1")
    }else {
      if (missing(sm)) {
        stop("Please provide an effect size measure (sm): sm should be stated as 'or' (odds ratio), 'rr' (risk ratio) or 'hr' (hazard ratio)")
      }else {
        if (sm=="or") {
          if (missing(t1)|missing(t2)|missing(t3)) {
            if (missing(t)) {
              stop("Please either provide (i) a positive value for argument t or (ii) a positive value for t1, t2 and t3")
            }else {
          t <- abs(t)
          bp <- -scale*(baseline.risk-((baseline.risk*sims)/(1-baseline.risk+(baseline.risk*sims))))
          sim_stud <- ifelse(bp>t,"higher",ifelse(bp< -t,"lower","trivial"))
          sim_studd <- as.data.frame(sim_stud)
          N <- ncol(sim_studd)

          sim_studa <- reshape(data=sim_studd,direction="long",varying=colnames(sim_studd)[1:N],v.names="effect")
          sim_studa <- sim_studa[,1:2]
          sim_studa1 <- table(sim_studa)
          nc <- dim(sim_studa1)[2]
          if (nc==1) {
            ds <- 0
          }else {
            if (nc==2 & N>2) {
              suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
              v <- ci_cramersv(X2)
              ds <- (sqrt((as.numeric(v$estimate)*as.numeric(v$estimate))/2))*100
            }else {
              suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
              v <- ci_cramersv(X2)
              ds <- (as.numeric(v$estimate))*100
            }}


          simsm <- as.matrix(sim_stud)
          cl00 <- c(simsm)
          cl0 <- data.frame("classification"=cl00)
          cl01 <- with(cl0, tapply(classification, classification, FUN = length))
          cl01 <- data.frame(classification = names(cl01), n = unname(cl01))
          cl_n <- data.frame(classification=c("higher","trivial","lower"),n=c(0,0,0))
          cll <- rbind(cl01,cl_n)
          cl1 <- with(cll, tapply(cll[,2], cll[,1], sum))
          cl1 <- data.frame(classification = names(cl1), n = unname(cl1))
          cl1$proportion <- round((cl1$n)/sum(cl1$n),4)
          cl1$d <- abs((1/3) - cl1$proportion)
          cl2 <- cl1[, c("classification", "proportion")]
          propnull <- round(length(bp[bp>1])/sum(cl1$n),3)
          di <- (1-((0.5*sum(cl1$d))/(2/3)))*100
      }}else {

         if (t3 <= t2|t3 <= t1|t2 <= t1) {
           stop("Please note that you should have t3 > t2 > t1")
              }else {
          t3 <- abs(t3)
          t2 <- abs(t2)
          t1 <- abs(t1)
          bp <- -scale*(baseline.risk-((baseline.risk*sims)/(1-baseline.risk+(baseline.risk*sims))))
          sim_stud <- ifelse(bp>=t3,"large (higher)",ifelse(bp>=t2,"moderate (higher)",ifelse(bp>=t1,"small (higher)",ifelse(bp> -t1,"not meaningful",ifelse(bp> -t2,"small (lower)",ifelse(bp> -t3,"moderate (lower)","large (lower)"))))))
          sim_studd <- as.data.frame(sim_stud)
          N <- ncol(sim_studd)

          sim_studa <- reshape(data=sim_studd,direction="long",varying=colnames(sim_studd)[1:N],v.names="effect")
          sim_studa <- sim_studa[,1:2]
          sim_studa1 <- table(sim_studa)
          nc <- dim(sim_studa1)[2]
          if (nc==1) {
            ds <- 0
          }else {
            if (nc<7 & N>=7) {
              suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
              v <- ci_cramersv(X2)
              ds <- (sqrt(((as.numeric(v$estimate))*(as.numeric(v$estimate)))/(1/((nc-1)/6))))*100
            }else {
              if (nc<7 & N<7 & nc<N) {
                suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
                v <- ci_cramersv(X2)
                ds <- (sqrt(((as.numeric(v$estimate))*(as.numeric(v$estimate)))/(1/((nc-1)/(N-1)))))*100
              }else {
                suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
                v <- ci_cramersv(X2)
                ds <- (as.numeric(v$estimate))*100
              }}}

          simsm <- as.matrix(sim_stud)
          cl00 <- c(simsm)
          cl0 <- data.frame("classification"=cl00)
          cl01 <- with(cl0, tapply(classification, classification, FUN = length))
          cl01 <- data.frame(classification = names(cl01), n = unname(cl01))
          cl_n <- data.frame(classification=c("large (higher)","moderate (higher)","small (higher)","not meaningful","small (lower)","moderate (lower)","large (lower)"),n=c(0,0,0,0,0,0,0))
          cll <- rbind(cl01,cl_n)
          cl1 <- with(cll, tapply(cll[,2], cll[,1], sum))
          cl1 <- data.frame(classification = names(cl1), n = unname(cl1))
          cl1$proportion <- round((cl1$n)/sum(cl1$n),4)
          cl1$d <- abs((1/7) - cl1$proportion)
          cl2 <- cl1[, c("classification", "proportion")]
          propnull <- round(length(bp[bp>1])/sum(cl1$n),3)
          di <- (1-((0.5*sum(cl1$d))/(6/7)))*100

        }}}else {
          if (sm=="rr") {
            if (missing(t1)|missing(t2)|missing(t3)) {
              if (missing(t)) {
                stop("Please either provide (i) a positive value for argument t or (ii) a positive value for t1, t2 and t3")
              }else {
            t <- abs(t)
            bp <- (baseline.risk*sims - baseline.risk)*scale
            sim_stud <- ifelse(bp>t,"higher",ifelse(bp< -t,"lower","trivial"))
            sim_studd <- as.data.frame(sim_stud)
            N <- ncol(sim_studd)

            sim_studa <- reshape(data=sim_studd,direction="long",varying=colnames(sim_studd)[1:N],v.names="effect")
            sim_studa <- sim_studa[,1:2]
            sim_studa1 <- table(sim_studa)
            nc <- dim(sim_studa1)[2]
            if (nc==1) {
              ds <- 0
            }else {
              if (nc==2 & N>2) {
                suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
                v <- ci_cramersv(X2)
                ds <- (sqrt((as.numeric(v$estimate)*as.numeric(v$estimate))/2))*100
              }else {
                suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
                v <- ci_cramersv(X2)
                ds <- (as.numeric(v$estimate))*100
              }}


            simsm <- as.matrix(sim_stud)
            cl00 <- c(simsm)
            cl0 <- data.frame("classification"=cl00)
            cl01 <- with(cl0, tapply(classification, classification, FUN = length))
            cl01 <- data.frame(classification = names(cl01), n = unname(cl01))
            cl_n <- data.frame(classification=c("higher","trivial","lower"),n=c(0,0,0))
            cll <- rbind(cl01,cl_n)
            cl1 <- with(cll, tapply(cll[,2], cll[,1], sum))
            cl1 <- data.frame(classification = names(cl1), n = unname(cl1))
            cl1$proportion <- round((cl1$n)/sum(cl1$n),4)
            cl1$d <- abs((1/3) - cl1$proportion)
            cl2 <- cl1[, c("classification", "proportion")]
            propnull <- round(length(bp[bp>1])/sum(cl1$n),3)
            di <- (1-((0.5*sum(cl1$d))/(2/3)))*100

          }}else {
                if (t3 <= t2|t3 <= t1|t2 <= t1) {
                  stop("Please note that you should have t3 > t2 > t1")
                }else {
                  t3 <- abs(t3)
                  t2 <- abs(t2)
                  t1 <- abs(t1)
                  bp <- (baseline.risk*sims - baseline.risk)*scale
                  sim_stud <- ifelse(bp>=t3,"large (higher)",ifelse(bp>=t2,"moderate (higher)",ifelse(bp>=t1,"small (higher)",ifelse(bp> -t1,"not meaningful",ifelse(bp> -t2,"small (lower)",ifelse(bp> -t3,"moderate (lower)","large (lower)"))))))
                  sim_studd <- as.data.frame(sim_stud)
                  N <- ncol(sim_studd)

                  sim_studa <- reshape(data=sim_studd,direction="long",varying=colnames(sim_studd)[1:N],v.names="effect")
                  sim_studa <- sim_studa[,1:2]
                  sim_studa1 <- table(sim_studa)
                  nc <- dim(sim_studa1)[2]
                  if (nc==1) {
                    ds <- 0
                  }else {
                    if (nc<7 & N>=7) {
                      suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
                      v <- ci_cramersv(X2)
                      ds <- (sqrt(((as.numeric(v$estimate))*(as.numeric(v$estimate)))/(1/((nc-1)/6))))*100
                    }else {
                      if (nc<7 & N<7 & nc<N) {
                        suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
                        v <- ci_cramersv(X2)
                        ds <- (sqrt(((as.numeric(v$estimate))*(as.numeric(v$estimate)))/(1/((nc-1)/(N-1)))))*100
                      }else {
                        suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
                        v <- ci_cramersv(X2)
                        ds <- (as.numeric(v$estimate))*100
                      }}}

                  simsm <- as.matrix(sim_stud)
                  cl00 <- c(simsm)
                  cl0 <- data.frame("classification"=cl00)
                  cl01 <- with(cl0, tapply(classification, classification, FUN = length))
                  cl01 <- data.frame(classification = names(cl01), n = unname(cl01))
                  cl_n <- data.frame(classification=c("large (higher)","moderate (higher)","small (higher)","not meaningful","small (lower)","moderate (lower)","large (lower)"),n=c(0,0,0,0,0,0,0))
                  cll <- rbind(cl01,cl_n)
                  cl1 <- with(cll, tapply(cll[,2], cll[,1], sum))
                  cl1 <- data.frame(classification = names(cl1), n = unname(cl1))
                  cl1$proportion <- round((cl1$n)/sum(cl1$n),4)
                  cl1$d <- abs((1/7) - cl1$proportion)
                  cl2 <- cl1[, c("classification", "proportion")]
                  propnull <- round(length(bp[bp>1])/sum(cl1$n),3)
                  di <- (1-((0.5*sum(cl1$d))/(6/7)))*100

          }}}else {
            if (sm=="hr") {
              if (missing(t1)|missing(t2)|missing(t3)) {
                if (missing(t)) {
                  stop("Please either provide (i) a positive value for argument t or (ii) a positive value for t1, t2 and t3")
                }else {
              t <- abs(t)
              bp <- -(baseline.risk-(1-(1-baseline.risk)^sims))*scale
              sim_stud <- ifelse(bp>t,"higher",ifelse(bp< -t,"lower","trivial"))
              sim_studd <- as.data.frame(sim_stud)
              N <- ncol(sim_studd)

              sim_studa <- reshape(data=sim_studd,direction="long",varying=colnames(sim_studd)[1:N],v.names="effect")
              sim_studa <- sim_studa[,1:2]
              sim_studa1 <- table(sim_studa)
              nc <- dim(sim_studa1)[2]
              if (nc==1) {
                ds <- 0
              }else {
                if (nc==2 & N>2) {
                  suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
                  v <- ci_cramersv(X2)
                  ds <- (sqrt((as.numeric(v$estimate)*as.numeric(v$estimate))/2))*100
                }else {
                  suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
                  v <- ci_cramersv(X2)
                  ds <- (as.numeric(v$estimate))*100
                }}


              simsm <- as.matrix(sim_stud)
              cl00 <- c(simsm)
              cl0 <- data.frame("classification"=cl00)
              cl01 <- with(cl0, tapply(classification, classification, FUN = length))
              cl01 <- data.frame(classification = names(cl01), n = unname(cl01))
              cl_n <- data.frame(classification=c("higher","trivial","lower"),n=c(0,0,0))
              cll <- rbind(cl01,cl_n)
              cl1 <- with(cll, tapply(cll[,2], cll[,1], sum))
              cl1 <- data.frame(classification = names(cl1), n = unname(cl1))
              cl1$proportion <- round((cl1$n)/sum(cl1$n),4)
              cl1$d <- abs((1/3) - cl1$proportion)
              cl2 <- cl1[, c("classification", "proportion")]
              propnull <- round(length(bp[bp>1])/sum(cl1$n),3)
              di <- (1-((0.5*sum(cl1$d))/(2/3)))*100
                }}else {
                  if (t3 <= t2|t3 <= t1|t2 <= t1) {
                    stop("Please note that you should have t3 > t2 > t1")
                  }else {
                    t3 <- abs(t3)
                    t2 <- abs(t2)
                    t1 <- abs(t1)
                    bp <- -(baseline.risk-(1-(1-baseline.risk)^sims))*scale
                    sim_stud <- ifelse(bp>=t3,"large (higher)",ifelse(bp>=t2,"moderate (higher)",ifelse(bp>=t1,"small (higher)",ifelse(bp> -t1,"not meaningful",ifelse(bp> -t2,"small (lower)",ifelse(bp> -t3,"moderate (lower)","large (lower)"))))))
                    sim_studd <- as.data.frame(sim_stud)
                    N <- ncol(sim_studd)

                    sim_studa <- reshape(data=sim_studd,direction="long",varying=colnames(sim_studd)[1:N],v.names="effect")
                    sim_studa <- sim_studa[,1:2]
                    sim_studa1 <- table(sim_studa)
                    nc <- dim(sim_studa1)[2]
                    if (nc==1) {
                      ds <- 0
                    }else {
                      if (nc<7 & N>=7) {
                        suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
                        v <- ci_cramersv(X2)
                        ds <- (sqrt(((as.numeric(v$estimate))*(as.numeric(v$estimate)))/(1/((nc-1)/6))))*100
                      }else {
                        if (nc<7 & N<7 & nc<N) {
                          suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
                          v <- ci_cramersv(X2)
                          ds <- (sqrt(((as.numeric(v$estimate))*(as.numeric(v$estimate)))/(1/((nc-1)/(N-1)))))*100
                        }else {
                          suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
                          v <- ci_cramersv(X2)
                          ds <- (as.numeric(v$estimate))*100
                        }}}

                    simsm <- as.matrix(sim_stud)
                    cl00 <- c(simsm)
                    cl0 <- data.frame("classification"=cl00)
                    cl01 <- with(cl0, tapply(classification, classification, FUN = length))
                    cl01 <- data.frame(classification = names(cl01), n = unname(cl01))
                    cl_n <- data.frame(classification=c("large (higher)","moderate (higher)","small (higher)","not meaningful","small (lower)","moderate (lower)","large (lower)"),n=c(0,0,0,0,0,0,0))
                    cll <- rbind(cl01,cl_n)
                    cl1 <- with(cll, tapply(cll[,2], cll[,1], sum))
                    cl1 <- data.frame(classification = names(cl1), n = unname(cl1))
                    cl1$proportion <- round((cl1$n)/sum(cl1$n),4)
                    cl1$d <- abs((1/7) - cl1$proportion)
                    cl2 <- cl1[, c("classification", "proportion")]
                    propnull <- round(length(bp[bp>1])/sum(cl1$n),3)
                    di <- (1-((0.5*sum(cl1$d))/(6/7)))*100
            }}}else {
              stop("Please provide a valid effect size measure (sm): sm should be stated as 'or' (odds ratio), 'rr' (risk ratio) or 'hr' (hazard ratio)")
            }
          }
        }
      }}}
  #
  output <- list(ASI = round(ds,3), DI = round(di,3), class_distribution = cl2,prop_over_null = propnull,stud_class=sim_studd, call = match.call())
  class(output) <- "inc_bin"
  output
}





#' @rdname inc_bin
#' @keywords print
#' @method print inc_bin
#' @export 

print.inc_bin <- function (x,...) {
  cat("Decision Inconsistency index (DI): ", round(x$DI,1), "%")
  cat("\n")
  cat("Across-Studies Inconsistency (ASI): ", round(x$ASI,1), "%")
  cat("\n")
  cat("\n")
  cat("Proportion of simulations by classification in relation to the decision threshold:", "\n")
  print(x$class_distribution)
  cat("\n")
  cat("\n")
  cat("Proportion of simulations higher than the null effect: ", round(x$prop_over_null,3))
  cat("\n")
  #
  invisible(NULL)
}





#' @rdname inc_bin
#' @method summary inc_bin
#' @export 

summary.inc_bin <- function(object, ...) {
  output <- object
  class(output) <- "summary.inc_bin"
  output
}





#' @rdname inc_bin
#' @keywords print
#' @method print summary.inc_bin
#' @export 

print.summary.inc_bin <- function (x,...) {
  cat("Decision Inconsistency index (DI): ", round(x$DI,1), "%")
  cat("\n")
  cat("Across-Studies Inconsistency (ASI): ", round(x$ASI,1), "%")
  cat("\n")
  cat("\n")
  cat("Proportion of simulations by classification in relation to the decision threshold:", "\n")
  print(x$class_distribution)
  cat("\n")
  cat("\n")
  cat("Proportion of simulations higher than the null effect: ", round(x$prop_over_null,3))
  cat("\n")
  #
  invisible(NULL)
}
