#' Decision Inconsistency index for individual studies and leave-one out sensitivity analysis for inconsistency measures
#'
#' @description
#' Computation of the Decision Inconsistency index for each individual primary study and leave-one-out sensitivity analysis for the Decision Inconsistency index and the Across-Studies Inconsistency index.
#'
#' @param x An R object of class \code{inc}.
#' @param object An R object of class \code{inf}.
#' @param \dots Additional arguments (ignored).
#'
#' @details
#' The \code{inf} function (i) computes the Decision Inconsistency index (DI)
#' for each individual primary study and (ii) performs leave-one-out
#' sensitivity analysis for the DI and for the Across-Studies Inconsistency
#' index. This function takes, as its single argument, an object of the class
#' \code{inc}.
#' 
#' @return
#' The \code{inf} function generates:
#' \item{di_stud}{A dataframe displaying the Decision Inconsistency Index for each individual primary study.}
#' \item{di_asi_sens}{A dataframe displaying the results of the leave-one-out sensitivity analysis for the Decision Inconsistency index and the Across-Studies Inconsistency index.}
#' 
#' @author Bernardo Sousa-Pinto \email{bernardo@@med.up.pt},
#'   Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @examples
#' data(anticoagulation)
#' inc_anticoagulation <-
#'   inc(log(anticoagulation), t1 = 20, t2 = 30, t3 = 40, br = 0.5, sm = "or")
#' sens <- inf(inc_anticoagulation)
#' sens
#'
#' @export inf

inf <- function(x) {
  
  #
  # (1) Check arguments
  #
  if (!inherits(x, "inc"))
    stop("Argument 'x' must be of class 'inc'.", call. = FALSE)
  #
  object <- x
  
  #
  # (2) Calculate DIS and AI indices
  #
  x <- as.data.frame(lapply(object$stud_class, as.character))
  N <- ncol(x)
  
  simsm <- as.matrix(x)
  
  if (is.na(object$t2)) {
    di_stud <- data.frame("study"=c(),"di"=c())
    
    for(i in 1:N) {
      cl00 <- c(simsm[, i])
      cl0 <- data.frame("classification"=cl00)
      cl01 <- with(cl0, tapply(classification, classification, FUN = length))
      cl01 <- data.frame(classification = names(cl01), n = unname(cl01))
      cl_n <- data.frame(classification=c("higher","trivial","lower"),n=c(0,0,0))
      cll <- rbind(cl01,cl_n)
      cl1 <- with(cll, tapply(cll[,2], cll[,1], sum))
      cl1 <- data.frame(classification = names(cl1), n = unname(cl1))
      cl1$proportion <- round((cl1$n)/sum(cl1$n),4)
      cl1$d <- abs((1/3) - cl1$proportion)
      di <- round((1-((0.5*sum(cl1$d))/(2/3)))*100,1)
      dii <- data.frame("Study"=i,"DI"=di)
      di_stud <- rbind(di_stud, dii)
    }
    
    di_sens <- data.frame("study"=c(),"di"=c())
    
    for(i in 1:N) {
      cl00 <- as.matrix(simsm[,-i])
      cl00 <- c(cl00)
      cl0 <- data.frame("classification"=cl00)
      cl01 <- with(cl0, tapply(classification, classification, FUN = length))
      cl01 <- data.frame(classification = names(cl01), n = unname(cl01))
      cl_n <- data.frame(classification=c("higher","trivial","lower"),n=c(0,0,0))
      cll <- rbind(cl01,cl_n)
      cl1 <- with(cll, tapply(cll[,2], cll[,1], sum))
      cl1 <- data.frame(classification = names(cl1), n = unname(cl1))
      cl1$proportion <- round((cl1$n)/sum(cl1$n),4)
      cl1$d <- abs((1/3) - cl1$proportion)
      di <- round((1-((0.5*sum(cl1$d))/(2/3)))*100,1)
      dii <- data.frame("Removed study"=i,"DI"=di)
      di_sens <- rbind(di_sens,dii)
    }
    
    ds_sens <- data.frame("study"=c(),"ds"=c())
    
    for(i in 1:N) {
      sim_studd <- as.data.frame(x[,-i])
      NN <- ncol(sim_studd)
      
      sim_studa <- reshape(data=sim_studd,direction="long",
                           varying=colnames(sim_studd)[1:NN],v.names="effect")
      sim_studa <- sim_studa[,1:2]
      sim_studa1 <- table(sim_studa)
      nc <- dim(sim_studa1)[2]
      if (nc == 1) {
        ds <- 0
      }
      else {
        if (nc == 2 & N>2) {
          suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
          v <- ci_cramersv(X2)
          ds <- (sqrt((as.numeric(v$estimate)*as.numeric(v$estimate))/2))*100
        }
        else {
          suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
          v <- ci_cramersv(X2)
          ds <- round((as.numeric(v$estimate))*100,1)
        }
      }
      #
      dsi <- data.frame("Removed study"=i,"ASI"=ds)
      ds_sens <- rbind(ds_sens,dsi)
    }
  }
  else {
    di_stud <- data.frame("study"=c(),"di"=c())
    
    for(i in 1:N) {
      cl00 <- c(simsm[,i])
      cl0 <- data.frame("classification"=cl00)
      cl01 <- with(cl0, tapply(classification, classification, FUN = length))
      cl01 <- data.frame(classification = names(cl01), n = unname(cl01))
      cl_n <- data.frame(classification=
                           c("large (higher)",
                             "moderate (higher)",
                             "small (higher)",
                             "not meaningful",
                             "small (lower)",
                             "moderate (lower)",
                             "large (lower)"),n=c(0,0,0,0,0,0,0))
      cll <- rbind(cl01,cl_n)
      cl1 <- with(cll, tapply(cll[,2], cll[,1], sum))
      cl1 <- data.frame(classification = names(cl1), n = unname(cl1))
      cl1$proportion <- round((cl1$n)/sum(cl1$n),4)
      cl1$d <- abs((1/7) - cl1$proportion)
      di <- round((1-((0.5*sum(cl1$d))/(6/7)))*100,1)
      dii <- data.frame("Study"=i,"DI"=di)
      di_stud <- rbind(di_stud,dii)
    }
    
    di_sens <- data.frame("study"=c(),"di"=c())
    
    for(i in 1:N) {
      cl00 <- as.matrix(simsm[,-i])
      cl00 <- c(cl00)
      cl0 <- data.frame("classification"=cl00)
      cl01 <- with(cl0, tapply(classification, classification, FUN = length))
      cl01 <- data.frame(classification = names(cl01), n = unname(cl01))
      cl_n <- data.frame(classification=c("large (higher)","moderate (higher)","small (higher)","not meaningful","small (lower)","moderate (lower)","large (lower)"),n=c(0,0,0,0,0,0,0))
      cll <- rbind(cl01,cl_n)
      cl1 <- with(cll, tapply(cll[,2], cll[,1], sum))
      cl1 <- data.frame(classification = names(cl1), n = unname(cl1))
      cl1$proportion <- round((cl1$n)/sum(cl1$n),4)
      cl1$d <- abs((1/7) - cl1$proportion)
      di <- round((1-((0.5*sum(cl1$d))/(6/7)))*100,1)
      dii <- data.frame("Removed study"=i,"DI"=di)
      di_sens <- rbind(di_sens,dii)
    }
    
    ds_sens <- data.frame("study"=c(),"ds"=c())
    
    for (i in 1:N) {
      sim_studd <- as.data.frame(x[,-i])
      NN <- ncol(sim_studd)
      
      sim_studa <- reshape(data=sim_studd,direction="long",varying=colnames(sim_studd)[1:NN],v.names="effect")
      sim_studa <- sim_studa[,1:2]
      sim_studa1 <- table(sim_studa)
      nc <- dim(sim_studa1)[2]
      if (nc == 1) {
        ds <- 0
      }
      else {
        if (nc < 7 & N >= 7) {
          suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
          v <- ci_cramersv(X2)
          ds <- (sqrt(((as.numeric(v$estimate))*(as.numeric(v$estimate)))/(1/((nc-1)/6))))*100
        }
        else{
          if (nc<7 & N<7 & nc<N) {
            suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
            v <- ci_cramersv(X2)
            ds <- (sqrt(((as.numeric(v$estimate))*(as.numeric(v$estimate)))/(1/((nc-1)/(N-1)))))*100
          }
          else{
            suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
            v <- ci_cramersv(X2)
            ds <- round((as.numeric(v$estimate))*100,1)
          }
        }
      }
      #
      dsi <- data.frame("Removed study"=i,"ASI"=ds)
      ds_sens <- rbind(ds_sens,dsi)
    }
  }
  #
  di_asi_sens <- data.frame(di_sens,"ASI"=ds_sens$ASI)
  output <- list(di_stud = di_stud, di_asi_sens = di_asi_sens, call = match.call())
  class(output) <- "inc_inf"
  output
}





#' @rdname inf
#' @keywords print
#' @method print inf
#' @export 

print.inf <- function (x,...) {
  cat("Decision Inconsistency index (DI) for each individual primary study:", "\n")
  cat("\n")
  print(x$di_stud)
  cat("\n")
  cat("\n")
  cat("Leave-one-out sensitivity analysis for the Decision Inconsistency index (DI) and the Across-Studies Inconsistency index (ASI): ", "\n")
  cat("\n")
  print(x$di_asi_sens)
}





#' @rdname inf
#' @method summary inf
#' @export 

summary.inf <- function(object, ...) {
  output <- object
  class(output) <- "summary.inf"
  output
}





#' @rdname inc
#' @keywords print
#' @method print summary.inf
#' @export 


print.summary.inf <- function (x,...) {
  cat("Decision Inconsistency index (DI) for each individual primary study:", "\n")
  cat("\n")
  print(x$di_stud)
  cat("\n")
  cat("\n")
  cat("Leave-one-out sensitivity analysis for the Decision Inconsistency index (DI) and the Across-Studies Inconsistency index (ASI): ", "\n")
  cat("\n")
  print(x$di_asi_sens)
}
