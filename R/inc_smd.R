#' Title line
#'
#' @description
#' 
#' Brief description (more details could be provided under Details).
#'
#' @param sims R object ...
#' @param t1 ...
#' @param t2 ...
#' @param t3 ...
#' @param x R object of class \code{inc_smd} or  \code{summary.inc_smd}.
#' @param object R object of class \code{inc_smd}.
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
#' @export inc_smd


inc_smd <- function(sims, t1 = 0.2, t2 = 0.5, t3 = 0.8) {
  sim_stud <- ifelse(sims>=t3,"large (higher)",ifelse(sims>=t2,"moderate (higher)",ifelse(sims>=t1,"small (higher)",ifelse(sims> -t1,"not meaningful",ifelse(sims> -t2,"small (lower)",ifelse(sims> -t3,"moderate (lower)","large (lower)"))))))
  sim_studd <- as.data.frame(sim_stud)
  N <- ncol(sim_studd)

  sim_studa <- reshape(data=sim_studd,direction="long",varying=colnames(sim_studd)[1:N],v.names="effect")
  sim_studa <- sim_studa[,1:2]
  sim_studa1 <- table(sim_studa)
  nc <- dim(sim_studa1)[2]
  if (nc==1) {
    ds <- 0
  }else{
    if (nc<7 & N>=7) {
      suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
      v <- ci_cramersv(X2)
      ds <- (sqrt(((as.numeric(v$estimate))*(as.numeric(v$estimate)))/(1/((nc-1)/6))))*100
    }else{
      if (nc<7 & N<7 & nc<N) {
        suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
        v <- ci_cramersv(X2)
        ds <- (sqrt(((as.numeric(v$estimate))*(as.numeric(v$estimate)))/(1/((nc-1)/(N-1)))))*100
      }else{
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
  propnull <- round(length(sims[sims>0])/sum(cl1$n),3)
  di <- (1-((0.5*sum(cl1$d))/(6/7)))*100

  output <- list(ASI = round(ds,3), DI = round(di,3), class_distribution = cl2,prop_over_null = propnull, stud_class=sim_studd, call = match.call())
  class(output) <- "inc_smd"
  output
}





#' @rdname inc_smd
#' @keywords print
#' @method print inc_smd
#' @export 

print.inc_smd <- function (x,...) {
  x <- x
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





#' @rdname inc_smd
#' @method summary inc_smd
#' @export 

summary.inc_smd <- function(object, ...) {
  output <- object
  class(output) <- "summary.inc_smd"
  output
}





#' @rdname inc_smd
#' @keywords print
#' @method print summary.inc_smd
#' @export 

print.summary.inc_smd <- function (x,...) {
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
