#' @export inc_subset

inc_subset <- function(object,data,variable,subset,exc=FALSE){
  
  catch <- function(argname, matchcall, data, encl)
    eval(matchcall[[match(argname, names(matchcall))]], data, enclos = encl)
  
  nulldata <- is.null(data)
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  
  if (nulldata)
    data <- sfsp

  variable <- catch("variable", mc, data, sfsp)


x <- object$stud_class
simsm <- as.matrix(x)

if(exc){
  subgroup_inc <- which(variable!=subset)
}else{
  subgroup_inc <- which(variable==subset)
}


simsm_i <- simsm[,c(subgroup_inc)]


if(is.null(object$call$t1)){
  
  cl00 <- c(simsm_i)
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
  di <- (1-((0.5*sum(cl1$d))/(2/3)))*100
  
  
  sim_studd <- as.data.frame(simsm_i)
  N <- ncol(sim_studd)
  
  sim_studa <- reshape(data=sim_studd,direction="long",varying=colnames(sim_studd)[1:N],v.names="effect")
  sim_studa <- sim_studa[,1:2]
  sim_studa1 <- table(sim_studa)
  nc <- dim(sim_studa1)[2]
  if(nc==1){
    ds <- 0
  }else{
    if(nc==2 & N>2){
      suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
      v <- ci_cramersv(X2)
      ds <- (sqrt((as.numeric(v$estimate)*as.numeric(v$estimate))/2))*100
    }else{
      suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
      v <- ci_cramersv(X2)
      ds <- (as.numeric(v$estimate))*100
    }}
  
  
  output <- list(ASI = round(ds,3), DI = round(di,3), class_distribution = cl2, n_stud=N, call = match.call())
  class(output) <- "inc_subset"
  output
  
}else{
  
  sim_studd <- as.data.frame(simsm_i)
  N <- ncol(sim_studd)
  
  sim_studa <- reshape(data=sim_studd,direction="long",varying=colnames(sim_studd)[1:N],v.names="effect")
  sim_studa <- sim_studa[,1:2]
  sim_studa1 <- table(sim_studa)
  nc <- dim(sim_studa1)[2]
  if(nc==1){
    ds <- 0
  }else{
    if(nc<7 & N>=7){
      suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
      v <- ci_cramersv(X2)
      ds <- (sqrt(((as.numeric(v$estimate))*(as.numeric(v$estimate)))/(1/((nc-1)/6))))*100
    }else{
      if(nc<7 & N<7 & nc<N){
        suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
        v <- ci_cramersv(X2)
        ds <- (sqrt(((as.numeric(v$estimate))*(as.numeric(v$estimate)))/(1/((nc-1)/(N-1)))))*100
      }else{
        suppressWarnings(X2 <- stats::chisq.test(sim_studa1))
        v <- ci_cramersv(X2)
        ds <- (as.numeric(v$estimate))*100
      }}}
  
  
  simsm <- as.matrix(simsm_i)
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
  di <- (1-((0.5*sum(cl1$d))/(6/7)))*100
  
  output <- list(ASI = round(ds,3), DI = round(di,3), class_distribution = cl2, n_stud=N, call = match.call())
  class(output) <- "inc_subset"
  output
}}


print.inc_subset <- function (x,...){
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
  cat("Number of primary studies in this subset: ", x$n_stud)
  cat("\n")
}

summary.inc_subset <- function (object,...){
  x <- object
  cat("Decision Inconsistency index (DI): ", round(x$DI,1), "%")
  cat("\n")
  cat("Across-Studies Inconsistency (ASI): ", round(x$ASI,1), "%")
  cat("\n")
  cat("\n")
  cat("Proportion of simulations by classification in relation to the decision threshold:", "\n")
  print(x$class_distribution)
  cat("\n")
  cat("\n")
  cat("Number of primary studies in this subset: ", x$n_stud)
  cat("\n")
}
