inc_br_sens <- function(sims,br1,br2,t,sm=c("or","rr","hr"),den=1000,y_min=0,y_max=100){
  seq_b <- seq(br1,br2,0.01)
  N <- length(seq_b)
  dis_idx <- c()
  asi_idx <- c()
  classif <- data.frame("higher"=c(),"lower"=c(),"trivial"=c())

  if(sm=="or"){
    for (i in 1:N) {
      bp <- -den*(seq_b[i]-((seq_b[i]*sims)/(1-seq_b[i]+(seq_b[i]*sims))))
      t <- abs(t)
      sim_stud <- ifelse(bp>t,"higher",ifelse(bp< -t,"lower","trivial"))
      sim_studd <- as.data.frame(sim_stud)
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
      asi_idx <- c(asi_idx,ds)


      simsm <- as.matrix(sim_stud)
      cl00 <- c(simsm)
      cl0 <- data.frame("classification"=cl00)
      cl01 <- with(cl0, tapply(classification, classification, FUN = length))
      cl01 <- data.frame(classification = names(cl01), n = unname(cl01))
      cl_n <- data.frame(classification=c("higher","trivial","lower"),n=c(0,0,0))
      cll <- rbind(cl01,cl_n)
      cl1 <- with(cll, tapply(cll[,2], cll[,1], sum))
      cl1 <- data.frame(classification = names(cl1), n = unname(cl1))
      cl1$p <- (cl1$n)/sum(cl1$n)
      cl1$d <- abs((1/3) - cl1$p)
      cl2 <- cl1[, c("classification", "p")]
      cl3 <- t(as.data.frame(cl2[,2]))
      id <- (1-((0.5*sum(cl1$d))/(2/3)))*100
      dis_idx <- c(dis_idx,id)
      classif <- rbind(classif,as.data.frame(cl3))
    }
  }else{
    if(sm=="rr"){
      for (i in 1:N) {
        bp <- (seq_b[i]*sims - seq_b[i])*den
        t <- abs(t)
        sim_stud <- ifelse(bp>t,"higher",ifelse(bp< -t,"lower","trivial"))
        sim_studd <- as.data.frame(sim_stud)
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
        asi_idx <- c(asi_idx,ds)


        simsm <- as.matrix(sim_stud)
        cl00 <- c(simsm)
        cl0 <- data.frame("classification"=cl00)
        cl01 <- with(cl0, tapply(classification, classification, FUN = length))
        cl01 <- data.frame(classification = names(cl01), n = unname(cl01))
        cl_n <- data.frame(classification=c("higher","trivial","lower"),n=c(0,0,0))
        cll <- rbind(cl01,cl_n)
        cl1 <- with(cll, tapply(cll[,2], cll[,1], sum))
        cl1 <- data.frame(classification = names(cl1), n = unname(cl1))
        cl1$p <- (cl1$n)/sum(cl1$n)
        cl1$d <- abs((1/3) - cl1$p)
        cl2 <- cl1[, c("classification", "p")]
        cl3 <- t(as.data.frame(cl2[,2]))
        id <- (1-((0.5*sum(cl1$d))/(2/3)))*100
        dis_idx <- c(dis_idx,id)
        classif <- rbind(classif,as.data.frame(cl3))
      }}else{
        if(sm=="hr"){
          for(i in 1:N){
            bp <- -(seq_b[i]-(1-(1-seq_b[i])^sims))*den
            t <- abs(t)
            sim_stud <- ifelse(bp>t,"higher",ifelse(bp< -t,"lower","trivial"))
            sim_studd <- as.data.frame(sim_stud)
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
            asi_idx <- c(asi_idx,ds)


            simsm <- as.matrix(sim_stud)
            cl00 <- c(simsm)
            cl0 <- data.frame("classification"=cl00)
            cl01 <- with(cl0, tapply(classification, classification, FUN = length))
            cl01 <- data.frame(classification = names(cl01), n = unname(cl01))
            cl_n <- data.frame(classification=c("higher","trivial","lower"),n=c(0,0,0))
            cll <- rbind(cl01,cl_n)
            cl1 <- with(cll, tapply(cll[,2], cll[,1], sum))
            cl1 <- data.frame(classification = names(cl1), n = unname(cl1))
            cl1$p <- (cl1$n)/sum(cl1$n)
            cl1$d <- abs((1/3) - cl1$p)
            cl2 <- cl1[, c("classification", "p")]
            cl3 <- t(as.data.frame(cl2[,2]))
            id <- (1-((0.5*sum(cl1$d))/(2/3)))*100
            dis_idx <- c(dis_idx,id)
            classif <- rbind(classif,as.data.frame(cl3))
          }}else{
            stop("Please provide a valid effect size measure (sm): sm should be stated as 'or' (odds ratio), 'rr' (risk ratio) or 'hr' (hazard ratio)")
          }
      }
  }
  plot(x=seq_b,y=dis_idx, type="l",ylim=c(y_min,y_max),xlab="Baseline risk",ylab="Decision Inconsistency index (%)")
  title('Decision Inconsistency index according to the baseline risk')

  plot(x=seq_b,y=asi_idx, type="l",ylim=c(y_min,y_max),xlab="Baseline risk",ylab="Across-Studies Inconsistency index (%)")
  title('Across-Studies Inconsistency index according to the baseline risk')

  plot(x=seq_b,y=1-classif[,1],type="l",ylim=c(0,1),xlab="Baseline risk",ylab="Proportion of simulations")
  polygon(c(seq_b, rev(seq_b)), c(seq_b*0, rev(classif[,2])),col = "#6BD7AF",lty=0)
  polygon(c(seq_b, rev(seq_b)), c(classif[,2], rev(1-classif[,1])),col = "gray85",lty=0)
  polygon(c(seq_b, rev(seq_b)), c(1-classif[,1], seq_b/seq_b),col = "red",lty=0)
  lines(x=seq_b,y=1-classif[,1],ylim=c(0,1))
  lines(x=seq_b,y=classif[,2],ylim=c(0,1))
  mtext("red: higher; grey: trivial; green: lower",side = 3, adj = 1)
  title('Proportion of simulations higher, lower and within the decision threshold values')


  inc_sens <- data.frame("BaselineRisk"=seq_b,"DI"=round(dis_idx,1),"ASI"=round(asi_idx,1))
  inc_sens
}
