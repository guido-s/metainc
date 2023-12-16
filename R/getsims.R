getsims <- function(obj,param,package=c("R2OpenBUGS","R2jags","rjags","brms","metafor","meta"),nsims=10000,t_exp=FALSE){
  if(missing(package)){
    if(substr(obj$call[1],1,3)=="rma"){
      N <- length(obj$vi.f)
      m0 <- blup(obj)
      yi <- m0$pred
      sei <- m0$se

      db <- data.frame("V1"=seq(1:nsims))
      for(i in 1:N){
        a <- rnorm(nsims,yi[i],sei[i])
        db <- cbind(db,a)
      }

      sims <- db[,2:(N+1)]
      sims <- as.matrix(sims)

      if(t_exp){return(exp(sims))} else {return(sims)}
    }else{
        if(class(obj)[2]=="meta"){
          obj1 <- metareg(obj,1)
          N <- length(obj1$vi.f)
          m0 <- blup(obj1)
          yi <- m0$pred
          sei <- m0$se

          db <- data.frame("V1"=seq(1:nsims))
          for(i in 1:N){
            a <- rnorm(nsims,yi[i],sei[i])
            db <- cbind(db,a)
          }

          sims <- db[,2:(N+1)]
          sims <- as.matrix(sims)

          if(t_exp){return(exp(sims))} else {return(sims)}
      }else{
        if(class(obj)[1]=="brmsfit"){
          sims_df <- as.data.frame(obj$fit)
          sims_l <- sims_df[, grepl(paste0("r_",param), names(sims_df))]
          sims <- as.matrix(sims_l)
          if(t_exp){return(exp(sims))} else {return(sims)}
        }else{
    stop("Please indicate the name of the package for meta-analysis")
  }}}}else{
  if(package=="R2OpenBUGS"){
    sims <- as.data.frame(obj$sims.list[param])
    sims <- as.matrix(sims)
    if(t_exp){return(exp(sims))} else {return(sims)}
  } else{
    if(package=="R2jags"){
      sims <- as.data.frame(obj$BUGSoutput$sims.list[param])
      sims <- as.matrix(sims)
      if(t_exp){return(exp(sims))} else {return(sims)}
    }else{
      if(package=="rjags"){
        N <- as.numeric(nrow(obj[[param]]))
        Np <- N+1
        sims2 <- data.frame(study=c(),sims0=c())
        for(i in 1:N){
          sims0 <- c(obj[[param]][i,,])
          NN <- length(sims0)
          sims1 <- data.frame("n_sim"=c(1:NN),"study"=i,"sims0"=sims0)
          sims2 <- rbind(sims2,sims1)
        }
        sims3 <- reshape(data=sims2,direction="wide",idvar="n_sim",timevar="study")
        sims3 <- sims3[,2:Np]
        sims <- as.matrix(sims3)
        if(t_exp){return(exp(sims))} else {return(sims)}
      }else{
        if(package=="brms"){
          sims_df <- as.data.frame(obj$fit)
          sims_l <- sims_df[, grepl(paste0("r_",param), names(sims_df))]
          sims <- as.matrix(sims_l)
          if(t_exp){return(exp(sims))} else {return(sims)}
        }else{
          if(package=="metafor"){
          N <- length(obj$vi.f)
          m0 <- blup(obj)
          yi <- m0$pred
          sei <- m0$se

          db <- data.frame("V1"=seq(1:nsims))
          for(i in 1:N){
            a <- rnorm(nsims,yi[i],sei[i])
            db <- cbind(db,a)
          }

          sims <- db[,2:(N+1)]
          sims <- as.matrix(sims)

          if(t_exp){return(exp(sims))} else {return(sims)}
        }else{
          if(package=="meta"){
            obj1 <- metareg(obj,1)
            N <- length(obj1$vi.f)
            m0 <- blup(obj1)
            yi <- m0$pred
            sei <- m0$se

            db <- data.frame("V1"=seq(1:nsims))
            for(i in 1:N){
              a <- rnorm(nsims,yi[i],sei[i])
              db <- cbind(db,a)
            }

            sims <- db[,2:(N+1)]
            sims <- as.matrix(sims)

            if(t_exp){return(exp(sims))} else {return(sims)}
          }else{
            stop("Please note that the package for meta-analysis should either be 'R2OpenBUGS', 'R2jags', 'rjags', 'brms', 'metafor' and 'meta'")
          }}}
        }
      }
    }}}
