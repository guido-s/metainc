inc_t_sens <- function(sims,br,t1_min,t1_max,t2_min,t2_max,sm=c("or","rr","hr"),den=1000,unit=1,null_effect=0){
  if(t1_max>=null_effect | t2_min<= null_effect | t1_min>=t1_max | t2_min>= t2_max) {
    stop("Invalid decision threshold values. When providing decision threshold values, please make sure that: (i) t1_min and t1_max are both lower than the null effect value, (ii) t2_min and t2_max are both higher than the null effect value, (iii) t1_max is higher than t1_min, and (iv) t2_max is higher than t2_min")
  }else{
  if(missing(br)){
    t1_seq <- seq(t1_min,t1_max,unit)
    t2_seq <- seq(t2_min,t2_max,unit)
    Nt1 <- length(t1_seq)
    Nt2 <- length(t2_seq)
    t1_seq_r <- rep(t1_seq,Nt1)
    t2_seq_r <- rep(t2_seq,Nt2)
    t2_seq_r <- t2_seq_r[order(t2_seq_r)]

    dsi=c()
    asi=c()
    prop_higher=c()
    prop_lower=c()
    prop_trivial=c()

    for(i in 1:length(t1_seq_r)){
      sim_stud <- ifelse(sims>t2_seq_r[i],"higher",ifelse(sims<t1_seq_r[i],"lower","trivial"))
      sim_studd <- as.data.frame(sim_stud)
      N <- ncol(sim_studd)

      cl <- c(sim_stud)
      cl0 <- data.frame("classification"=cl)
      cl01 <- with(cl0, tapply(classification, classification, FUN = length))
      cl01 <- data.frame(classification = names(cl01), n = unname(cl01))
      cl_n <- data.frame(classification=c("higher","trivial","lower"),n=c(0,0,0))
      cll <- rbind(cl01,cl_n)
      cl1 <- with(cll, tapply(cll[,2], cll[,1], sum))
      cl1 <- data.frame(classification = names(cl1), n = unname(cl1))
      cl1$p <- (cl1$n)/sum(cl1$n)
      cl1$d <- abs((1/3) - cl1$p)
      cl2 <- cl1[, c("classification", "p")]
      id <- (1-((0.5*sum(cl1$d))/(2/3)))*100
      list(paste0("Decision Inconsistency index: ",id),cl2)
      prop_higher <- c(prop_higher,as.numeric(cl2[1,2]))
      prop_lower <- c(prop_lower,as.numeric(cl2[2,2]))
      prop_trivial <- c(prop_trivial,as.numeric(cl2[3,2]))
      dsi <- rbind(dsi,id)


      sim_stud1=data.frame("b"=c(),"h"=c(),"t"=c())
      for(j in 1:N){
        b <- as.numeric(length(sim_studd[,j][sim_studd[,j]=="lower"]))
        h <- as.numeric(length(sim_studd[,j][sim_studd[,j]=="higher"]))
        t <- as.numeric(length(sim_studd[,j][sim_studd[,j]=="trivial"]))

        class_stud=data.frame("b"=b,"h"=h,"t"=t)
        sim_stud1 <- rbind(sim_stud1,class_stud)
      }

      sim_stud2=data.frame("e_b"=c(),"e_h"=c(),"e_t"=c())
      for(k in 1:N){
        e_b <- as.numeric((sum(sim_stud1$b) * sum(sim_stud1[k,]))/sum(sim_stud1))
        e_h <- as.numeric((sum(sim_stud1$h) * sum(sim_stud1[k,]))/sum(sim_stud1))
        e_t <- as.numeric((sum(sim_stud1$t) * sum(sim_stud1[k,]))/sum(sim_stud1))

        class_stud1=data.frame("e_b"=e_b,"e_h"=e_h,"e_t"=e_t)
        sim_stud2 <- rbind(sim_stud2,class_stud1)
      }

      sim_stud3 <- cbind(sim_stud1,sim_stud2)
      sim_stud3$dif_b <- if(sum(sim_stud3$e_b)>0){((sim_stud3$b-sim_stud3$e_b)^2)/sim_stud3$e_b}else{0}
      sim_stud3$dif_h <- if(sum(sim_stud3$e_h)>0){((sim_stud3$h-sim_stud3$e_h)^2)/sim_stud3$e_h}else{0}
      sim_stud3$dif_t <- if(sum(sim_stud3$e_t)>0){((sim_stud3$t-sim_stud3$e_t)^2)/sim_stud3$e_t}else{0}

      dsi_n <- sum(sim_stud3$dif_b) + sum(sim_stud3$dif_h) + sum(sim_stud3$dif_t)

      dsi_d <- (sum(sim_stud1)*(min(3,N)-1))

      ds <- sqrt((dsi_n/dsi_d))*100
      asi <- rbind(asi,ds)
    }

    p1 <- ggplot(mapping=aes(x=t1_seq_r,y=t2_seq_r))+geom_raster(aes(fill = dsi),interpolate=T)+scale_fill_gradient(high="red", low="green",limits=c(0,100),name="DI") + labs(x = "Lower threshold of appreciable effect", y="Upper threshold of appreciable effect",title="Decision Inconsistency index")
    p2 <- ggplot(mapping=aes(x=t1_seq_r,y=t2_seq_r))+geom_raster(aes(fill = asi),interpolate=T)+scale_fill_gradient(high="red", low="green",limits=c(0,100),name="V") + labs(x = "Lower threshold of appreciable effect", y="Upper threshold of appreciable effect",title="Across-studies Inconsistency index")
    p3 <- ggplot(mapping=aes(x=t1_seq_r,y=t2_seq_r))+geom_raster(aes(fill = prop_higher),interpolate=T)+scale_fill_gradient(low="white", high="black",limits=c(0,1),name="Proportion") + labs(x = "Lower threshold of appreciable effect", y="Upper threshold of appreciable effect",title="Proportion of simulations higher than \n the highest decision threshold value")
    p4 <- ggplot(mapping=aes(x=t1_seq_r,y=t2_seq_r))+geom_raster(aes(fill = prop_lower),interpolate=T)+scale_fill_gradient(low="white", high="black",limits=c(0,1),name="Proportion") + labs(x = "Lower threshold of appreciable effect", y="Upper threshold of appreciable effect",title="Proportion of simulations lower than \n the lowest decision threshold value")
    p5 <- ggplot(mapping=aes(x=t1_seq_r,y=t2_seq_r))+geom_raster(aes(fill = prop_trivial),interpolate=T)+scale_fill_gradient(low="white", high="black",limits=c(0,1),name="Proportion") + labs(x = "Lower threshold of appreciable effect", y="Upper threshold of appreciable effect",title="Proportion of simulations with trivial effect")

  }else{

    t1_seq <- seq(t1_min,t1_max,1)
    t2_seq <- seq(t2_min,t2_max,1)
    Nt1 <- length(t1_seq)
    Nt2 <- length(t2_seq)
    t1_seq_r <- rep(t1_seq,Nt1)
    t2_seq_r <- rep(t2_seq,Nt2)
    t2_seq_r <- t2_seq_r[order(t2_seq_r)]

    dsi=c()
    asi=c()
    prop_higher=c()
    prop_lower=c()
    prop_trivial=c()

    if(sm=="or"){
      for(i in 1:length(t1_seq_r)){
        bp <- -den*(br-((br*sims)/(1-br+(br*sims))))
        sim_stud <- ifelse(bp>t2_seq_r[i],"higher",ifelse(bp<t1_seq_r[i],"lower","trivial"))
        sim_studd <- as.data.frame(sim_stud)
        N <- ncol(sim_studd)

        cl <- c(sim_stud)
        cl0 <- data.frame("classification"=cl)
        cl01 <- with(cl0, tapply(classification, classification, FUN = length))
        cl01 <- data.frame(classification = names(cl01), n = unname(cl01))
        cl_n <- data.frame(classification=c("higher","trivial","lower"),n=c(0,0,0))
        cll <- rbind(cl01,cl_n)
        cl1 <- with(cll, tapply(cll[,2], cll[,1], sum))
        cl1 <- data.frame(classification = names(cl1), n = unname(cl1))
        cl1$p <- (cl1$n)/sum(cl1$n)
        cl1$d <- abs((1/3) - cl1$p)
        cl2 <- cl1[, c("classification", "p")]
        id <- (1-((0.5*sum(cl1$d))/(2/3)))*100
        list(paste0("Decision Inconsistency index: ",id),cl2)
        prop_higher <- c(prop_higher,as.numeric(cl2[1,2]))
        prop_lower <- c(prop_lower,as.numeric(cl2[2,2]))
        prop_trivial <- c(prop_trivial,as.numeric(cl2[3,2]))
        dsi <- rbind(dsi,id)


        sim_stud1=data.frame("b"=c(),"h"=c(),"t"=c())
        for(j in 1:N){
          b <- as.numeric(length(sim_studd[,j][sim_studd[,j]=="lower"]))
          h <- as.numeric(length(sim_studd[,j][sim_studd[,j]=="higher"]))
          t <- as.numeric(length(sim_studd[,j][sim_studd[,j]=="trivial"]))

          class_stud=data.frame("b"=b,"h"=h,"t"=t)
          sim_stud1 <- rbind(sim_stud1,class_stud)
        }

        sim_stud2=data.frame("e_b"=c(),"e_h"=c(),"e_t"=c())
        for(k in 1:N){
          e_b <- as.numeric((sum(sim_stud1$b) * sum(sim_stud1[k,]))/sum(sim_stud1))
          e_h <- as.numeric((sum(sim_stud1$h) * sum(sim_stud1[k,]))/sum(sim_stud1))
          e_t <- as.numeric((sum(sim_stud1$t) * sum(sim_stud1[k,]))/sum(sim_stud1))

          class_stud1=data.frame("e_b"=e_b,"e_h"=e_h,"e_t"=e_t)
          sim_stud2 <- rbind(sim_stud2,class_stud1)
        }

        sim_stud3 <- cbind(sim_stud1,sim_stud2)
        sim_stud3$dif_b <- ((sim_stud3$b-sim_stud3$e_b)^2)/sim_stud3$e_b
        sim_stud3$dif_h <- ((sim_stud3$h-sim_stud3$e_h)^2)/sim_stud3$e_h
        sim_stud3$dif_t <- ((sim_stud3$t-sim_stud3$e_t)^2)/sim_stud3$e_t

        dsi_n <- sum(sim_stud3$dif_b) + sum(sim_stud3$dif_h) + sum(sim_stud3$dif_t)

        dsi_d <- (sum(sim_stud1)*(min(3,N)-1))

        ds <- sqrt((dsi_n/dsi_d))*100
        asi <- rbind(asi,ds)

      }} else {
        if(sm=="rr"){
          for(i in 1:length(t1_seq_r)){
            bp <- (br*sims - br)*den
            sim_stud <- ifelse(bp>t2_seq_r[i],"higher",ifelse(bp<t1_seq_r[i],"lower","trivial"))
            sim_studd <- as.data.frame(sim_stud)
            N <- ncol(sim_studd)

            cl <- c(sim_stud)
            cl0 <- data.frame("classification"=cl)
            cl01 <- with(cl0, tapply(classification, classification, FUN = length))
            cl01 <- data.frame(classification = names(cl01), n = unname(cl01))
            cl_n <- data.frame(classification=c("higher","trivial","lower"),n=c(0,0,0))
            cll <- rbind(cl01,cl_n)
            cl1 <- with(cll, tapply(cll[,2], cll[,1], sum))
            cl1 <- data.frame(classification = names(cl1), n = unname(cl1))
            cl1$p <- (cl1$n)/sum(cl1$n)
            cl1$d <- abs((1/3) - cl1$p)
            cl2 <- cl1[, c("classification", "p")]
            id <- (1-((0.5*sum(cl1$d))/(2/3)))*100
            list(paste0("Decision Inconsistency index: ",id),cl2)
            prop_higher <- c(prop_higher,as.numeric(cl2[1,2]))
            prop_lower <- c(prop_lower,as.numeric(cl2[2,2]))
            prop_trivial <- c(prop_trivial,as.numeric(cl2[3,2]))
            dsi <- rbind(dsi,id)


            sim_stud1=data.frame("b"=c(),"h"=c(),"t"=c())
            for(j in 1:N){
              b <- as.numeric(length(sim_studd[,j][sim_studd[,j]=="lower"]))
              h <- as.numeric(length(sim_studd[,j][sim_studd[,j]=="higher"]))
              t <- as.numeric(length(sim_studd[,j][sim_studd[,j]=="trivial"]))

              class_stud=data.frame("b"=b,"h"=h,"t"=t)
              sim_stud1 <- rbind(sim_stud1,class_stud)
            }

            sim_stud2=data.frame("e_b"=c(),"e_h"=c(),"e_t"=c())
            for(k in 1:N){
              e_b <- as.numeric((sum(sim_stud1$b) * sum(sim_stud1[k,]))/sum(sim_stud1))
              e_h <- as.numeric((sum(sim_stud1$h) * sum(sim_stud1[k,]))/sum(sim_stud1))
              e_t <- as.numeric((sum(sim_stud1$t) * sum(sim_stud1[k,]))/sum(sim_stud1))

              class_stud1=data.frame("e_b"=e_b,"e_h"=e_h,"e_t"=e_t)
              sim_stud2 <- rbind(sim_stud2,class_stud1)
            }

            sim_stud3 <- cbind(sim_stud1,sim_stud2)
            sim_stud3$dif_b <- ((sim_stud3$b-sim_stud3$e_b)^2)/sim_stud3$e_b
            sim_stud3$dif_h <- ((sim_stud3$h-sim_stud3$e_h)^2)/sim_stud3$e_h
            sim_stud3$dif_t <- ((sim_stud3$t-sim_stud3$e_t)^2)/sim_stud3$e_t

            dsi_n <- sum(sim_stud3$dif_b) + sum(sim_stud3$dif_h) + sum(sim_stud3$dif_t)

            dsi_d <- (sum(sim_stud1)*(min(3,N)-1))

            ds <- sqrt((dsi_n/dsi_d))*100
            asi <- rbind(asi,ds)
          }}else{
            if(sm=="hr"){
              for(i in 1:length(t1_seq_r)){
                bp <- -(br-(1-(1-br)^sims))*den
                sim_stud <- ifelse(bp>t2_seq_r[i],"higher",ifelse(bp<t1_seq_r[i],"lower","trivial"))
                sim_studd <- as.data.frame(sim_stud)
                N <- ncol(sim_studd)

                cl <- c(sim_stud)
                cl0 <- data.frame("classification"=cl)
                cl01 <- with(cl0, tapply(classification, classification, FUN = length))
                cl01 <- data.frame(classification = names(cl01), n = unname(cl01))
                cl_n <- data.frame(classification=c("higher","trivial","lower"),n=c(0,0,0))
                cll <- rbind(cl01,cl_n)
                cl1 <- with(cll, tapply(cll[,2], cll[,1], sum))
                cl1 <- data.frame(classification = names(cl1), n = unname(cl1))
                cl1$p <- (cl1$n)/sum(cl1$n)
                cl1$d <- abs((1/3) - cl1$p)
                cl2 <- cl1[, c("classification", "p")]
                id <- (1-((0.5*sum(cl1$d))/(2/3)))*100
                list(paste0("Decision Inconsistency index: ",id),cl2)
                prop_higher <- c(prop_higher,as.numeric(cl2[1,2]))
                prop_lower <- c(prop_lower,as.numeric(cl2[2,2]))
                prop_trivial <- c(prop_trivial,as.numeric(cl2[3,2]))
                dsi <- rbind(dsi,id)


                sim_stud1=data.frame("b"=c(),"h"=c(),"t"=c())
                for(j in 1:N){
                  b <- as.numeric(length(sim_studd[,j][sim_studd[,j]=="lower"]))
                  h <- as.numeric(length(sim_studd[,j][sim_studd[,j]=="higher"]))
                  t <- as.numeric(length(sim_studd[,j][sim_studd[,j]=="trivial"]))

                  class_stud=data.frame("b"=b,"h"=h,"t"=t)
                  sim_stud1 <- rbind(sim_stud1,class_stud)
                }

                sim_stud2=data.frame("e_b"=c(),"e_h"=c(),"e_t"=c())
                for(k in 1:N){
                  e_b <- as.numeric((sum(sim_stud1$b) * sum(sim_stud1[k,]))/sum(sim_stud1))
                  e_h <- as.numeric((sum(sim_stud1$h) * sum(sim_stud1[k,]))/sum(sim_stud1))
                  e_t <- as.numeric((sum(sim_stud1$t) * sum(sim_stud1[k,]))/sum(sim_stud1))

                  class_stud1=data.frame("e_b"=e_b,"e_h"=e_h,"e_t"=e_t)
                  sim_stud2 <- rbind(sim_stud2,class_stud1)
                }

                sim_stud3 <- cbind(sim_stud1,sim_stud2)
                sim_stud3$dif_b <- ((sim_stud3$b-sim_stud3$e_b)^2)/sim_stud3$e_b
                sim_stud3$dif_h <- ((sim_stud3$h-sim_stud3$e_h)^2)/sim_stud3$e_h
                sim_stud3$dif_t <- ((sim_stud3$t-sim_stud3$e_t)^2)/sim_stud3$e_t

                dsi_n <- sum(sim_stud3$dif_b) + sum(sim_stud3$dif_h) + sum(sim_stud3$dif_t)

                dsi_d <- (sum(sim_stud1)*(min(3,N)-1))

                ds <- sqrt((dsi_n/dsi_d))*100
                asi <- rbind(asi,ds)
              }} else{
              stop("Please provide a valid effect size measure (sm): sm should be stated as 'or' (odds ratio), 'rr' (risk ratio) or 'hr' (hazard ratio)")
            }
          }
      }
    p1 <- ggplot(mapping=aes(x=t1_seq_r,y=t2_seq_r))+geom_raster(aes(fill = dsi),interpolate=T)+scale_fill_gradient(high="red", low="green",limits=c(0,100),name="DI") + labs(x = "Lower threshold of appreciable effect", y="Upper threshold of appreciable effect",title="Decision Inconsistency index")
    p2 <- ggplot(mapping=aes(x=t1_seq_r,y=t2_seq_r))+geom_raster(aes(fill = asi),interpolate=T)+scale_fill_gradient(high="red", low="green",limits=c(0,100),name="V") + labs(x = "Lower threshold of appreciable effect", y="Upper threshold of appreciable effect",title="Across-studies Inconsistency index")
    p3 <- ggplot(mapping=aes(x=t1_seq_r,y=t2_seq_r))+geom_raster(aes(fill = prop_higher),interpolate=T)+scale_fill_gradient(low="white", high="black",limits=c(0,1),name="Proportion") + labs(x = "Lower threshold of appreciable effect", y="Upper threshold of appreciable effect",title="Proportion of simulations higher than \n the highest decision threshold value")
    p4 <- ggplot(mapping=aes(x=t1_seq_r,y=t2_seq_r))+geom_raster(aes(fill = prop_lower),interpolate=T)+scale_fill_gradient(low="white", high="black",limits=c(0,1),name="Proportion") + labs(x = "Lower threshold of appreciable effect", y="Upper threshold of appreciable effect",title="Proportion of simulations lower than \n the lowest decision threshold value")
    p5 <- ggplot(mapping=aes(x=t1_seq_r,y=t2_seq_r))+geom_raster(aes(fill = prop_trivial),interpolate=T)+scale_fill_gradient(low="white", high="black",limits=c(0,1),name="Proportion") + labs(x = "Lower threshold of appreciable effect", y="Upper threshold of appreciable effect",title="Proportion of simulations with trivial effect")
  }

  inc_sens <- data.frame("LowerThreshold"=t1_seq_r,"UpperThreshold"=t2_seq_r,"DI"=round(dsi,1),"ASI"=round(asi,1))

  output <- list(inc_sens = inc_sens, p1=p1,p2=p2,p3=p3,p4=p4,p5=p5)
  class(output) <- "inc_t_sens"
  output
  }}

print.inc_t_sens <- function (x,...){
  x <- x
  plot(x$p1)
  plot(x$p2)
  plot(x$p3)
  plot(x$p4)
  plot(x$p5)
}
