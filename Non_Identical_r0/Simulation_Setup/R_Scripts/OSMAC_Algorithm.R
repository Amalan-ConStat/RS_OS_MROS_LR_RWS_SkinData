# Cordeiro Bias Estimation ----
Cordeiro<-function(XData,With_bias)
{
  p <- as.vector(invlogit(XData%*%as.vector(With_bias)))
  W <- diag(p*(1-p))
  inverse_term <- solve(t(XData)%*%W%*%XData)
  
  Term1 <- inverse_term%*%t(XData)%*%W
  Term2 <- diag(diag(XData%*%(inverse_term)%*%t(XData))) %*% as.matrix(p-0.5)
  
  bias <- as.vector(Term1%*%Term2)
  return(bias)
} 

# Newton Rhapson Algorithm for MLE ----
getMLE <- function(x, y, w) {
  beta <- rep(0, ncol(x))
  loop  <- 1
  Loop  <- 200
  msg <- "NA"
  while (loop <= Loop) {
    pr <- c(1 - 1 / (1 + exp(x %*% beta)))
    H <- t(x) %*% (pr * (1 - pr) * w * x)
    S <- colSums((y - pr) * w * x)
    tryCatch(
      {shs <- NA
      shs <- solve(H, S) }, 
      error=function(e){
        cat("\n ERROR :", loop, conditionMessage(e), "\n")})
    if (is.na(shs[1])) {
      msg <- "Not converge"
      beta <- loop <- NA
      break
    }
    beta.new <- beta + shs
    tlr  <- sum((beta.new - beta)^2)
    beta  <- beta.new
    if(tlr < 0.000001) {
      msg <- "Successful convergence"
      break
    }
    if (loop == Loop)
      warning("Maximum iteration reached")
    loop  <- loop + 1
  }
  list(par=beta, message=msg, iter=loop)
}

# Two step OSMAC ----
AlgTwoStp <- function(r1=r1, r2=r2,Y,X,Real_Data,n,Model) {
  if(Model=="Real_Model")
  {
    n1 <- sum(Y)
    n0 <- n - n1
    PI.prop <- rep(1/(2*n0), n)
    PI.prop[Y==1] <- 1/(2*n1)
    idx.prop <- sample(1:n, r1, T, PI.prop)
    x.prop <- X[idx.prop,]
    y.prop <- Y[idx.prop]
    pinv.prop <- n
    pinv.prop <- 1/PI.prop[idx.prop]
    fit.prop <- getMLE(x=x.prop, y=y.prop, w=pinv.prop)
    beta.prop <- fit.prop$par
    if (is.na(beta.prop[1]))
      return(list(opt=NA, msg="first stage not converge"))
    P.prop  <- 1 - 1 / (1 + exp(X %*% beta.prop))
    
    beta.mVc<-matrix(nrow = length(r2),ncol = ncol(X) )
    Utility_mVc<-matrix(nrow = length(r2),ncol = 3 )
    Bias_mVc<-matrix(nrow = length(r2),ncol = ncol(X) )
    
    beta.mMSE<-matrix(nrow = length(r2),ncol = ncol(X) )
    Utility_mMSE<-matrix(nrow = length(r2),ncol = 3 )
    Bias_mMSE<-matrix(nrow = length(r2),ncol = ncol(X) )
    
    #Sample.mMSE<-list()
    #Sample.mVc<-list()
    
    ## mVc
    PI.mVc <- sqrt((Y - P.prop)^2 * rowSums(X^2))
    PI.mVc <- PI.mVc / sum(PI.mVc)
    
    ## mMSE
    p.prop <- P.prop[idx.prop]
    w.prop <- p.prop * (1 - p.prop)
    W.prop <- solve(t(x.prop) %*% (x.prop * w.prop * pinv.prop))
    PI.mMSE <- sqrt((Y - P.prop)^2 * rowSums((X%*%W.prop)^2))
    PI.mMSE <- PI.mMSE / sum(PI.mMSE)
    
    for (i in 1:length(r2)) 
    {
      ## mVc
      idx.mVc <- sample(1:n, r2[i]-r1, T, PI.mVc)
      x.mVc <- X[c(idx.mVc, idx.prop),]
      y.mVc <- Y[c(idx.mVc, idx.prop)]
      pinv.mVc <- c(1 / PI.mVc[idx.mVc], pinv.prop)
      fit.mVc <- getMLE(x=x.mVc, y=y.mVc, w=pinv.mVc)
      
      #Sample.mVc[[i]]<-cbind(r2[i],y.mVc,1/pinv.mVc,x.mVc)
      
      ru <- length(y.mVc)
      beta.mVc[i,] <- fit.mVc$par
      
      pi<- invlogit(x.mVc %*% beta.mVc[i,])
      W<-diag(as.vector(pi*(1-pi)*pinv.mVc))
      Mx<-(t(x.mVc) %*% W %*% x.mVc)
      Mx<-solve(Mx)
      Middle<-diag(((as.vector(y.mVc)-as.vector(pi))*as.vector(pinv.mVc))^2)
      V_Temp<-(t(x.mVc)%*%Middle%*%x.mVc)
      
      V_Final<-Mx %*% V_Temp %*% Mx
      
      Utility_mVc[i,]<-cbind(r2[i],tr(V_Final),det(solve(V_Final)))
      Bias_mVc[i,]<-Cordeiro(XData=x.mVc,With_bias = beta.mVc[i,])
      
      ## mMSE
      idx.mMSE <- sample(1:n, r2[i]-r1, T, PI.mMSE)
      x.mMSE <- X[c(idx.mMSE, idx.prop),]
      y.mMSE <- Y[c(idx.mMSE, idx.prop)]
      pinv.mMSE <- c(1 / PI.mMSE[idx.mMSE], pinv.prop)
      fit.mMSE <- getMLE(x=x.mMSE, y=y.mMSE, w=pinv.mMSE)
      
      #Sample.mMSE[[i]]<-cbind(r2[i],y.mMSE,1/pinv.mMSE,x.mMSE)
      
      ru <- length(y.mMSE)
      beta.mMSE[i,] <- fit.mMSE$par
      
      pi<-invlogit(x.mMSE %*% beta.mMSE[i,])
      W<-diag(as.vector(pi*(1-pi)*pinv.mMSE))
      Mx<-(t(x.mMSE) %*% W %*% x.mMSE)
      Mx<-solve(Mx)
      Middle<-diag(((as.vector(y.mMSE)-as.vector(pi))*as.vector(pinv.mMSE))^2)
      V_Temp<-(t(x.mMSE)%*%Middle%*%x.mMSE) 
      
      V_Final<-Mx %*% V_Temp %*% Mx
      
      Utility_mMSE[i,]<-cbind(r2[i],tr(V_Final),det(solve(V_Final)))
      Bias_mMSE[i,]<-Cordeiro(XData=x.mMSE,With_bias = beta.mMSE[i,])
    }
    
    #Sample_mVc<-do.call(rbind,Sample.mVc)
    #Sample_mMSE<-do.call(rbind,Sample.mMSE)
    
    opt <- list("Est_Theta_mMSE"=cbind(r2,beta.mMSE),"Utility_mMSE"=Utility_mMSE,"Bias_mMSE"=cbind(r2,Bias_mMSE),
                #"Sample_mMSE"=Sample_mMSE,
                "Est_Theta_mVc"=cbind(r2,beta.mVc),"Utility_mVc"=Utility_mVc,"Bias_mVc"=cbind(r2,Bias_mVc)#,
                #"Sample_mVc"=Sample_mVc
                )
    msg <- c(fit.prop$message, fit.mMSE$message, fit.mVc$message)
    return(list(opt=opt, msg=msg)) 
  }
  else{
    x_Real<-Real_Data[,-1]
    
    n1 <- sum(Y)
    n0 <- n - n1
    PI.prop <- rep(1/(2*n0), n)
    PI.prop[Y==1] <- 1/(2*n1)
    idx.prop <- sample(1:n, r1, T, PI.prop)
    x.prop <- X[idx.prop,]
    y.prop <- Y[idx.prop]
    x_Real.prop <- x_Real[idx.prop,] #
    pinv.prop <- n
    pinv.prop <- 1/PI.prop[idx.prop]
    fit.prop <- getMLE(x=x.prop, y=y.prop, w=pinv.prop)
    fit_Real.prop <- getMLE(x=x_Real.prop, y=y.prop, w=pinv.prop) #
    beta.prop <- fit.prop$par
    beta_Real.prop <- fit_Real.prop$par #
    if (is.na(beta.prop[1]) & is.na(beta_Real.prop[1]))
      return(list(opt=NA, msg="first stage not converge"))
    P.prop  <- 1 - 1 / (1 + exp(X %*% beta.prop))
    P_Real.prop  <- 1 - 1 / (1 + exp(x_Real %*% beta_Real.prop)) #
    
    beta.mVc<-matrix(nrow = length(r2),ncol = ncol(x_Real) )
    Utility_mVc<-matrix(nrow = length(r2),ncol = 3 )
    Bias_mVc<-matrix(nrow = length(r2),ncol = ncol(x_Real) )
    
    beta.mMSE<-matrix(nrow = length(r2),ncol = ncol(x_Real) )
    Utility_mMSE<-matrix(nrow = length(r2),ncol = 3 )
    Bias_mMSE<-matrix(nrow = length(r2),ncol = ncol(x_Real) )
    
    #Sample.mMSE<-list()
    #Sample.mVc<-list()
    
    ## mVc
    PI.mVc <- sqrt((Y - P.prop)^2 * rowSums(X^2))
    PI.mVc <- PI.mVc / sum(PI.mVc)
    
    ## mMSE
    p.prop <- P.prop[idx.prop]
    w.prop <- p.prop * (1 - p.prop)
    W.prop <- solve(t(x.prop) %*% (x.prop * w.prop * pinv.prop))
    PI.mMSE <- sqrt((Y - P.prop)^2 * rowSums((X%*%W.prop)^2))
    PI.mMSE <- PI.mMSE / sum(PI.mMSE)
    
    for (i in 1:length(r2)) 
    {
      ## mVc
      idx.mVc <- sample(1:n, r2[i]-r1, T, PI.mVc)
      
      x_Real.mVc <- x_Real[c(idx.mVc, idx.prop),]
      y_Real.mVc <- Y[c(idx.mVc, idx.prop)]
      pinv_Real.mVc <- c(1 / PI.mVc[idx.mVc], pinv.prop)
      fit_Real.mVc <- getMLE(x=x_Real.mVc, y=y_Real.mVc, w=pinv_Real.mVc)
      
      #Sample.mVc[[i]]<-cbind(r2[i],y_Real.mVc,1/pinv_Real.mVc,x_Real.mVc)
      
      ru <- length(y_Real.mVc)
      beta.mVc[i,] <- fit_Real.mVc$par
      
      pi<- invlogit(x_Real.mVc %*% beta.mVc[i,])
      W<-diag(as.vector(pi*(1-pi)*pinv_Real.mVc))
      Mx<-(t(x_Real.mVc) %*% W %*% x_Real.mVc)
      Mx<-solve(Mx)
      Middle<-diag(((as.vector(y_Real.mVc)-as.vector(pi))*as.vector(pinv_Real.mVc))^2)
      V_Temp<-(t(x_Real.mVc)%*%Middle%*%x_Real.mVc)
      
      V_Final<-Mx %*% V_Temp %*% Mx
      
      Utility_mVc[i,]<-cbind(r2[i],tr(V_Final),det(solve(V_Final)))
      Bias_mVc[i,]<-Cordeiro(XData=x_Real.mVc,With_bias = beta.mVc[i,])
      
      ## mMSE
      idx.mMSE <- sample(1:n, r2[i]-r1, T, PI.mMSE)
      x_Real.mMSE <- x_Real[c(idx.mMSE, idx.prop),]
      y_Real.mMSE <- Y[c(idx.mMSE, idx.prop)]
      pinv_Real.mMSE <- c(1 / PI.mMSE[idx.mMSE], pinv.prop)
      fit_Real.mMSE <- getMLE(x=x_Real.mMSE, y=y_Real.mMSE, w=pinv_Real.mMSE)
      
      #Sample.mMSE[[i]]<-cbind(r2[i],y_Real.mMSE,1/pinv_Real.mMSE,x_Real.mMSE)
      
      ru <- length(y_Real.mMSE)
      beta.mMSE[i,] <- fit_Real.mMSE$par
      
      pi<-invlogit(x_Real.mMSE %*% beta.mMSE[i,])
      W<-diag(as.vector(pi*(1-pi)*pinv_Real.mMSE))
      Mx<-(t(x_Real.mMSE) %*% W %*% x_Real.mMSE)
      Mx<-solve(Mx)
      Middle<-diag(((as.vector(y_Real.mMSE)-as.vector(pi))*as.vector(pinv_Real.mMSE))^2)
      V_Temp<-(t(x_Real.mMSE)%*%Middle%*%x_Real.mMSE)
      
      V_Final<-Mx %*% V_Temp %*% Mx
      
      Utility_mMSE[i,]<-cbind(r2[i],tr(V_Final),det(solve(V_Final)))
      Bias_mMSE[i,]<-Cordeiro(XData=x_Real.mMSE,With_bias = beta.mMSE[i,])
    }
    
    if(anyNA(beta.mMSE) || anyNA(beta.mVc) || anyNA(Utility_mMSE) || anyNA(Utility_mVc) || anyNA(Bias_mMSE) || anyNA(Bias_mVc) )
    {
      stop("There are NA or NaN values")
    }
    
    #Sample_mVc<-do.call(rbind,Sample.mVc)
    #Sample_mMSE<-do.call(rbind,Sample.mMSE)
    
    opt <- list("Est_Theta_mMSE"=cbind(r2,beta.mMSE),"Utility_mMSE"=Utility_mMSE,"Bias_mMSE"=cbind(r2,Bias_mMSE),
                #"Sample_mMSE"=Sample_mMSE,
                "Est_Theta_mVc"=cbind(r2,beta.mVc),"Utility_mVc"=Utility_mVc,"Bias_mVc"=cbind(r2,Bias_mVc)#,
                #"Sample_mVc"=Sample_mVc
                )
    msg <- c(fit_Real.prop$message, fit_Real.mMSE$message, fit_Real.mVc$message)
    return(list(opt=opt, msg=msg)) 
  }
}
