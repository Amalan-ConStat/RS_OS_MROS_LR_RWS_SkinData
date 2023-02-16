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

# Two Step OSMAC ----
OSMAC_MF <- function(r1,r2,Y,X,n,alpha,combs,All_Covariates) {
  n1 <- sum(Y)
  n0 <- n - n1
  PI.prop <- rep(1/(2*n0), n)
  PI.prop[Y==1] <- 1/(2*n1)
  idx.prop <- sample(1:n, r1, T, PI.prop)
  
  y.prop <- Y[idx.prop]
  pinv.prop <- n
  pinv.prop <- 1/PI.prop[idx.prop]
  
  beta.prop<-list()
  beta.prop<-lapply(1:length(combs),function(j) {
    getMLE(x=X[idx.prop,colnames(X) %in% combs[[j]]  ], y=y.prop, w=pinv.prop)})
   
  for (j in 1:length(combs)) 
  {
    if (anyNA(beta.prop[[j]]$par))
    {
      stop("There are NA or NaN values")
    }
  }
  
  P.prop<-lapply(1:length(combs), function(j)
    1 - 1 / (1 + exp(X[,colnames(X) %in% combs[[j]]] %*% beta.prop[[j]]$par)))
  
  beta.mVc<-Utility_mVc<-Bias_mVc<-list()
  beta.mMSE<-Utility_mMSE<-Bias_mMSE<-list()
  
  for (a in 1:length(combs)) 
  {
    beta.mVc[[a]]<-matrix(nrow = length(r2),ncol = length(combs[[a]])) #
    Utility_mVc[[a]]<-matrix(nrow = length(r2),ncol = 3 ) #
    Bias_mVc[[a]]<-matrix(nrow = length(r2),ncol = length(combs[[a]])) #
    
    beta.mMSE[[a]]<-matrix(nrow = length(r2),ncol = length(combs[[a]])) #
    Utility_mMSE[[a]]<-matrix(nrow = length(r2),ncol = 3) #
    Bias_mMSE[[a]]<-matrix(nrow = length(r2),ncol = length(combs[[a]]))  #
  }
  
  Sample.mMSE<-list()
  Sample.mVc<-list()
  
  ## mVc
  PI.mVc<-lapply(1:length(combs),function(j){
    PI.mVc <- sqrt((Y - P.prop[[j]])^2 * rowSums(X[,All_Covariates %in% combs[[j]] ]^2)) #
    return(PI.mVc / sum(PI.mVc))})
  
  PI.mVc<-matrix(unlist(PI.mVc),nrow = n,byrow = FALSE)
  pjoin.mVc<-rowWeightedMeans(PI.mVc,w=alpha)
  
  ## mMSE
  p.prop <-lapply(1:length(combs),function(a) P.prop[[a]][idx.prop])  #
  w.prop <-lapply(1:length(combs),function(a) p.prop[[a]] * (1 - p.prop[[a]]))  #
  
  W.prop <-lapply(1:length(combs),function(a) {
    solve(t(X[idx.prop,All_Covariates %in% combs[[a]] ]) %*% 
            (X[idx.prop,All_Covariates %in% combs[[a]] ] * w.prop[[a]] * pinv.prop)) })   #
  
  PI.mMSE <-lapply(1:length(combs), function(a){
    PI.mMSE<-sqrt((Y - P.prop[[a]])^2 * 
                    rowSums((X[,All_Covariates %in% combs[[a]] ]%*%W.prop[[a]])^2)) 
    return(PI.mMSE / sum(PI.mMSE)) } )   #
  
  PI.mMSE<-matrix(unlist(PI.mMSE),nrow = n,byrow = FALSE)
  pjoin.mMSE<-rowWeightedMeans(PI.mMSE,w=alpha)
  
  for (i in 1:length(r2)) 
  {
    ## mVc
    idx.mVc <- sample(1:n, r2[i]-r1, T, pjoin.mVc ) # ##
    
    fit.mVc<-lapply(1:length(combs),function(a){
      x.mVc <- X[c(idx.mVc, idx.prop),All_Covariates %in% combs[[a]] ] #
      y.mVc <- Y[c(idx.mVc, idx.prop)]
      pinv.mVc <- c(1 / pjoin.mVc[idx.mVc], pinv.prop)
      return(getMLE(x=x.mVc, y=y.mVc, w=pinv.mVc))
    })
    
    Sample.mVc[[i]]<-cbind(r2=r2[i],Y=Y[c(idx.mVc, idx.prop)],
                           SP=1/c(1 / pjoin.mVc[idx.mVc], pinv.prop),
                           X[c(idx.mVc, idx.prop),])
    
    for (j in 1:length(combs))
    {
      ru <- length(Y[c(idx.mVc, idx.prop)])
      beta.mVc[[j]][i,] <- fit.mVc[[j]]$par
      Bias_mVc[[j]][i,]<-Cordeiro(XData=X[c(idx.mVc, idx.prop),All_Covariates %in% combs[[j]] ],
                                  With_bias = beta.mVc[[j]][i,])
      
      pi <- invlogit(X[c(idx.mVc, idx.prop),All_Covariates %in% combs[[j]] ] %*% beta.mVc[[j]][i,]) #
      W <- diag(as.vector(pi*(1-pi)*c(1 / pjoin.mVc[idx.mVc], pinv.prop))) #
      Mx <- (t(X[c(idx.mVc, idx.prop),All_Covariates %in% combs[[j]] ]) %*% 
               W %*% X[c(idx.mVc, idx.prop),All_Covariates %in% combs[[j]] ])
      Mx <- solve(Mx) #
      Middle <- diag(((as.vector(Y[c(idx.mVc, idx.prop)])-as.vector(pi))*
                        as.vector(c(1 / pjoin.mVc[idx.mVc], pinv.prop)))^2) #
      
      V_Temp<-(t(X[c(idx.mVc, idx.prop),All_Covariates %in% combs[[j]] ]) %*% Middle %*%
                 X[c(idx.mVc, idx.prop),All_Covariates %in% combs[[j]] ]) #
      V_Final<-Mx %*% V_Temp %*% Mx #
      
      Utility_mVc[[j]][i,]<-cbind(r2[i],tr(V_Final),det(solve(V_Final))) #
    }
    
    ## mMSE
    idx.mMSE <- sample(1:n, r2[i]-r1, T, pjoin.mMSE)
    
    fit.mMSE<-lapply(1:length(combs),function(a){
      x.mMSE <- X[c(idx.mMSE, idx.prop),All_Covariates %in% combs[[a]] ] #
      y.mMSE <- Y[c(idx.mMSE, idx.prop)]
      pinv.mMSE <- c(1 / pjoin.mMSE[idx.mMSE], pinv.prop)
      return(getMLE(x=x.mMSE, y=y.mMSE, w=pinv.mMSE)$par)
    })
    
    Sample.mMSE[[i]]<-cbind(r2=r2[i],Y=Y[c(idx.mMSE, idx.prop)],
                            SP=1/c(1 / pjoin.mMSE[idx.mMSE], pinv.prop),
                            X[c(idx.mMSE, idx.prop),])
    
    for (j in 1:length(combs)) 
    {
      ru <- length(Y[c(idx.mMSE, idx.prop)])
      beta.mMSE[[j]][i,] <- fit.mMSE[[j]]
      Bias_mMSE[[j]][i,]<-Cordeiro(XData=X[c(idx.mMSE, idx.prop),All_Covariates %in% combs[[j]] ],
                                   With_bias = beta.mMSE[[j]][i,])
      
      pi <- invlogit(X[c(idx.mMSE, idx.prop),All_Covariates %in% combs[[j]] ] %*% beta.mMSE[[j]][i,]) #
      W <- diag(as.vector(pi*(1-pi)*c(1 / pjoin.mMSE[idx.mMSE], pinv.prop))) #
      Mx <- (t(X[c(idx.mMSE, idx.prop),All_Covariates %in% combs[[j]] ]) %*% 
               W %*% X[c(idx.mMSE, idx.prop),All_Covariates %in% combs[[j]] ])
      Mx <- solve(Mx) #
      Middle <- diag(((as.vector(Y[c(idx.mMSE, idx.prop)])-as.vector(pi))*
                        as.vector(c(1 / pjoin.mMSE[idx.mMSE], pinv.prop)))^2) #
      
      V_Temp<-(t(X[c(idx.mMSE, idx.prop),All_Covariates %in% combs[[j]] ]) %*% Middle %*%
                 X[c(idx.mMSE, idx.prop),All_Covariates %in% combs[[j]] ]) #
      V_Final<-Mx %*% V_Temp %*% Mx #
      
      Utility_mMSE[[j]][i,]<-cbind(r2[i],tr(V_Final),det(solve(V_Final))) #
    }
  }
  
  Sample_mVc<-do.call(rbind,Sample.mVc) #
  Sample_mMSE<-do.call(rbind,Sample.mMSE) #
  opt<-list()
  
  for (i in 1:length(combs)) 
  {
    opt[[i]] <- list("Est_Theta_mMSE"=cbind(r2,beta.mMSE[[i]]),"Utility_mMSE"=Utility_mMSE[[i]],
                     "Bias_mMSE"=cbind(r2,Bias_mMSE[[i]]),
                     "Est_Theta_mVc"=cbind(r2,beta.mVc[[i]]),"Utility_mVc"=Utility_mVc[[i]],
                     "Bias_mVc"=cbind(r2,Bias_mVc[[i]]))
    if(anyNA(opt[[i]]$Est_Theta_mMSE) || anyNA(opt[[i]]$Est_Theta_mVc) )
    {
      stop("There are NA or NaN values")
    }
    
    if(anyNA(opt[[i]]$Utility_mMSE) || anyNA(opt[[i]]$Utility_mVc) )
    {
      stop("There are NA or NaN values")
    }
    
    if(anyNA(opt[[i]]$Bias_mMSE) || anyNA(opt[[i]]$Bias_mVc) )
    {
      stop("There are NA or NaN values")
    }
  }
  if(anyNA(opt$Sample_mMSE) || anyNA(opt$Sample_mVc) )
  {
    stop("There are NA or NaN values")
  }
  return(list("opt"=opt,"Sample_mMSE"=Sample_mMSE,"Sample_mVc"=Sample_mVc))
}
