#' Function to write JABBA model in JAGS
#'
#' Writes JAGS code for JABBA into a temporary directory
#' @param jbinput JABBA input object from build_jabba()
#' @export

jabba2jags = function(jbinput, dir){
  
  # JAGS MODEL Standard
  sink(file.path(dir,"JABBA.jags"))
  cat("

    model {

    # Prior specifications
    eps <- 0.0000000000000000000000000000000001 # small constant")

if ("index_type" %in% names(jbinput$settings) & length(grep("absolute", jbinput$settings$index_type, ignore.case = TRUE)) > 0 ){
  if(jbinput$settings$nran.q == 1){
    cat("
  #Camera effective radius prior
    rad ~ dlnorm(rad.pr[1], pow(rad.pr[2],-2))T(7.5,60.6)

    q[1] ~ dunif(q_bounds[1], q_bounds[2])
    q[2] <- a.grid/(rad*rad*3.14159) 

  ", append = TRUE)
  }
  if(jbinput$settings$nran.q == 2){
cat("
   #Camera effective radius prior
    rad ~ dlnorm(rad.pr[1], pow(rad.pr[2],-2))T(7.5,60.6)

    q[1] ~ dunif(q_bounds[1], q_bounds[2])
    q[2] ~ dunif(q_bounds[1], q_bounds[2])
    q[3] <- a.grid/(rad*rad*3.14159) 
  
  ", append = TRUE)
  }

} else{
    cat("
    #Catchability coefficient
    for(i in 1:nq)
    {
      q[i] ~ dunif(q_bounds[1], q_bounds[2])
    }
    ", append = TRUE)
}
  
  cat("
    # Process variance prior
    isigma2.est ~ dgamma(igamma[1],igamma[2])

    # Carrying Capacity SB0
    K ~ dlnorm(log(K.pr[1]),pow(K.pr[2], -2))

    # informative priors for Hmsy as a function of r
    r ~ dlnorm(log(r.pr[1]),pow(r.pr[2],-2))

    ", append = TRUE)
  
  if(jbinput$settings$model.id==4){
    cat("
      # Shape m prior
      m ~ dlnorm(log(mu.m),pow(m.CV,-2))
      ",append=TRUE)
  }else{ cat("
      m <- mu.m
    ",append=TRUE)}
  
  if(jbinput$settings$psi.dist =="beta"){
    cat("
      # Beta Prior for Biomass depletion at the start (deteministic)
      psi ~ dbeta(psi.pr[1],psi.pr[2])
      ",append=TRUE)
  } else {
    cat("
      # Lognormal for Biomass depletion at the start (deteministic)
      psi ~ dlnorm(log(psi.pr[1]),pow(psi.pr[2],-2)) #I(0.1,1.1)
      ",append=TRUE)
  }
  
  if(jbinput$settings$sigma.proc==TRUE){
    cat("
      # Process variance
      isigma2 <- isigma2.est
      sigma2 <- pow(isigma2,-1)
      sigma <- sqrt(sigma2)
      fakesigma.fixed <- sigma.fixed # Prevent unused variable error msg
      ",append=TRUE)
  }else{ cat("
      isigma2 <- pow(sigma.fixed+eps,-2)
           sigma2 <- pow(isigma2,-1)
           sigma <- sqrt(sigma2)

           ",append=TRUE)}
  
  if(jbinput$settings$sigma.est==TRUE){
    if(jbinput$settings$nsig.off.ind == 0){
      cat("
      # Observation variance
      for(i in 1:(nvar))  #JS added -1 #MO removed
      {
      # Observation error
      itau2[i]~ dgamma(0.2,0.1)   #These are (0.001,0.001) OR (0.2,0.1) in previous assessment
      tau2[i] <- 1/itau2[i]
      }
  
      for(i in 1:(nvar)) #JS added -1 #MO changed from nI-1
      {
      for(t in 1:N)
      {
      var.obs[t,i] <- SE2[t,i]+tau2[sets.var[i]] 
      ivar.obs[t,i] <- (cpue_lambda[i]*cpue_lambda[i])/var.obs[t,i]
      # note total observation error (TOE)
      TOE[t,i] <- sqrt(var.obs[t,i]) # Total observation variance

      }}

      ",append=TRUE)

    }else if(jbinput$settings$nsig.off.ind == 1){
      cat("
      # Observation variance
      for(i in 1:(nvar))  #JS added -1 #MO removed
      {
      # Observation error
      itau2[i]~ dgamma(0.2,0.1)   #These are (0.001,0.001) OR (0.2,0.1) in previous assessment
      tau2[i] <- 1/itau2[i]
      }
  
      for(i in 1:(nvar)) #JS added -1 #MO changed from nI-1
      {
      for(t in 1:N)
      {
      var.obs[t,i] <- SE2[t,i]+tau2[sets.var[i]] 
      ivar.obs[t,i] <- (cpue_lambda[i]*cpue_lambda[i])/var.obs[t,i]
      # note total observation error (TOE)
      TOE[t,i] <- sqrt(var.obs[t,i]) # Total observation variance

      }}

     for(t in 1:N)  #MO added. This is assuming that the last index is the one that obs error is not being estimated for
      {
      var.obs[t,nq] <- SE2[t,nq]
      ivar.obs[t,nq] <- (cpue_lambda[nq]*cpue_lambda[nq])/var.obs[t,nq]
      # note total observation error (TOE)
      TOE[t,nq] <- sqrt(var.obs[t,nq]) # Total observation variance
      }
      ",append=TRUE)

    }
    
  }else{ cat("
      # Observation variance
           for(i in 1:(nvar)) #MO removed
           {
           # Observation error
           itau2[i]~ dgamma(0.2,0.1)
           tau2[i] <- 1/itau2[i]
           }

           for(i in 1:(nvar)) #MO changed from nI-1
           {
           for(t in 1:N)
           {
           var.obs[t,i] <- SE2[t,i] # drop tau2
           fake.tau[t,i] <- tau2[sets.var[i]]

           ivar.obs[t,i] <- (cpue_lambda[i]*cpue_lambda[i])/var.obs[t,i]  #JS added CPUE_lambda (cpue_lambda[i]*cpue_lambda[i])
           # note total observation error (TOE)
           TOE[t,i] <- sqrt(var.obs[t,i])

           }}

           ",append=TRUE)}
  
  # Run standard JABBA
  if(jbinput$settings$add.catch.CV==FALSE){
    cat("
  for(t in 1:N){
  estC[t] <- TC[t]
  }
",append=TRUE)} else if(jbinput$settings$catch.error=="random"){
  cat("
  for(t in 1:N){
      estC[t] ~ dlnorm(log(TC[t]),pow(CV.C[t],-2))
  }

",append=TRUE)} else if(jbinput$settings$catch.error=="deep7"){

  cat("
  #JS added uniform component to early years
  #Unrep is ~1/2 of catch, which had +/- 40% in previous assessment - reduce to 0.8,1.2? 
  #40% looks more consistent regarding 95% CIs on plot 
  for(t in 1:55){
      estC[t] ~ dunif((TC[t]*catch.adj[1]),(TC[t]*catch.adj[2]))
  }

  for(t in 56:N){
      estC[t] ~ dlnorm(log(TC[t]),pow(CV.C[t],-2))
  }

  ", append=TRUE)

} else {
  cat("
  for(t in 1:N){
  cdev[t] ~ dnorm(0,pow(CV.C[t],-2))
  cpos[t] <- sqrt(pow(cdev[t],2))
  estC[t] <- exp(log(TC[t])+cpos[t]) 
  }
",append=TRUE)}  
  cat("

    #Process equation base prediction
    Pmean[1] <- log(psi) 

   # penB[1]  <- ifelse( P.mean2.base[1]<P_bound[1],log(K* P.mean2.base[1])-log(K*P_bound[1]),ifelse( P.mean2.base[1]>P_bound[2],log(K* P.mean2.base[1])-log(K*P_bound[2]),0)) # penalty if Pmean is outside viable biomass 
    penBK[1] <- 0

    # Process equation base prediction
    for (t in 2:(N+1))
    {
    Pmean[t] <- ifelse(Pmean[t-1] > Plim,
    log(max(Pmean[t-1] +  r/(m-1)*Pmean[t-1]*(1-pow(Pmean[t-1],m-1)) - estC[t-1]/K,0.001)),
    log(max(Pmean[t-1] +  r/(m-1)*Pmean[t-1]*(1-pow(Pmean[t-1],m-1))*Pmean[t-1]*slope.HS - estC[t-1]/K,0.001)))
    }

   #Initialize uncorrelated residuals
   for (t in 1:(N+1))
    {
   P.base[t] <- dlnorm(P.mean2.base,isigma2)   
   log.resid.base[t]<-log(P.base[t]-log(Pmean[t]))
    }

   #Initialize prediction that accounts for autocorrelation
    P.mean2.base[1]<-Pmean[1] + phi * log.resid.0
   
    for (t in 2:(N+1))
    {
    P.mean2.base[t] <- Pmean[t] + phi * log.resid.base[t-1]
    }


  # P.mean1 <- Pmean[1:N]   ##JS added jags is particular about indexing so had to add these
  # P <- P.mean2.base[1:(N+1)]  ##JS added jags is particular about indexing so had to add these
  # log.resid<-log.resid.base[1:N]  ##JS added jags is particular about indexing so had to add these



    phi ~ dnorm(0,1.0E-4) I(-1,1)   #serial autocorrelation coefficient  
    log.resid.0 ~ dnorm(0,tau.red) T(-5,5)
    tau.red <- isigma2 * (1-phi*phi)  #this is precision of total biomass variation that includes white noise and correlated error
    sigma.red <- 1 / sqrt(tau.red)
    
   # for (t in 2:N)
   # {
   # penB[t]  <- ifelse(P.mean2.base[t]<(P_bound[1]),log(K*P.mean2.base[t])-log(K*(P_bound[1])),ifelse(P.mean2.base[t]>P_bound[2],log(K*P.mean2.base[t])-log(K*(P_bound[2])),0)) # penalty if Pmean is outside viable biomass
    # Depletion prior
   # penBK[t] <- ifelse(b.yr[t] < 1,0,log(ifelse(b.pr[4]<1,P.mean2.base[t],ifelse(b.pr[4]>1,HtoHmsy[t],BtoBmsy[t])))-log(b.pr[1]))
   # }


   #Previous process error incorporation
   # iPV[1] <- ifelse(1<(stI),10000,isigma2) # inverse process variance
   # P[1] ~ dlnorm(Pmean[1],iPV[1]) # set to small noise instead of isigma2
    

   # iPV[t] <- ifelse(t<(stI),10000,isigma2) # inverse process variance
   # P[t] ~ dlnorm(Pmean[t],iPV[t])



    # Process error deviation
    for(t in 1:N){
    Proc.Dev[t] <- log(P.mean2.base[t]*K)-log(exp(Pmean[t])*K)}

    # Enforce soft penalties on bounds for P
    for(t in 1:N){
    pen.P[t] ~ dnorm(penB[t],1000) # enforce penalty with CV = 0.1
    pen.bk[t] ~ dnorm(penBK[t],pow(b.pr[2],-2))
    }

    Hmsy <- r*pow(m-1,-1)*(1-1/m)
    SBmsy_K <- (m)^(-1/(m-1))
    SBmsy <- SBmsy_K*K
    MSY <- SBmsy*Hmsy
    
    for (t in 1:N)
    {
    
    H[t] <- TC[t]/SB[t]
    HtoHmsy[t] <- H[t]/(Hmsy)
    }


    
    
    for (t in 1:(N+1)) # One step ahead biomass
    {
    SB[t] <- K*P.mean2.base[t]   
    BtoBmsy[t] <- SB[t]/SBmsy
    }

    for (t in 1:N)  # Overfishing status given state-dependent rule
    {
    Overfishing_ind[t]<-ifelse(BtoBmsy[t]>0.866, H[t]/(Hmsy),  H[t]/((Hmsy*SB[t])/(0.866*SBmsy)))
    }


    ",append=TRUE)

  if(jbinput$settings$CatchOnly==FALSE){
    if (length(grep("absolute", jbinput$settings$index_type, ignore.case = TRUE)) > 0){

     cat("
     # Observation equation in related to EB

    for(i in 1:nran.q)
    {
    for (t in 1:N)
    {
    Imean[t,i] <- log(q[sets.q[i]]*P.mean2.base[t]*K);
    I[t,i] ~ dlnorm(Imean[t,i],(ivar.obs[t,i]));
    CPUE[t,i] ~ dlnorm(Imean[t,i],(ivar.obs[t,i]))   ####q[[i]]*P[t]*SB0*EBtoSB[t,i]
    Ihat[t,i]  <- exp(Imean[t,i])

    }}",append=TRUE)
     if (length(grep("absolute", jbinput$settings$index_type, ignore.case = TRUE)) > 0 & jbinput$settings$nran.q == 1){
  cat("
      ## added for BFISH index
      for (t in 1:N){ 
        
        Imean[t,2] <- log(P.mean2.base[t] * K/((q[sets.q[nq]])*n.grid))
        survey_precision[t] <- (s_lambda*s_lambda)/((ivar.obs[t,2]))  #( ivar.obs[t,2]) #MO changed SE2 to ivar.obs
        I[t, 2] ~ dlnorm(Imean[t,2], survey_precision[t])
        CPUE[t,2] <- P.mean2.base[t]*K/(q[sets.q[nq]]*n.grid)
        Ihat[t,2]  <- exp(Imean[t,2])

        }

  ",append=TRUE)
  }else if(length(grep("absolute", jbinput$settings$index_type, ignore.case = TRUE)) > 0 & jbinput$settings$nran.q == 2){
      cat("
      
      ## added for BFISH index
      for (t in 1:N){ 
        
        Imean[t,3] <- log(PP.mean2.base[t] * K/((q[sets.q[nq]])*n.grid))
        survey_precision[t] <- (s_lambda*s_lambda)/((ivar.obs[t,3]))  #( ivar.obs[t,3])
        I[t, 3] ~ dlnorm(Imean[t,3], survey_precision[t])
        CPUE[t,3] <- P.mean2.base[t]*K/(q[sets.q[nq]]*n.grid)
        Ihat[t,3]  <- exp(Imean[t,3])

        }
      
      ",append=TRUE)}

  }else{
    cat("
    # Observation equation in related to EB

    for(i in 1:(nI))
    {
    for (t in 1:N)
    {
    Imean[t,i] <- log(q[sets.q[i]]*P.mean2.base[t]*K);
    I[t,i] ~ dlnorm(Imean[t,i],ivar.obs[t,i]);
    CPUE[t,i] ~ dlnorm(Imean[t,i],(ivar.obs[t,i]))   ####q[[i]]*P[t]*SB0*EBtoSB[t,i]
    Ihat[t,i]  <- exp(Imean[t,i])

    }}",append=TRUE)

  }

    }
  
    if(jbinput$settings$CatchOnly==TRUE){

      if(length(grep("absolute", jbinput$settings$index_type, ignore.case = TRUE)) > 0){
        cat("
            # Observation equation in related to EB
              one <- sets.q
              for(i in 1:(nran.q))
              {
              for (t in 1:N)
              {
              Imean[t,i] <- log(pow(P.mean2.base[1],-1)*P.mean2.base[t]);
              I[t,i] ~ dlnorm(Imean[t,i],(ivar.obs[t,i]));
              CPUE[t,i] ~ dlnorm(Imean[t,i],(ivar.obs[t,i]))   ####q[[i]]*P[t]*SB0*EBtoSB[t,i]
              Ihat[t,i]  <- exp(Imean[t,i])

              }}
        ",append=TRUE)

        if (length(grep("absolute", jbinput$settings$index_type, ignore.case = TRUE)) > 0 & jbinput$settings$nran.q == 1){

          cat("
          ## added for BFISH index
              for (t in 1:N){ 
              
              Imean[t,2] <- log(P.mean2.base[t] * K/((q[sets.q[nq]])*n.grid))
              survey_precision[t] <- (s_lambda*s_lambda)/(SE2[t,2])  #( ivar.obs[t,2]) #JS deleted sqrt of SE2
              I[t, 2] ~ dlnorm(Imean[t,2], survey_precision[t])
              CPUE[t,2] <- P.mean2.base[t]*K/(q[sets.q[nq]]*n.grid)
              Ihat[t,2]  <- exp(Imean[t,2])

              }
          
          ",append=TRUE)


        }else if(length(grep("absolute", jbinput$settings$index_type, ignore.case = TRUE)) > 0 & jbinput$settings$nran.q == 2){

          cat("
             ## added for BFISH index
              for (t in 1:N){ 
              
              Imean[t,3] <- log(P.mean2.base[t] * K/((q[sets.q[nq]])*n.grid))
              survey_precision[t] <- (s_lambda*s_lambda)/(SE2[t,3])  #( ivar.obs[t,3]) #JS deleted sqrt of SE2
              I[t, 3] ~ dlnorm(Imean[t,3], survey_precision[t])
              CPUE[t,3] <- P.mean2.base[t]*K/(q[sets.q[nq]]*n.grid)
              Ihat[t,3]  <- exp(Imean[t,3])

              }
          
          ",append=TRUE)


        }
      }else{
          cat("
           # Observation equation in related to EB
            one <- sets.q
            for(i in 1:(nI))
            {
            for (t in 1:N)
            {
            Imean[t,i] <- log(pow(P.mean2.base[1],-1)*P.mean2.base[t]);
            I[t,i] ~ dlnorm(Imean[t,i],(ivar.obs[t,i]));
            CPUE[t,i] ~ dlnorm(Imean[t,i],(ivar.obs[t,i]))   ####q[[i]]*P[t]*SB0*EBtoSB[t,i]
            Ihat[t,i]  <- exp(Imean[t,i])

            }}
          
          ",append=TRUE)
      }
} 
    
  cat("
      
    # Enforce soft penalty on K if < K_bounds >
    K.pen ~ dnorm(penK,1000) # enforce penalty
    penK  <- ifelse(K<(K_bounds[1]),log(K)-log(K_bounds[1]),ifelse(K>K_bounds[2],log(K)-log(K_bounds[2]),0)) # penalty if Pmean is outside viable biomass


    # Enforce soft penalty on process deviance if sigma.proc > 0.2
    proc.pen ~ dnorm(penProc,1000) # enforce penalty
    penProc  <- ifelse(sigma>sigmaproc_bound,log(sigma)-log(sigmaproc_bound),0)


    # Enforce soft penalty on observation error if sigma.obs > sigma_bound
    for(i in 1:(nvar)){    #JS - may need to add -1 here... #MO removed
    obs.pen[i] ~ dnorm(penObs[i],1000) # enforce penalty
    penObs[i]  <- ifelse(pow(tau2[i],0.5)>sigmaobs_bound,log(pow(tau2[i],0.5))-log(sigmaobs_bound),0)
    }

    ", append=TRUE)
  # Add Auxiliary observation likelihood
  if(jbinput$settings$Auxiliary==TRUE){

    if(jbinput$settings$auxiliary.sigma==TRUE){
      cat("
      # Observation variance
      for(i in 1:nAvar)
      {
      # Observation error
      ieta2[i]~ dgamma(0.001,0.001)
      eta2[i] <- 1/ieta2[i]
      }

      for(i in 1:nA)
      {
      for(t in 1:N)
      {
      varA.obs[t,i] <- A.SE2[t,i]+eta2[sets.varA[i]]
      ivarA.obs[t,i] <- 1/varA.obs[t,i]
      # note total auxiliary observation error (TAE)
      TAE[t,i] <- sqrt(varA.obs[t,i]) # Total observation variance
      }}
      ",append=TRUE) }
      else{ cat("
      # Observation variance
          
         
          
      for(i in 1:nA)
      {
      for(t in 1:N)
      {
      varA.obs[t,i] <- A.SE2[t,i] # Fixed only
      ivarA.obs[t,i] <- 1/varA.obs[t,i]
      # note total auxiliary observation error (TAE)
      TAE[t,i] <- sqrt(varA.obs[t,i]) # Total observation variance
      }}

           ",append=TRUE)}
    
      if(jbinput$settings$auxiliary.type=="effort"){
      cat("
  
      for(i in 1:nA){ 
      qA[i] ~ dunif(qA_bounds[1],qA_bounds[2])} 
      for(t in 1:N){
      Ax[t] <- H[t]}
      ",append=TRUE)}
      
    if(jbinput$settings$auxiliary.type=="z"){
      cat("
      for(i in 1:nA){
      qA[i] ~ dlnorm(0,pow(qA.cv,-2))T(qA_bounds[1],qA_bounds[2])
      }  
      for(t in 1:N){
      #Ax[t] <- -log(max(1-H[t],0.01))-log(1-max(1-Hmsy,0.01))
      Ax[t] <- max(1-(1-Hmsy)*(1-H[t]),0.001)
      }
      ",append=TRUE)}
    
      if(jbinput$settings$auxiliary.type=="f"){
      cat("
      for(i in 1:nA){
      qA[i] ~ dlnorm(0,pow(qA.cv,-2))T(qA_bounds[1],qA_bounds[2])
      }  
      for(t in 1:N){
      Ax[t] <- -log(1-H[t])
      }
      ",append=TRUE)}
     
    if(jbinput$settings$auxiliary.type=="ffmsy"){
      cat("
      for(i in 1:nA){
      qA[i] ~ dlnorm(0,pow(qA.cv,-2))T(qA_bounds[1],qA_bounds[2])
      }  
      for(t in 1:N){
      Ax[t] <- HtoHmsy[t]
      }
      ",append=TRUE)}
    
    if(jbinput$settings$auxiliary.type=="bbmsy"){
      cat("
      for(i in 1:nA){
      qA[i] ~ dlnorm(0,pow(qA.cv,-2))T(qA_bounds[1],qA_bounds[2])
      }  
      for(t in 1:N){
      Ax[t] <- BtoBmsy[t]
      }
      ",append=TRUE)}
       
    if(jbinput$settings$auxiliary.type=="bk"){
      cat("
      for(i in 1:nA){
      qA[i] ~ dlnorm(0,pow(qA.cv,-2))T(qA_bounds[1],qA_bounds[2])
      }  
      for(t in 1:N){
      Ax[t] <- P[t]
      }
      ",append=TRUE)}
    
    if(jbinput$settings$auxiliary.lag>0){
      cat("
      for(i in 1:nA){ 
      for(t in 1:A.lag){
      
      #Amean[t,i] <- log(qA[i])+log(mean(Ax[1:t]))
      Amean[t,i] <- log(qA[i])+log(mean(Ax[1]))
      
      A[t,i] ~ dlnorm(Amean[t,i],ivarA.obs[t,i])
      AUXI[t,i] ~ dlnorm(Amean[t,i],ivarA.obs[t,i]) 
      Ahat[t,i]  <- exp(Amean[t,i])}}
      
      for(i in 1:nA){
      for(t in (A.lag+1):N){
      Amean[t,i] <- log(qA[i])+log(mean(Ax[(t-A.lag):(t-1)]))
      #Amean[t,i] <- log(qA[i])+log(mean(Ax[(t-A.lag)]))
      
      A[t,i] ~ dlnorm(Amean[t,i],ivarA.obs[t,i])
      AUXI[t,i] ~ dlnorm(Amean[t,i],ivarA.obs[t,i]) 
      Ahat[t,i]  <- exp(Amean[t,i])
      }}
  
    ",append=TRUE)}
    # No Lag
    if(jbinput$settings$auxiliary.lag==0){
     cat("
      for(i in 1:nA){
      for(t in 1:N){
      Amean[t,i] <- log(qA[i])+log(mean(Ax[t]))
      A[t,i] ~ dlnorm(Amean[t,i],ivarA.obs[t,i])
      AUXI[t,i] ~ dlnorm(Amean[t,i],ivarA.obs[t,i]) 
      Ahat[t,i]  <- exp(Amean[t,i])
      }}
      
      ",append=TRUE)
    }
  } # End of auxiliary
  
  # PROJECTION
  if(jbinput$settings$projection==TRUE){
    cat("
      for(i in 1:nTAC){
      # Project first year into the future
      prPmean[1,i] <- ifelse(P[N] > Plim,
      log(max(P[N] +  Hmsy/(1-1/m)*P[N]*(1-pow(P[N],m-1)) - TC[N]/K,0.005)),
      log(max(P[N] +  Hmsy/(1-1/m)*P[N]*(1-pow(P[N],m-1))*4*P[N] - TC[N]/K,0.005)))
      prP[1,i] ~ dlnorm(prPmean[1,i],isigma2)
      # Project all following years
      for(t in 2:pyrs){
      prPmean[t,i] <- ifelse(prP[t-1,i] > Plim,
      log(max(prP[t-1,i] +  Hmsy/(1-1/m)*prP[t-1,i]*(1-pow(prP[t-1,i],m-1)) - TAC[t-1,i]/K,0.005)),
      log(max(prP[t-1,i] +  Hmsy/(1-1/m)*prP[t-1,i]*(1-pow(prP[t-1,i],m-1))*slope.HS*prP[t-1,i] - TAC[t-1,i]/K,0.005)))
      # process error (as monte-carlo simular)
      prP[t,i] ~ dlnorm(prPmean[t,i],isigma2)}
      for(t in 1:pyrs){
      prB[t,i] <- prP[t,i]*K
      prH[t,i] <- TAC[t,i]/prB[t,i]
      prHtoHmsy[t,i] <- prH[t,i]/Hmsy
      prBtoBmsy[t,i] <- prB[t,i]/SBmsy
      prOverfishing_ind[t,i]<-ifelse(prBtoBmsy[t,i]>0.866, prH[t,i]/(Hmsy),  prH[t,i]/((Hmsy*prB[t,i])/(0.866*SBmsy)))
      }}
      ",append=TRUE)} else {
        cat("
            #Prevent error for unused input
            fakeTAC <-  TAC
            fakepyrs <- pyrs
            fakenTAC <- nTAC
            fakenranq <- nran.q
            prHtoHmsy <- 1
            prP <- 1
            prBtoBmsy <- 1
            ", append=TRUE)}
  
  cat("
} # END OF MODEL
    ",append=TRUE,fill = TRUE)
  sink()
} # #END jabba JAGS
