###Function TSS_mclapply_ie  evaluate model fit by cross-validation
## Applies to binomial and gaussian family using default link function.
# Can be used with mix or standard models; gam_type="mixte" or gam_type="standard"#
#Works with GAM. Can be easily modified for GLM
#Run in parallel using ncore-2. If you are not running other software you can use ncore at line 49
#Data are splitted 70/30 for model calibration/validation


#arguments: 
#mcomp: name of the gam object followed by $gam for mixted models
#Data: data.frame used to fit the models. NA should be remove before model. Data should have the same number of observations as used in the model
#gam_type= "mixte" OR "standard", for GAMM OR GAM
#fam= family "binomial OR "gaussian"
#nrem= N type to repeat the function. Try a small number first such as 10 to see how much time it takes. Can be time consuming with gamm. Try to use a number that fits the number of cores of your computer.
#Sp= name of the response variable as named in Data.


#output: 1 single table for binomial and gaussian models. Is NA when not applicable to the type of model

#suffixe e: external validation, using 70% for model fittinh and 30% for model validation
#suffixe i: internal validation, using 100% for fitting and validation, resubstitution

#ObsPrev_abund: Prevalence (proportion of 0-1) observed for the binomial family OR  mean of the response variable for the gaussian family

#Fo binomial models:
#Tresh_opt_e_i: Mean of the optimal treshold to determine 0-1. Using method : 3	MaxSens+Spec	 maximizes (sensitivity+specificity)/2
#TSS_e_i: Mean of TSS= (sensitivity+specificity)-1 
#TSSsd_e_i: Standard deviation du TSS (internal valdiation is 0)
#AUC_e_i : Mean Area under the curve ROC

#For gaussian models
#corr_r_e_i: Mean Spearman correlation coefficient between prediction and observation
#p.corr_e_i: Mean p-value of the correlation coefficient
#intercept_e_i: Mean intercept of the linear regression between observation and prediction: optimal =0
#slope_e_i: Mean slope of the linear regression between observation and prediction, optimal =1
#r2_e_i: Mean R2 of the same linear regression
#r2sd_e_i: Mean R2 of the same linear regression



mclapply.hack <- function(...) {
  library(parallel)
  ## Create a cluster
  ## ... How many workers do you need?
  ## ... N.B. list(...)[[1]] returns the first 
  ##          argument passed to the function. In
  ##          this case it is the list to iterate over
  size.of.list <- length(list(...)[[1]])
  cl <- makeCluster(min(size.of.list, detectCores())-2)#removing to core to use other softwares at the same time
  
  ## Find out the names of the loaded packages 
  loaded.package.names <- c(
    ## Base packages
    sessionInfo()$basePkgs,
    ## Additional packages
    names( sessionInfo()$otherPkgs ))
  
  ## N.B. tryCatch() allows us to properly shut down the 
  ##      cluster if an error in our code halts execution
  ##      of the function. For details see: help(tryCatch)
  tryCatch( {
    
    ## Copy over all of the objects within scope to
    ## all clusters. 
    ## 
    ## The approach is as follows: Beginning with the 
    ## current environment, copy over all objects within
    ## the environment to all clusters, and then repeat
    ## the process with the parent environment. 
    ##
    this.env <- environment()
    while( identical( this.env, globalenv() ) == FALSE ) {
      clusterExport(cl,
                    ls(all.names=TRUE, env=this.env),
                    envir=this.env)
      this.env <- parent.env(environment())
    }
    ## repeat for the global environment
    clusterExport(cl,
                  ls(all.names=TRUE, env=globalenv()),
                  envir=globalenv())
    
    ## Load the libraries on all the clusters
    ## N.B. length(cl) returns the number of clusters
    parLapply( cl, 1:length(cl), function(xx){
      lapply(loaded.package.names, function(yy) {
        ## N.B. the character.only option of 
        ##      require() allows you to give the 
        ##      name of a package as a string. 
        require(yy , character.only=TRUE)})
    })
    
    ## Run the lapply in parallel 
    return( parLapply( cl, ...) )
  }, finally = {        
    ## Stop the cluster
    stopCluster(cl)
  })
}



TSS_mclapply_ie <-function(mcomp, Data, gam_type, fam, nrem,Sp){   
 
    library(PresenceAbsence)
    library(gamm4)
    library(verification)
    res <- mclapply.hack(1:nrem,function(w){
    ObsPrev_abund<-NA
    Tresh_opt<-NA
    Tresh_opti<-NA
    TSS<-NA
    TSSi<-NA
    AUC<-NA
    AUCi<-NA
    corr_coeff<-NA
    corr_coeffi<-NA
    p.corr<-NA 
    p.corri<-NA
    intercept<-NA
    intercepti<-NA
    slope<-NA
    slopei<-NA
    r2<-NA 
    r2i<-NA
    TSSsd_e_i<-NA
    rssd_e_i<-NA
    
    nRem <- round(0.3 * nrow(Data))
    sel <- sample(1:nrow(Data), nRem)
    dataV <- Data[sel,]
    dataC <- Data[-sel,]  
    
    if(gam_type=="mixte"){
      mod <- gamm4(summary(mcomp)$formula, random=~(1|fYear), data=dataC,family=fam)
      Pred_mod <- predict(mod$gam,dataV,type="response")
      Data_valid <- cbind(plotID= seq(1,dim(dataV)[1],1),Observed=dataV[,Sp],Predicted=Pred_mod)
      
      Pred_mod1 <- predict(mcomp,Data,type="response")
      Data_valid1 <- cbind(plotID= seq(1,dim(Data)[1],1),Observed=Data[,Sp],Predicted=Pred_mod1)
      
    }
    if(gam_type=="standard"){
      mod <- gam(summary(mcomp)$formula, data=dataC,family=fam)
      Pred_mod <- predict(mod,dataV,type="response")
      Data_valid <- cbind(plotID= seq(1,dim(dataV)[1],1),Observed=dataV[,Sp],Predicted=Pred_mod)
      
      Pred_mod1 <- predict(mcomp,Data,type="response")
      Data_valid1 <- cbind(plotID= seq(1,dim(Data)[1],1),Observed=Data[,Sp],Predicted=Pred_mod1)
    }
    
    
    
    if(fam=="binomial"){
      Tresh_opt <- optimal.thresholds(Data_valid,na.rm=TRUE,opt.methods=3)[1,2]
      CMX1 <- cmx(Data_valid, na.rm=TRUE, threshold=Tresh_opt)
      sens1 <- PresenceAbsence::sensitivity(CMX1, st.dev=FALSE)
      spec1 <- PresenceAbsence::specificity(CMX1, st.dev=FALSE)
      TSS <- (sens1+spec1)-1 
      AUC<-roc.area(Data_valid[,2], Pred_mod)$A
      
      
      Tresh_opti <- optimal.thresholds(Data_valid1,na.rm=TRUE,opt.methods=3)[1,2]
      CMX1b <- cmx(Data_valid1, na.rm=TRUE, threshold=Tresh_opti)
      sens1b <-  PresenceAbsence::sensitivity(CMX1b, st.dev=FALSE)
      spec1b <- PresenceAbsence::specificity(CMX1b, st.dev=FALSE)
      TSSi <- (sens1b+spec1b)-1 
      AUCi<-roc.area(Data_valid1[,2], Pred_mod1)$A
      
      
    }
    if(fam=="gaussian"){
      corr<-cor.test(Data_valid[,2],Data_valid[,3], method="spearman")
      p.corr<-corr[3]
      corr_coeff <-corr[4]$estimate
      lm1<-lm(Data_valid[,3]~Data_valid[,2])
      coeff_lm<-as.data.frame(t(lm1[1]))
      intercept<-as.data.frame(lm1$coefficient[1])[1,1]
      slope<-as.data.frame(lm1$coefficient[2])[1,1]
      r2 <- summary(lm1)[8]
      p.lm<-anova(lm1)[1,5]
      
      corri<-cor.test(Data_valid1[,2],Data_valid1[,3], method="spearman")
      p.corri<-corri[3]
      corr_coeffi <-corri[4]$estimate
      lm1b<-lm(Data_valid1[,3]~Data_valid1[,2])
      coeff_lmi<-as.data.frame(t(lm1b[1]))
      intercepti<-as.data.frame(lm1b$coefficient[1])[1,1]
      slopei<-as.data.frame(lm1b$coefficient[2])[1,1]
      r2i <- summary(lm1b)[8]
      p.lmi <-anova(lm1b)[1,5]
      
    }
    
    Ti<-cbind(Tresh_opt,Tresh_opti,TSS, TSSi,AUC,AUCi,corr_coeff,corr_coeffi, p.corr, p.corri,intercept,intercepti,slope, slopei, r2, r2i)
    Ti
  }
  )
  lres<-unlist(res)
  Tresh_opt_e_i<-paste(round(mean(lres[seq(1,16*nrem,16)]),2),round(mean(lres[seq(2,16*nrem,16)]),2), sep="/")
  TSS_e_i<-paste(round(mean(lres[seq(3,16*nrem,16)]),2),round(mean(lres[seq(4,16*nrem,16)]),2), sep="/")
  TSSsd_e_i<-paste(round(sd(lres[seq(3,16*nrem,16)]),2),round(sd(lres[seq(4,16*nrem,16)]),2), sep="/")
  AUC_e_i<-paste(round(mean(lres[seq(5,16*nrem,16)]),2),round(mean(lres[seq(6,16*nrem,16)]),2), sep="/")
  corr_r_e_i<-paste(round(mean(lres[seq(7,16*nrem,16)]),2),round(mean(lres[seq(8,16*nrem,16)]),2), sep="/")
  p.corr_e_i<-paste(round(mean(lres[seq(9,16*nrem,16)]),2),round(mean(lres[seq(10,16*nrem,16)]),2), sep="/")
  intercept_e_i<-paste(round(mean(lres[seq(11,16*nrem,16)]),2),round(mean(lres[seq(12,16*nrem,16)]),2), sep="/")
  slope_e_i<-paste(round(mean(lres[seq(13,16*nrem,16)]),2),round(mean(lres[seq(14,16*nrem,16)]),2), sep="/")
  r2_e_i<-paste(round(mean(lres[seq(15,16*nrem,16)]),2),round(mean(lres[seq(16,16*nrem,16)]),2), sep="/")
  r2sd_e_i<-paste(round(sd(lres[seq(15,16*nrem,16)]),2),round(sd(lres[seq(16,16*nrem,16)]),2), sep="/")
  ObsPrev_abund<-round(mean(Data[,Sp]),2)
  
  fres<-cbind(ObsPrev_abund,Tresh_opt_e_i,TSS_e_i,TSSsd_e_i,AUC_e_i,corr_r_e_i,p.corr_e_i,intercept_e_i, slope_e_i, r2_e_i, r2sd_e_i)
  fres
}

