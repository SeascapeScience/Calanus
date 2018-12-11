mclapply.hack <- function(...) {
  library(parallel)
  ## Create a cluster
  ## ... How many workers do you need?
  ## ... N.B. list(...)[[1]] returns the first 
  ##          argument passed to the function. In
  ##          this case it is the list to iterate over
  size.of.list <- length(list(...)[[1]])
  cl <- makeCluster(min(size.of.list, detectCores())-2)#j'enlève un core pour pouvoir fonctionner ailleurs
  
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
################fonctions pour modèles mixte##############

automate_gamm4 <- function(RVar,family,Variables,Data,nb_simu,nvar){
  
  res <- mclapply.hack(1:nb_simu,function(w){
    cat("w=",w,"\n")#permet d'inscrire du texte pour te dire oÃ¹ en est ton script
    
    if(family=="gaussian"| family=="Gamma(link=log)"){
      samp <- sample(1:dim(Data)[1],round(dim(Data)[1]*70/100),replace=FALSE)
      dataC<- Data[samp,]; dataV<- Data[-samp,]#####data calibaration et data  Validation
    }
    
    if(family=="binomial"){ 
      Perc_nb_zero <- round(length(which(Data[,RVar]==0))/dim(Data)[1],2)#determination de la prÃ©valence de l'espÃ¨ce dans le cas d'une matrice binomiale
      Perc_nb_one  <- round(length(which(Data[,RVar]==1))/dim(Data)[1],2)# permet de ne pas créer de biais dans l'échantillonnage aléatoire
      
      samp_zero <- sample(which(Data[,RVar]==0),round(dim(Data)[1]*70/100)*Perc_nb_zero,replace=FALSE)
      samp_one  <- sample(which(Data[,RVar]==1),round(dim(Data)[1]*70/100)*Perc_nb_one ,replace=FALSE)
      
      sample_prev_ok <- c(samp_zero,samp_one)
      dataC<- Data[sample_prev_ok,]; dataV<- Data[-sample_prev_ok,]#data calibration et validation avec la prÃ©valence ok
    }
    
    global_formule <- paste(RVar,Variables,sep="")
    mS<-uGamm(as.formula(global_formule), random =~(1|fYear),data= dataC,family=family, lme4 = TRUE)
    
    cat("Computing dredge function ...","\n")
    d1 <- dredge(mS,rank="AIC",m.lim=c(2,10))#dredge permet la selection des variables en comparant les AIC 
    data_mod <- d1[which(d1[,"delta"]<=2)] #2e valeur est le nimbre de variables du modèle #selectionne seulement les modÃ¨les pour delta < 2
    
    cat("Best model selection ...","\n")
    best_mod<-var_best_mod<- list()
    for(i in 1: dim(data_mod)[1]){
      var_select <- colnames(data_mod)[which(data_mod[i,]=="+")]
      formule_C <- paste(RVar,"~",var_select[1],sep="")
      
      for(j in 2:length(var_select)) formule_C= paste(formule_C,var_select[j],sep="+")
      var_best_mod[[i]]  <- var_select
      best_mod[[i]] <- formule_C
    } # end of i/ cette Ã©tape permet de retranscrire le tableau de sortie de dredge en formule
    # end of i/ cette Ã©tape permet de retranscrire le tableau de sortie de dredge en formule
    
    ### Gam validation
    dataF <- dataV[,RVar]; dataG <- dataV[,c(1:nvar)]
    data_end <- cbind(dataF,dataG)
    colnames(data_end)[1] <- RVar # je recombine un tableau avec seulement les colonnes qui m'interessent
    
    cat("Predicting test ...","\n")
    
    if(family=="binomial"){
      vect_TSS <- vector()
      for(z in 1:length(best_mod)){
        mC <- gamm4(as.formula(best_mod[[z]]), random=~(1|fYear),data= dataC,family=family)
        Pred_mC <- predict(mC$gam,data_end,type="response")
        data_valid <- cbind(plotID= seq(1,dim(data_end)[1],1),Observed=data_end[,1],Predicted=Pred_mC)
        Tresh_opt <- optimal.thresholds(data_valid,na.rm=TRUE,opt.methods=3)[1,2]#option du package presence absence
        
        CMX1 <- cmx(data_valid, na.rm=TRUE, threshold=Tresh_opt)
        sens1 <- PresenceAbsence::sensitivity(CMX1, st.dev=FALSE)
        spec1 <-PresenceAbsence::specificity(CMX1, st.dev=FALSE)
        vect_TSS[z] <- (sens1+spec1)-1     
        
      }# end of z
      cat("Computing resultats ...","\n")
      if(length(vect_TSS)>1){
        res_var<- var_best_mod[[which.max(vect_TSS)]]# si plusieurs modÃ¨les, choix du meilleur TSS
      } else {
        res_var<- var_best_mod[[1]]
      } # end of if
      return(res_var)  
    } # end of if 
    
    if(family=="gaussian"| family=="Gamma(link=log)"){
      
      cor_data <- vector()
      for(z in 1:length(best_mod)){
        mC <- gamm4(as.formula(best_mod[[z]]), random=~(1|fYear), data= dataC,family=family)
        Pred_mC <- predict(mC$gam,data_end,type="response")
        cor_data[z]  <- cor.test(Pred_mC,data_end[,RVar], method="spearman")$estimate#donne le coef 
      } # end of z
      cat("Computing resultats ...","\n")
      if(length(cor_data)>1){
        res_var<- var_best_mod[[which.max(cor_data)]]# si plusieurs modeles choix du plus fort coef se spearman
      } else {
        res_var<- var_best_mod[[1]]
      } # end of if
      return(res_var)  
      
    } # end of if
    
    
    
  })
  
  return(res)
  
} # end of function



mumin_fct_gamm4 = function(RVar,family,Variables,nb_simu,Data,nvar,Sp, step){
  library(mgcv)
  library(plyr)
  library(PresenceAbsence)
  library(parallel)
  library(MuMIn)
  library(doBy)
  library(gamm4)
  
  inter = automate_gamm4(RVar=RVar,family=family,Variables=Variables,Data=Data,nb_simu=nb_simu,nvar=nvar)
  
  nb_var_mean = round(mean(sapply(inter,length)))#dÃ©termine le nombre moyen de selection de chaque variable
  selected_var = names(sort(table(unlist(inter))*100/nb_simu,decreasing=TRUE)[1:nb_var_mean])#determine les variables selectionnees sur l ensemble des iterations en tenant compte du nb total de vairables selectionnees Ã  chauqe simu
  
  formule_end <- paste(RVar,"~",selected_var[1],sep="")#Ã©dite la formule globale du gam selectionnÃ© sur les x itÃ©rations
  for(j in 2:length(selected_var)) formule_end= paste(formule_end,selected_var [j],sep="+")
  
  mC <- gamm4(as.formula(formule_end), random=~(1|fYear), data=Data,family=family)#refait tourner le gam sur 100% des data
  
  eval(parse(text=paste("Resmod_",RVar,Sp,step,"= mC",sep="")))
  eval(parse(text=paste("save(Resmod_",RVar,Sp,step,",file='Resmod_",RVar,Sp,step,"')",sep="")))#sauvegarde ton objet gam
  
  
  tiff(file=paste("Resmod_",RVar,Sp,step,".tiff",sep=""))#sauvegarde ta sortie graphique 
  plot(mC$gam,page=1,shade=TRUE,all.terms=TRUE, scale=0)
  title(round(summary(mC$gam)$r.sq,digits=3))
  dev.off()
  
}

