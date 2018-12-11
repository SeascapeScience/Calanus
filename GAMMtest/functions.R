
#############Fonction pour sauvegarder les GAMS##################

gam.table<- function(mod=mod, Region, Stock){
  library(mgcv)
  
  p.val<-as.data.frame(summary(mod)$s.table[, 4]) 
  colnames(p.val) <- "P"
  p.val$var<-rownames(p.val)
  
  
  if(length(summary(mod)$p.coeff)==1){
    intercept<-round(summary(mod)$p.coeff,4)
    intercept.p<-summary(mod)$p.pv
    if(intercept.p > 0.05){
      star <- "ns"
    }
    if(intercept.p < 0.05 & intercept.p > 0.01){
      star <- "*"
    } 
    if(intercept.p < 0.01 & intercept.p > 0.001){
      star <- "**"
    }
    if(intercept.p < 0.001){
      star <- "***"
    }
    star
    
    int<-paste(intercept, star, sep="")
    int
  }
  
  if(length(summary(mod)$p.coeff)==2){
    interceptD<-round(summary(mod)$p.coeff[1],4)
    interceptD.p<-summary(mod)$p.pv[1]
    interceptN<-round(summary(mod)$p.coeff[2],4)
    interceptN.p<-summary(mod)$p.pv[2]
    
    if(interceptD.p > 0.05){
      starD <- "ns"
    }
    if(interceptD.p < 0.05 & interceptD.p > 0.01){
      starD <- "*"
    } 
    if(interceptD.p < 0.01 & interceptD.p > 0.001){
      starD <- "**"
    }
    if(interceptD.p < 0.001){
      starD <- "***"
    }
    starD
    
    intD<-paste(interceptD, starD, sep="")
    intD
    
    if(interceptN.p > 0.05){
      starN <- "ns"
    }
    if(interceptN.p < 0.05 & interceptN.p > 0.01){
      starN <- "*"
    } 
    if(interceptN.p < 0.01 & interceptN.p > 0.001){
      starN <- "**"
    }
    if(interceptN.p < 0.001){
      starN <- "***"
    }
    starN
    
    intN<-paste(interceptN, starN, sep="")
    intN
    int <- paste("D:",intD," N:",intN, sep="")
  }
  
  GCV <- round(summary(mod)$sp.criterion,7)
  nvar<-summary(mod)$m
  DEV<-round(summary(mod)$dev.expl *100,2)
  R2<-round(summary(mod)$r.sq,2)
  
  
  
  library(plyr)
  
  pval<-ddply(p.val, "var", function(p.val){
    if(p.val$P > 0.05){
      star2 <- "ns"
    }
    if(p.val$P < 0.05 & p.val$P > 0.01){
      star2 <- "*"
    } 
    if(p.val$P < 0.01 & p.val$P > 0.001){
      star2 <- "**"
    }
    if(p.val$P < 0.001){
      star2 <- "***"
    }
    data.frame(star=star2)
  }
  )
  variables <- vector()
  for(i in 1:nrow(pval)){
    variables[i] <-paste(pval[i,1], pval[i,2], sep="")
  }
  formul <-paste(variables, collapse=" + ")
  formula <-paste("Var", formul, sep=" ~ ")
  
  
  N<-summary(mod)$n
  table<-as.data.frame(cbind(Stock,Region, formula, N, int, GCV,DEV,R2))
  table
}


