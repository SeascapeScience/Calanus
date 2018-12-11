source('PMZA_Data_load.R')
#library(lattice)
#xyplot(Cfin_CI_CIV ~ Month|REGION, data=pmza)


panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
pairs(pmza[,c("T_NB","T0_50","S0_50","lDEPTH","Year")], lower.panel=panel.cor)

############Exploration###########
#relation entre STRAT et MONTH
#STRAT est enlev? de tous les mod?les
#library(lattice)
#xyplot(lSTRAT~Month|REGION, data=pmza)
source('E:/Projects/Scripts_R/functions.R')
library(doBy)
library(MuMIn)
library(parallel)


source('E:/Projects/Scripts_R/fct_selMODEL.R')
anyNA(pmza)
options(na.action="na.fail")
setwd("E:/Projects/ACCASP_FINAL/NEW_ACCASP_pheno/PMZA_Models/")

formul = "~ s(T_NB, k=4) + s(T0_50, k=5) + s(S0_50, k=5) + s(lDEPTH, k=5) + s(Month, k=6, by=REGION)"
dat2b<-pmza[,c("lCfin_CV_CVI", "T_NB", "T0_50", "S0_50", "lDEPTH", "Month", "REGION", "fYear")]
mumin_fct_gamm4(RVar="lCfin_CV_CVI",family="gaussian",Variables=formul,nb_simu=10,Data=dat2b,nvar=ncol(dat2b),Sp="", step="")


formul = "~ s(T_NB, k=5) + s(T0_50, k=5) + s(S0_50, k=5) + s(lDEPTH, k=5) + s(Month, k=8, by=REGION)"
dat3<-pmza[,c("Cgla_CI_CIII_occ", "T_NB", "T0_50", "S0_50", "lDEPTH", "Month", "REGION", "fYear")]
mumin_fct_gamm4(RVar="Cgla_CI_CIII_occ",family="binomial",Variables=formul,nb_simu=10,Data=dat3,nvar=ncol(dat3),Sp="", step="")

dat4<-pmza_cglac1_pres[,c("lCgla_CI_CIII", "T_NB", "T0_50", "S0_50", "lDEPTH", "Month", "REGION", "fYear")]
mumin_fct_gamm4(RVar="lCgla_CI_CIII",family="gaussian",Variables=formul,nb_simu=10,Data=dat4,nvar=ncol(dat4),Sp="", step="")

dat5<-pmza[,c("Chyp_CI_CIII_occ", "T_NB", "T0_50", "S0_50", "lDEPTH", "Month", "REGION", "fYear")]
mumin_fct_gamm4(RVar="Chyp_CI_CIII_occ",family="binomial",Variables=formul,nb_simu=10,Data=dat5,nvar=ncol(dat5),Sp="", step="")

dat6<-pmza_chypc1_pres[,c("lChyp_CI_CIII", "T_NB", "T0_50", "S0_50", "lDEPTH", "Month", "REGION", "fYear")]
mumin_fct_gamm4(RVar="lChyp_CI_CIII",family="gaussian",Variables=formul,nb_simu=10,Data=dat6,nvar=ncol(dat6),Sp="", step="")

dat7<-pmza[,c("Cgla_CIV_CVI_occ", "T_NB", "T0_50", "S0_50", "lDEPTH", "Month", "REGION", "fYear")]
mumin_fct_gamm4(RVar="Cgla_CIV_CVI_occ",family="binomial",Variables=formul,nb_simu=10,Data=dat7,nvar=ncol(dat7),Sp="", step="")

formul = "~ s(T_NB, k=4) + s(T0_50, k=5) + s(S0_50, k=5) + s(lDEPTH, k=5) + s(Month, k=6, by=REGION)"
dat8<-pmza_cglac4_pres[,c("lCgla_CIV_CVI", "T_NB", "T0_50", "S0_50", "lDEPTH", "Month", "REGION", "fYear")]
dat8<-subset(dat8, dat8$T_NB <10)
mumin_fct_gamm4(RVar="lCgla_CIV_CVI",family="gaussian",Variables=formul,nb_simu=10,Data=dat8,nvar=ncol(dat8),Sp="", step="")

formul = "~ s(T_NB, k=5) + s(T0_50, k=5) + s(S0_50, k=5) + s(lDEPTH, k=5) + s(Month, k=8, by=REGION)"
dat9<-pmza[,c("Chyp_CIV_CVI_occ", "T_NB", "T0_50", "S0_50", "lDEPTH", "Month", "REGION", "fYear")]
mumin_fct_gamm4(RVar="Chyp_CIV_CVI_occ",family="binomial",Variables=formul,nb_simu=10,Data=dat9,nvar=ncol(dat9),Sp="", step="")

formul = "~ s(T_NB, k=5) + s(T0_50, k=5) + s(S0_50, k=5) + s(lDEPTH, k=5) + s(Month, k=8, by=REGION)"
dat10<-pmza_chypc4_pres[,c("lChyp_CIV_CVI", "T_NB", "T0_50", "S0_50", "lDEPTH", "Month", "REGION", "fYear")]
mumin_fct_gamm4(RVar="lChyp_CIV_CVI",family="gaussian",Variables=formul,nb_simu=10,Data=dat10,nvar=ncol(dat10),Sp="", step="")

formul = "~ s(T_NB, k=5) + s(T0_50, k=5) + s(S0_50, k=5) + s(lDEPTH, k=5) + s(Month, k=8, by=REGION)"
dat11<-pmza[,c("Para_occ", "T_NB", "T0_50", "S0_50", "lDEPTH", "Month", "REGION", "fYear")]
mumin_fct_gamm4(RVar="Para_occ",family="binomial",Variables=formul,nb_simu=10,Data=dat11,nvar=ncol(dat11),Sp="", step="")

dat12<-pmza_para_pres[,c("lPara", "T_NB", "T0_50", "S0_50", "lDEPTH", "Month", "REGION", "fYear")]
mumin_fct_gamm4(RVar="lPara",family="gaussian",Variables=formul,nb_simu=10,Data=dat12,nvar=ncol(dat12),Sp="", step="")


####remove overfitting#####
formul = "~ s(T_NB, k=5) + s(T0_50, k=5) + s(S0_50, k=5) + s(lDEPTH, k=4) + s(Month, k=8, by=REGION)"
dat3<-pmza[,c("Cgla_CI_CIII_occ", "T_NB", "T0_50", "S0_50", "lDEPTH", "Month", "REGION", "fYear")]
mumin_fct_gamm4(RVar="Cgla_CI_CIII_occ",family="binomial",Variables=formul,nb_simu=10,Data=dat3,nvar=ncol(dat3),Sp="", step="klimit")


#j'ai enlev? T_NB parce qu'effet significatif positif
formul = "~ s(T0_50, k=5) + s(S0_50, k=5) + s(lDEPTH, k=4) + s(Month, k=8, by=REGION)"
dat7<-pmza[,c("Cgla_CIV_CVI_occ", "T0_50", "S0_50", "lDEPTH", "Month", "REGION", "fYear")]
mumin_fct_gamm4(RVar="Cgla_CIV_CVI_occ",family="binomial",Variables=formul,nb_simu=10,Data=dat3,nvar=ncol(dat3),Sp="", step="klimit")

formul = "~ s(T_NB, k=5) + s(T0_50, k=5) + s(S0_50, k=5) + s(lDEPTH, k=4) + s(Month, k=6, by=REGION)"
dat5<-pmza[,c("Chyp_CI_CIII_occ", "T_NB", "T0_50", "S0_50", "lDEPTH", "Month", "REGION", "fYear")]
mumin_fct_gamm4(RVar="Chyp_CI_CIII_occ",family="binomial",Variables=formul,nb_simu=10,Data=dat5,nvar=ncol(dat5),Sp="", step="klimit")

formul = "~ s(T_NB, k=5) + s(T0_50, k=5) + s(S0_50, k=5) + s(lDEPTH, k=4) + s(Month, k=8, by=REGION)"
dat9<-pmza[,c("Chyp_CIV_CVI_occ", "T_NB", "T0_50", "S0_50", "lDEPTH", "Month", "REGION", "fYear")]
mumin_fct_gamm4(RVar="Chyp_CIV_CVI_occ",family="binomial",Variables=formul,nb_simu=10,Data=dat9,nvar=ncol(dat9),Sp="", step="klimit")

formul = "~ s(T_NB, k=5) + s(T0_50, k=5) + s(S0_50, k=5) + s(lDEPTH, k=4) + s(Month, k=6, by=REGION)"
dat10<-pmza_chypc4_pres[,c("lChyp_CIV_CVI", "T_NB", "T0_50", "S0_50", "lDEPTH", "Month", "REGION", "fYear")]
mumin_fct_gamm4(RVar="lChyp_CIV_CVI",family="gaussian",Variables=formul,nb_simu=10,Data=dat10,nvar=ncol(dat10),Sp="", step="klimit")

formul = "~ s(T_NB, k=5) + s(T0_50, k=5) + s(S0_50, k=5) + s(lDEPTH, k=6) + s(Month, k=4, by=REGION, bs='cc')"
dat2<-pmza[,c("lCfin_CI_CIV", "T_NB", "T0_50", "S0_50", "lDEPTH", "Month", "REGION", "fYear","STATION")]
#remove outliers
dat2<-subset(dat2, round(dat2$lCfin_CI_CIV,4)!=13.1691 &round(dat2$lCfin_CI_CIV,4)!=13.0756 &round(dat2$lCfin_CI_CIV,4)!=9.9259 & dat2$T_NB <10.2)
mumin_fct_gamm4(RVar="lCfin_CI_CIV",family="gaussian",Variables=formul,nb_simu=10,Data=dat2,nvar=ncol(dat2),Sp="", step="klimit")

formul = "~  s(T0_50, k=5) + s(S0_50, k=5) + s(lDEPTH, k=5) + s(Month, k=5, by=REGION)"
dat2b<-pmza[,c("lCfin_CV_CVI", "T_NB", "T0_50", "S0_50", "lDEPTH", "Month", "REGION", "fYear")]
mumin_fct_gamm4(RVar="lCfin_CV_CVI",family="gaussian",Variables=formul,nb_simu=10,Data=dat2b,nvar=ncol(dat2b),Sp="", step="klimit")


##########################EXPORT########################
library(mgcv)
lfile<-list.files("E:/Projects/ACCASP_FINAL/NEW_ACCASP_pheno/PMZA_Models", pattern="Resmod")
lfin<-lfile[-grep(lfile, pattern=".tiff")]

setwd("E:/Projects/ACCASP_FINAL/NEW_ACCASP_pheno/PMZA_Models")
source("E:/Projects/Scripts_R/TSS_mclapply_ie.R")


load(lfin[1])
load(lfin[2])
load(lfin[3])
load(lfin[4])
load(lfin[5])
load(lfin[6])
load(lfin[7])
load(lfin[8])
load(lfin[9])
load(lfin[10])
load(lfin[11])
load(lfin[12])
load(lfin[13])
load(lfin[14])
########################CFIN############################

#######################CI-CIV###########################
summary(Resmod_lCfin_CI_CIVklimit$gam)
gam.check(Resmod_lCfin_CI_CIVklimit$gam)
dat2<-pmza[,c("lCfin_CI_CIV", "T_NB", "T0_50", "S0_50", "lDEPTH", "Month", "REGION", "fYear","STATION")]
#remove outliers
dat2<-subset(dat2, round(dat2$lCfin_CI_CIV,4)!=13.1691 &round(dat2$lCfin_CI_CIV,4)!=13.0756 &round(dat2$lCfin_CI_CIV,4)!=9.9259 & dat2$T_NB <10.2)
Resmod_lCfin_CI_CIVklimitTSS<-cbind(gam.table(mod=Resmod_lCfin_CI_CIVklimit$gam, Region="Cfin", Stock="CI_CIV"),TSS_mclapply_ie(mcomp =Resmod_lCfin_CI_CIVklimit$gam ,Data=dat2, Sp="lCfin_CI_CIV", gam_type ="mixte" ,fam="gaussian",nrem=10))
write.table(Resmod_lCfin_CI_CIVklimitTSS, file="PMZA_Models_pheno_Results.txt", append=F, sep="\t", dec=".", row.names=F, col.names=T)#col.names=TRUE seulement pour le premier; append=F seulement pour le premier

lab<-"Abundance"
png("GAMM_fig/Resmod_lCfin_CI_CIVklimit.png",height=7,width=8, units="in", res=600)
par(mfrow=c(3,3), cex.lab=1.6,cex.axis=1.3,mar=c(4,5,1,0), yaxs="r")
plot(Resmod_lCfin_CI_CIVklimit$gam, select=2, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE, xlim=c(4,12))
abline(h=0, lty=2, col="azure4")
title("NL", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_lCfin_CI_CIVklimit$gam, select=3, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
title("GSL", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_lCfin_CI_CIVklimit$gam, select=1, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
title("SS", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_lCfin_CI_CIVklimit$gam, select=4, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="S0_50",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_lCfin_CI_CIVklimit$gam, select=5, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="T0_50",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_lCfin_CI_CIVklimit$gam, select=6, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE,xlab="lDEPTH", ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
#plot(Resmod_lCfin_CI_CIVklimit$gam, select=7, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="T_NB",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE)
#abline(h=0, lty=2, col="azure4")
dev.off()  



#######################CV-CVI###########################
summary(Resmod_lCfin_CV_CVIklimit$gam)
gam.check(Resmod_lCfin_CV_CVIklimit$gam)
Resmod_lCfin_CV_CVIklimitTSS<-cbind(gam.table(mod=Resmod_lCfin_CV_CVIklimit$gam, Region="Cfin", Stock="CV_CVI"),TSS_mclapply_ie(mcomp =Resmod_lCfin_CV_CVIklimit$gam ,Data=pmza, Sp="lCfin_CV_CVI", gam_type ="mixte" ,fam="gaussian",nrem=10))
write.table(Resmod_lCfin_CV_CVIklimitTSS, file="PMZA_Models_pheno_Results.txt", append=T, sep="\t", dec=".", row.names=F, col.names=F)#col.names=TRUE seulement pour le premier; append=F seulement pour le premier

lab<-"Abundance"
png("GAMM_fig/Resmod_lCfin_CV_CVIklimit.png",height=7,width=8, units="in", res=600)
par(mfrow=c(3,3), cex.lab=1.6,cex.axis=1.3,mar=c(4,5,1,0), yaxs="r")
plot(Resmod_lCfin_CV_CVIklimit$gam, select=3, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE, xlim=c(4,12))
abline(h=0, lty=2, col="azure4")
title("NL", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_lCfin_CV_CVIklimit$gam, select=4, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
title("GSL", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_lCfin_CV_CVIklimit$gam, select=2, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
title("SS", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_lCfin_CV_CVIklimit$gam, select=1, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="lDEPTH",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_lCfin_CV_CVIklimit$gam, select=5, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="S0_50",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_lCfin_CV_CVIklimit$gam, select=6, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE,xlab="T0_50", ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
dev.off()  



#######################CVI_occ###########################
summary(Resmod_Cfin_CVI_occ$gam)
gam.check(Resmod_Cfin_CVI_occ$gam)
Resmod_Cfin_CVI_occTSS<-cbind(gam.table(mod=Resmod_Cfin_CVI_occ$gam, Region="Cfin", Stock="CVI_occ"),TSS_mclapply_ie(mcomp =Resmod_Cfin_CVI_occ$gam ,Data=pmza, Sp="Cfin_CVI_occ", gam_type ="mixte" ,fam="binomial",nrem=10))
write.table(Resmod_Cfin_CVI_occTSS, file="PMZA_Models_pheno_Results.txt", append=T, sep="\t", dec=".", row.names=F, col.names=F)#col.names=TRUE seulement pour le premier; append=F seulement pour le premier

lab<-"Abundance"
png("GAMM_fig/Resmod_Cfin_CVI_occ.png",height=7,width=8, units="in", res=600)
par(mfrow=c(3,3), cex.lab=1.6,cex.axis=1.3,mar=c(4,5,1,0), yaxs="r")
plot(Resmod_Cfin_CVI_occ$gam, select=3, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE, xlim=c(4,12))
abline(h=0, lty=2, col="azure4")
title("NL", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_Cfin_CVI_occ$gam, select=4, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
title("GSL", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_Cfin_CVI_occ$gam, select=2, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
title("SS", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_Cfin_CVI_occ$gam, select=1, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="lDEPTH",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_Cfin_CVI_occ$gam, select=5, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="S0_50",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_Cfin_CVI_occ$gam, select=6, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE,xlab="T_NB", ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_Cfin_CVI_occ$gam, select=7, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="T0_50",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")

dev.off()  








########################CGLAC############################

#######################CI-CIII###########################
#####################Occurrence########################
summary(Resmod_Cgla_CI_CIII_occklimit$gam)
gam.check(Resmod_Cgla_CI_CIII_occklimit$gam)
Resmod_Cgla_CI_CIII_occklimitTSS<-cbind(gam.table(mod=Resmod_Cgla_CI_CIII_occklimit$gam, Region="Cgla", Stock="CI_CIII"),TSS_mclapply_ie(mcomp =Resmod_Cgla_CI_CIII_occklimit$gam ,Data=pmza, Sp="Cgla_CI_CIII_occ", gam_type ="mixte" ,fam="binomial",nrem=10))
write.table(Resmod_Cgla_CI_CIII_occklimitTSS, file="PMZA_Models_pheno_Results.txt", append=T, sep="\t", dec=".", row.names=F, col.names=F)#col.names=TRUE seulement pour le premier; append=F seulement pour le premier

lab<-"Occurrence"
png("GAMM_fig/Resmod_Cgla_CI_CIII_occklimit.png",height=7,width=8, units="in", res=600)
par(mfrow=c(3,3), cex.lab=1.6,cex.axis=1.3,mar=c(4,5,1,0), yaxs="r")
plot(Resmod_Cgla_CI_CIII_occklimit$gam, select=3, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE, xlim=c(4,12))
abline(h=0, lty=2, col="azure4")
title("NL", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_Cgla_CI_CIII_occklimit$gam, select=4, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
title("GSL", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_Cgla_CI_CIII_occklimit$gam, select=2, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
title("SS", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_Cgla_CI_CIII_occklimit$gam, select=1, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="lDEPTH",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_Cgla_CI_CIII_occklimit$gam, select=5, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="S0_50",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_Cgla_CI_CIII_occklimit$gam, select=6, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE,xlab="T_NB", ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_Cgla_CI_CIII_occklimit$gam, select=7, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="T0_50",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
dev.off()  


#####################Abundance########################
summary(Resmod_lCgla_CI_CIII$gam)
gam.check(Resmod_lCgla_CI_CIII$gam)
Resmod_lCgla_CI_CIIITSS<-cbind(gam.table(mod=Resmod_lCgla_CI_CIII$gam, Region="Cgla", Stock="CI_CIII"),TSS_mclapply_ie(mcomp =Resmod_lCgla_CI_CIII$gam ,Data=pmza, Sp="lCgla_CI_CIII", gam_type ="mixte" ,fam="gaussian",nrem=10))
write.table(Resmod_lCgla_CI_CIIITSS, file="PMZA_Models_pheno_Results.txt", append=T, sep="\t", dec=".", row.names=F, col.names=F)#col.names=TRUE seulement pour le premier; append=F seulement pour le premier

lab<-"Abundance"
png("GAMM_fig/Resmod_lCgla_CI_CIII.png",height=7,width=8, units="in", res=600)
par(mfrow=c(3,3), cex.lab=1.6,cex.axis=1.3,mar=c(4,5,1,0), yaxs="r")
plot(Resmod_lCgla_CI_CIII$gam, select=3, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE, xlim=c(4,12))
abline(h=0, lty=2, col="azure4")
title("NL", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_lCgla_CI_CIII$gam, select=4, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
title("GSL", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_lCgla_CI_CIII$gam, select=2, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
title("SS", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_lCgla_CI_CIII$gam, select=1, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="lDEPTH",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_lCgla_CI_CIII$gam, select=5, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="S0_50",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_lCgla_CI_CIII$gam, select=6, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE,xlab="T_NB", ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_lCgla_CI_CIII$gam, select=7, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="T0_50",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
dev.off()  



#######################CIV-CVI###########################
#####################Occurrence########################
summary(Resmod_Cgla_CIV_CVI_occklimit$gam)
gam.check(Resmod_Cgla_CIV_CVI_occklimit$gam)
Resmod_Cgla_CIV_CVI_occklimitTSS<-cbind(gam.table(mod=Resmod_Cgla_CIV_CVI_occklimit$gam, Region="Cgla", Stock="CIV_CVI"),TSS_mclapply_ie(mcomp =Resmod_Cgla_CIV_CVI_occklimit$gam ,Data=pmza, Sp="Cgla_CIV_CVI_occ", gam_type ="mixte" ,fam="binomial",nrem=10))
write.table(Resmod_Cgla_CIV_CVI_occklimitTSS, file="PMZA_Models_pheno_Results.txt", append=T, sep="\t", dec=".", row.names=F, col.names=F)#col.names=TRUE seulement pour le premier; append=F seulement pour le premier

lab<-"Occurrence"
png("GAMM_fig/Resmod_Cgla_CIV_CVI_occklimit.png",height=7,width=8, units="in", res=600)
par(mfrow=c(3,3), cex.lab=1.6,cex.axis=1.3,mar=c(4,5,1,0), yaxs="r")
plot(Resmod_Cgla_CIV_CVI_occklimit$gam, select=3, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE, xlim=c(4,12))
abline(h=0, lty=2, col="azure4")
title("NL", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_Cgla_CIV_CVI_occklimit$gam, select=4, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
title("GSL", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_Cgla_CIV_CVI_occklimit$gam, select=2, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
title("SS", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_Cgla_CIV_CVI_occklimit$gam, select=1, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="lDEPTH",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_Cgla_CIV_CVI_occklimit$gam, select=5, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="S0_50",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_Cgla_CIV_CVI_occklimit$gam, select=6, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE,xlab="T0_50", ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
#plot(Resmod_Cgla_CIV_CVI_occklimit$gam, select=7, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="T_NB",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE)
#abline(h=0, lty=2, col="azure4")
dev.off()  


#####################Abundance########################
summary(Resmod_lCgla_CIV_CVI$gam)
gam.check(Resmod_lCgla_CIV_CVI$gam)

dat8<-pmza_cglac4_pres[,c("lCgla_CIV_CVI", "T_NB", "T0_50", "S0_50", "lDEPTH", "Month", "REGION", "fYear")]
dat8<-subset(dat8, dat8$T_NB <10)

Resmod_lCgla_CIV_CVITSS<-cbind(gam.table(mod=Resmod_lCgla_CIV_CVI$gam, Region="Cgla", Stock="CIV_CVI"),TSS_mclapply_ie(mcomp =Resmod_lCgla_CIV_CVI$gam ,Data=dat8, Sp="lCgla_CIV_CVI", gam_type ="mixte" ,fam="gaussian",nrem=10))
write.table(Resmod_lCgla_CIV_CVITSS, file="PMZA_Models_pheno_Results.txt", append=T, sep="\t", dec=".", row.names=F, col.names=F)#col.names=TRUE seulement pour le premier; append=F seulement pour le premier

lab<-"Abundance"
png("GAMM_fig/Resmod_lCgla_CIV_CVI.png",height=7,width=8, units="in", res=600)
par(mfrow=c(3,3), cex.lab=1.6,cex.axis=1.3,mar=c(4,5,1,0), yaxs="r")
plot(Resmod_lCgla_CIV_CVI$gam, select=3, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE, xlim=c(4,12))
abline(h=0, lty=2, col="azure4")
title("NL", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_lCgla_CIV_CVI$gam, select=4, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
title("GSL", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_lCgla_CIV_CVI$gam, select=2, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
title("SS", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_lCgla_CIV_CVI$gam, select=1, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="lDEPTH",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_lCgla_CIV_CVI$gam, select=5, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="S0_50",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_lCgla_CIV_CVI$gam, select=6, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE,xlab="T_NB", ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_lCgla_CIV_CVI$gam, select=7, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="T0_50",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
dev.off() 



########################Chyp############################

#######################CI-CIII###########################
#####################Occurrence########################
summary(Resmod_Chyp_CI_CIII_occklimit$gam)
gam.check(Resmod_Chyp_CI_CIII_occklimit$gam)
Resmod_Chyp_CI_CIII_occklimitTSS<-cbind(gam.table(mod=Resmod_Chyp_CI_CIII_occklimit$gam, Region="Chyp", Stock="CI_CIII"),TSS_mclapply_ie(mcomp =Resmod_Chyp_CI_CIII_occklimit$gam ,Data=pmza, Sp="Chyp_CI_CIII_occ", gam_type ="mixte" ,fam="binomial",nrem=10))
write.table(Resmod_Chyp_CI_CIII_occklimitTSS, file="PMZA_Models_pheno_Results.txt", append=T, sep="\t", dec=".", row.names=F, col.names=F)#col.names=TRUE seulement pour le premier; append=F seulement pour le premier

lab<-"Occurrence"
png("GAMM_fig/Resmod_Chyp_CI_CIII_occklimit.png",height=7,width=8, units="in", res=600)
par(mfrow=c(3,3), cex.lab=1.6,cex.axis=1.3,mar=c(4,5,1,0), yaxs="r")
plot(Resmod_Chyp_CI_CIII_occklimit$gam, select=3, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE, xlim=c(4,12))
abline(h=0, lty=2, col="azure4")
title("NL", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_Chyp_CI_CIII_occklimit$gam, select=4, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
title("GSL", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_Chyp_CI_CIII_occklimit$gam, select=2, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
title("SS", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_Chyp_CI_CIII_occklimit$gam, select=1, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="lDEPTH",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_Chyp_CI_CIII_occklimit$gam, select=5, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="S0_50",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_Chyp_CI_CIII_occklimit$gam, select=6, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE,xlab="T_NB", ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_Chyp_CI_CIII_occklimit$gam, select=7, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="T0_50",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
dev.off()  


#####################Abundance########################
summary(Resmod_lChyp_CI_CIII$gam)
gam.check(Resmod_lChyp_CI_CIII$gam)
Resmod_lChyp_CI_CIIITSS<-cbind(gam.table(mod=Resmod_lChyp_CI_CIII$gam, Region="Chyp", Stock="CI_CIII"),TSS_mclapply_ie(mcomp =Resmod_lChyp_CI_CIII$gam ,Data=pmza, Sp="lChyp_CI_CIII", gam_type ="mixte" ,fam="gaussian",nrem=10))
write.table(Resmod_lChyp_CI_CIIITSS, file="PMZA_Models_pheno_Results.txt", append=T, sep="\t", dec=".", row.names=F, col.names=F)#col.names=TRUE seulement pour le premier; append=F seulement pour le premier

lab<-"Abundance"
png("GAMM_fig/Resmod_lChyp_CI_CIII.png",height=7,width=8, units="in", res=600)
par(mfrow=c(3,3), cex.lab=1.6,cex.axis=1.3,mar=c(4,5,1,0), yaxs="r")
plot(Resmod_lChyp_CI_CIII$gam, select=2, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE, xlim=c(4,12))
abline(h=0, lty=2, col="azure4")
title("NL", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_lChyp_CI_CIII$gam, select=3, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
title("GSL", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_lChyp_CI_CIII$gam, select=1, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
title("SS", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_lChyp_CI_CIII$gam, select=4, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="S0_50",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_lChyp_CI_CIII$gam, select=5, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="T0_50",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
dev.off()  



#######################CIV-CVI###########################
#####################Occurrence########################
summary(Resmod_Chyp_CIV_CVI_occklimit$gam)
gam.check(Resmod_Chyp_CIV_CVI_occklimit$gam)
Resmod_Chyp_CIV_CVI_occklimitTSS<-cbind(gam.table(mod=Resmod_Chyp_CIV_CVI_occklimit$gam, Region="Chyp", Stock="CIV_CVI"),TSS_mclapply_ie(mcomp =Resmod_Chyp_CIV_CVI_occklimit$gam ,Data=pmza, Sp="Chyp_CIV_CVI_occ", gam_type ="mixte" ,fam="binomial",nrem=10))
write.table(Resmod_Chyp_CIV_CVI_occklimitTSS, file="PMZA_Models_pheno_Results.txt", append=T, sep="\t", dec=".", row.names=F, col.names=F)#col.names=TRUE seulement pour le premier; append=F seulement pour le premier

lab<-"Occurrence"
png("GAMM_fig/Resmod_Chyp_CIV_CVI_occklimit.png",height=7,width=8, units="in", res=600)
par(mfrow=c(3,3), cex.lab=1.6,cex.axis=1.3,mar=c(4,5,1,0), yaxs="r")
plot(Resmod_Chyp_CIV_CVI_occklimit$gam, select=3, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE, xlim=c(4,12))
abline(h=0, lty=2, col="azure4")
title("NL", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_Chyp_CIV_CVI_occklimit$gam, select=4, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
title("GSL", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_Chyp_CIV_CVI_occklimit$gam, select=2, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
title("SS", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_Chyp_CIV_CVI_occklimit$gam, select=1, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="lDEPTH",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_Chyp_CIV_CVI_occklimit$gam, select=5, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="S0_50",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_Chyp_CIV_CVI_occklimit$gam, select=6, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE,xlab="T_NB", ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_Chyp_CIV_CVI_occklimit$gam, select=7, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="T0_50",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
dev.off()  


#####################Abundance########################
summary(Resmod_lChyp_CIV_CVIklimit$gam)
gam.check(Resmod_lChyp_CIV_CVIklimit$gam)
Resmod_lChyp_CIV_CVIklimitTSS<-cbind(gam.table(mod=Resmod_lChyp_CIV_CVIklimit$gam, Region="Chyp", Stock="CIV_CVI"),TSS_mclapply_ie(mcomp =Resmod_lChyp_CIV_CVIklimit$gam ,Data=pmza, Sp="lChyp_CIV_CVI", gam_type ="mixte" ,fam="gaussian",nrem=10))
write.table(Resmod_lChyp_CIV_CVIklimitTSS, file="PMZA_Models_pheno_Results.txt", append=T, sep="\t", dec=".", row.names=F, col.names=F)#col.names=TRUE seulement pour le premier; append=F seulement pour le premier

lab<-"Abundance"
png("GAMM_fig/Resmod_lChyp_CIV_CVIklimit.png",height=7,width=8, units="in", res=600)
par(mfrow=c(3,3), cex.lab=1.6,cex.axis=1.3,mar=c(4,5,1,0), yaxs="r")
plot(Resmod_lChyp_CIV_CVIklimit$gam, select=3, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE, xlim=c(4,12))
abline(h=0, lty=2, col="azure4")
title("NL", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_lChyp_CIV_CVIklimit$gam, select=4, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
title("GSL", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_lChyp_CIV_CVIklimit$gam, select=2, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
title("SS", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_lChyp_CIV_CVIklimit$gam, select=1, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="lDEPTH",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_lChyp_CIV_CVIklimit$gam, select=5, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="S0_50",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_lChyp_CIV_CVIklimit$gam, select=6, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE,xlab="T_NB", ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_lChyp_CIV_CVIklimit$gam, select=7, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="T0_50",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
dev.off() 




##############Paraalanus################
#######################CIV-CVI###########################
#####################Occurrence########################
summary(Resmod_Para_occ$gam)
gam.check(Resmod_Para_occ$gam)
Resmod_Para_occTSS<-cbind(gam.table(mod=Resmod_Para_occ$gam, Region="Para", Stock="CI_CVI"),TSS_mclapply_ie(mcomp =Resmod_Para_occ$gam ,Data=pmza, Sp="Para_occ", gam_type ="mixte" ,fam="binomial",nrem=10))
write.table(Resmod_Para_occTSS, file="PMZA_Models_pheno_Results.txt", append=T, sep="\t", dec=".", row.names=F, col.names=F)#col.names=TRUE seulement pour le premier; append=F seulement pour le premier

lab<-"Occurrence"
png("GAMM_fig/Resmod_Para_occ.png",height=7,width=8, units="in", res=600)
par(mfrow=c(3,3), cex.lab=1.6,cex.axis=1.3,mar=c(4,5,1,0), yaxs="r")
plot(Resmod_Para_occ$gam, select=3, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE, xlim=c(4,12))
abline(h=0, lty=2, col="azure4")
title("NL", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_Para_occ$gam, select=4, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
title("GSL", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_Para_occ$gam, select=2, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
title("SS", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_Para_occ$gam, select=1, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="lDEPTH",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_Para_occ$gam, select=5, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="S0_50",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_Para_occ$gam, select=6, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE,xlab="T_NB", ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_Para_occ$gam, select=7, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="T0_50",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
dev.off()  


#####################Abundance########################
summary(Resmod_lPara$gam)
gam.check(Resmod_lPara$gam)
Resmod_lParaTSS<-cbind(gam.table(mod=Resmod_lPara$gam, Region="Para", Stock="CI_CVI"),TSS_mclapply_ie(mcomp =Resmod_lPara$gam ,Data=pmza, Sp="lPara", gam_type ="mixte" ,fam="gaussian",nrem=10))
write.table(Resmod_lParaTSS, file="PMZA_Models_pheno_Results.txt", append=T, sep="\t", dec=".", row.names=F, col.names=F)#col.names=TRUE seulement pour le premier; append=F seulement pour le premier

lab<-"Abundance"
png("GAMM_fig/Resmod_lPara.png",height=7,width=8, units="in", res=600)
par(mfrow=c(3,3), cex.lab=1.6,cex.axis=1.3,mar=c(4,5,1,0), yaxs="r")
plot(Resmod_lPara$gam, select=2, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE, xlim=c(4,12))
abline(h=0, lty=2, col="azure4")
title("NL", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_lPara$gam, select=3, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
title("GSL", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_lPara$gam, select=1, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="Month",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
title("SS", line = -1.4, adj=0.05, cex.main=1.4)
plot(Resmod_lPara$gam, select=4, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="S0_50",ylab=paste0("Effect on ",lab), pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_lPara$gam, select=5, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE, xlab="T_NB",ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
plot(Resmod_lPara$gam, select=6, scale=0,residuals=FALSE, shade=TRUE, seWithMean=TRUE,xlab="T0_50", ylab="", pch=1, cex=1, all.terms=TRUE)
abline(h=0, lty=2, col="azure4")
dev.off() 


                
        



