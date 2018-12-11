#les donn.es manquantes sont en partie enlev?es
library(doBy)
pmza<-read.csv("ACCASP_all_stations_Env-Copdata.txt")

pmza<-pmza[,c(1,2,3,4,6,7,8,9,10,12,13,17,67,32:37,39:50,52:59)]
pmza<-renameCol(pmza, c("YEAR","MONTH", "DAY"), c("Year","Month", "Day"))


pmza$lSTRAT <- log1p(pmza$STRAT_corrNL)
pmza$lDEPTH<-log1p(pmza$DEPTH)
pmza$fYear<-as.factor(pmza$Year)

pmza<-na.omit(pmza)

#generate a unique ID to select data
pmza$ID<-seq(1,nrow(pmza),1)

pmza$test_occ <- 1


library(doBy)
library(plyr)
#data duplicate
################Remove duplicate################
tabl_dupl<-summaryBy(REGION~ REGION  + Year+Month + Day +LATITUDE+LONGITUDE + STATION + DEPTH, FUN=length, data=pmza)
#which metadata have more than one entry
duplicate<-subset(tabl_dupl,tabl_dupl$REGION.length!=1)

#select data that have more than one entry
aberr<-join(duplicate, pmza, type="inner")#subset avec pmza$STATION %in% duplicate$STATION... ne fonctionne pas car trop de combinaisons possibles

#verify data
oaberr<-orderBy(~ REGION+Year+Month+Day+STATION+LATITUDE, data=aberr)

#add sum of data. If one entry >0 and one ==0
pmza_sum<-summaryBy(Cfin_CI + Cfin_CII + Cfin_CIII +  Cfin_CIV + Cfin_CV + Cfin_CVI + Cgla_CI + Cgla_CII + Cgla_CIII + Cgla_CIV + Cgla_CV + Cgla_CVI + Chyp_CI + Chyp_CII + Chyp_CIII + Chyp_CIV + Chyp_CV + Chyp_CVI + ParaC_CI + ParaC_CII + ParaC_CIII + ParaC_CIV + ParaC_CV + ParaC_CVI + ParaC_CI_CV + ParaC_none + test_occ ~ REGION+Year+Month+Day+SEASON+STATION+TRANSECT+LATITUDE+LONGITUDE, data=aberr, FUN = sum, na.rm=T, keep.names=T)


#Sum =2 do not use sum, use mean. means that its the same entry on all variables.
occ2<-subset(pmza_sum, pmza_sum$test_occ>1)
occ1<-subset(pmza_sum, pmza_sum$test_occ<2)#corrected OK

#occ2 is not OK use mean
occ2<-occ2[,c(1:8)]
aberr_occ2<-join(pmza,occ2, type="inner", by=c("REGION", "Year", "Month", "Day", "STATION", "TRANSECT", "LATITUDE", "LONGITUDE"))#subset avec pmza$STATION %in% duplicate$STATION... ne fonctionne pas car trop de combinaisons possibles
#verify data
oaberr2<-orderBy(~ REGION+Year+Month+Day+STATION+LATITUDE, data=aberr_occ2)
aberr_occ2_mean<-summaryBy(Cfin_CI + Cfin_CII + Cfin_CIII +  Cfin_CIV + Cfin_CV + Cfin_CVI + Cgla_CI + Cgla_CII + Cgla_CIII + Cgla_CIV + Cgla_CV + Cgla_CVI + Chyp_CI + Chyp_CII + Chyp_CIII + Chyp_CIV + Chyp_CV + Chyp_CVI + ParaC_CI + ParaC_CII + ParaC_CIII + ParaC_CIV + ParaC_CV + ParaC_CVI + ParaC_CI_CV + ParaC_none + test_occ~ REGION+Year+Month+Day+SEASON+STATION+TRANSECT+LATITUDE+LONGITUDE, data=aberr_occ2, FUN = mean, na.rm=T)
colnames(aberr_occ2_mean)<- colnames(occ1)

#4entry 2x 0 and 2x 1 =0.5
aberr_occ2_05<-subset(aberr_occ2_mean, aberr_occ2_mean$test_occ==0.5)
aberr_occ2_1<-subset(aberr_occ2_mean, aberr_occ2_mean$test_occ!=0.5)

occ2<-join(occ2,aberr_occ2_1, type="inner")

pmza_sum2<-rbind(occ1, occ2)


#join data with metadata
pmza_aberr1<-join(pmza_sum2, aberr[,c(1:8,10:14,41,42,43,44)], type="left",match="first")
#which are not duplicated
pmza_noaberr<-subset(pmza, !pmza$ID%in% aberr$ID)
#subset columns
pmza_noaberr2<-pmza_noaberr[,colnames(pmza_aberr1)]
#join the duplicates and non duplicates
pmza<-as.data.frame(rbind(pmza_aberr1, pmza_noaberr2))


#data duplicate
tabl_dupl<-summaryBy(REGION~ REGION  + Year+Month + Day +LATITUDE+LONGITUDE + STATION + DEPTH, FUN=length, data=pmza)
#which metadata have more than one entry
duplicate<-subset(tabl_dupl,tabl_dupl$REGION.length!=1)
nrow(duplicate)
###no more duplicate


pmza$Date<-as.Date(paste(pmza$Year, pmza$Month, pmza$Day, sep="-"), format="%Y-%m-%d")
pmza$doy<-as.numeric(format(pmza$Date, format="%j"))


#je ne souvient pas pourquoi pour aberr2. je pense que c'?tait un r?sidus tr?s loin des autres
aberr2<-subset(pmza,pmza$Year==2011 & pmza$Month==7  & pmza$T_NB <5 & pmza$T0_50<10 & pmza$T0_50 >9.5 & pmza$lDEPTH <3)

pmza<-subset(pmza, pmza$ID!=aberr2[1,"ID"])

pmza<-subset(pmza, select= -test_occ)


####groups######

pmza$Cfin_CI_CIV <- apply(pmza[,c("Cfin_CI","Cfin_CII","Cfin_CIII","Cfin_CIV")], 1, FUN=sum, na.rm=T)
pmza$Cfin_CV_CVI <- apply(pmza[,c("Cfin_CV","Cfin_CVI")], 1, FUN=sum, na.rm=T)

pmza$Cgla_CI_CIII <- apply(pmza[,c("Cgla_CI","Cgla_CII","Cgla_CIII")], 1, FUN=sum, na.rm=T)
pmza$Cgla_CIV_CVI <- apply(pmza[,c("Cgla_CIV","Cgla_CV","Cgla_CVI")], 1, FUN=sum, na.rm=T)

pmza$Chyp_CI_CIII <- apply(pmza[,c("Chyp_CI","Chyp_CII","Chyp_CIII")], 1, FUN=sum, na.rm=T)
pmza$Chyp_CIV_CVI <- apply(pmza[,c("Chyp_CIV","Chyp_CV","Chyp_CVI")], 1, FUN=sum, na.rm=T)

pmza$Para <- apply(pmza[,c("ParaC_CI","ParaC_CII","ParaC_CIII","ParaC_CIV","ParaC_CV","ParaC_CVI","ParaC_CI_CV","ParaC_none")], 1, FUN=sum, na.rm=T)


pmza$Cfin_CI_CIV_occ <- ifelse(pmza$Cfin_CI_CIV>0,1,0)
pmza$Cfin_CV_CVI_occ <- ifelse(pmza$Cfin_CV_CVI>0,1,0)

pmza$Chyp_CI_CIII_occ <- ifelse(pmza$Chyp_CI_CIII>0, 1,0)
pmza$Chyp_CIV_CVI_occ <- ifelse(pmza$Chyp_CIV_CVI>0, 1,0)

pmza$Cgla_CI_CIII_occ <- ifelse(pmza$Cgla_CI_CIII>0, 1,0)
pmza$Cgla_CIV_CVI_occ <- ifelse(pmza$Cgla_CIV_CVI>0, 1,0)

pmza$Para_occ <- ifelse(pmza$Para>0, 1,0)
pmza$Para_occ <- ifelse(pmza$Para>0, 1,0)


pmza$lCfin_CI_CIV<-log1p(pmza$Cfin_CI_CIV)
pmza$lCfin_CV_CVI<-log1p(pmza$Cfin_CV_CVI)

pmza$lChyp_CI_CIII<-log1p(pmza$Chyp_CI_CIII)
pmza$lChyp_CIV_CVI<-log1p(pmza$Chyp_CIV_CVI)

pmza$lCgla_CI_CIII<-log1p(pmza$Cgla_CI_CIII)
pmza$lCgla_CIV_CVI<-log1p(pmza$Cgla_CIV_CVI)

pmza$lPara <- log1p(pmza$Para)


pmza_cfinc5_pres<-subset(pmza, pmza$Cfin_CV_CVI_occ ==1)
pmza_chypc4_pres<-subset(pmza, pmza$Chyp_CIV_CVI_occ ==1)
pmza_cglac4_pres<-subset(pmza, pmza$Cgla_CIV_CVI_occ ==1)

pmza_cfinc1_pres<-subset(pmza, pmza$Cfin_CI_CIV_occ ==1)
pmza_chypc1_pres<-subset(pmza, pmza$Chyp_CI_CIII_occ ==1)
pmza_cglac1_pres<-subset(pmza, pmza$Cgla_CI_CIII_occ ==1)

pmza_para_pres<-subset(pmza, pmza$Para_occ==1)

rm(aberr,aberr_occ2,aberr_occ2_05,aberr_occ2_1,aberr_occ2_mean,aberr2,duplicate,oaberr,oaberr2,occ1,occ2,pmza_aberr1,pmza_noaberr,pmza_noaberr2,pmza_sum,pmza_sum2,tabl_dupl)

#write.table(pmza, "E:/Projects/ACCASP_FINAL/NEW_ACCASP_pheno/Data/ACCASP_all_stations_Env-Copdata_complete_without_duplicate.txt", row.names=F,sep="\t", dec=".")
       