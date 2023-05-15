######################################
###### prepare dataframe #############
######################################
{BaselineCohort <- rbind(NApatientAgree[,c(1,2,3,9,11,13,15)],PairedWithInTimeSelected[,c(1,2,3,9,11,13,12)])
BaselineCohort <- merge(BaselineCohort,final_check[,c("PROCHI","LVEF")],by="PROCHI", all.x=TRUE)
BaselineCohort$FinalDiagnosis <- as.character(BaselineCohort$FinalDiagnosis)
for (i in 1:dim(final_check)[1]) {
  p <- as.character(final_check[i,"PROCHI"]) 
  BaselineCohort[BaselineCohort$PROCHI==p,"FinalDiagnosis"] <- as.character(final_check[i,"ClinicalDiagnosis"]) 
}
BaselineCohort$FinalDiagnosis <- as.factor(BaselineCohort$FinalDiagnosis)
colnames(BaselineCohort)[8] <- "LVEFDiagnosis"
BaselineCohort[BaselineCohort$FinalDiagnosis=="HFmrEF","FinalDiagnosis"] <- "HFrEF"

BaselineCohort<- merge(Demography[,c(1,3,5,7)], BaselineCohort, by = c("PROCHI"))


}
BaselineCohort_08 <- BaselineCohort
########################################################
###### prepare dataframe for clinical echo #############
########################################################

# combine clinical EKO
ESC2014 <- Echo_TaysideFrom2014[,c(1,5)]
colnames(ESC2014) <- c("PROCHI","Date")
ESC2014$LVEF <- round((Echo_TaysideFrom2014[,c(31)]+Echo_TaysideFrom2014[,c(34)])/2,2)
ESC2014$patientHeight <- Echo_TaysideFrom2014[,c(160)]
ESC2014$patientWeight <- Echo_TaysideFrom2014[,c(161)]
ESC2014$LVMassIndex <- Echo_TaysideFrom2014[,c(163)]
ESC2014$PASystolicPressure <- NA
ESC2014$TRVelocity <- Echo_TaysideFrom2014[,c(45)]
ESC2014$RWT <- round(2*Echo_TaysideFrom2014[,c(71)]/Echo_TaysideFrom2014[,c(26)],2)

ESC2013 <- Echo_TaysideTo2013[,c(1,9,49)]
colnames(ESC2013) <- c("PROCHI","Date","LVEF")
ESC2013$patientHeight <- NA
ESC2013$patientWeight <- NA
ESC2013$LVMassIndex <- NA
ESC2013$PASystolicPressure <- Echo_TaysideTo2013[,c(42)]
ESC2013$TRVelocity <- Echo_TaysideTo2013[,c(25)]
ESC2013$RWT <- round(2*Echo_TaysideTo2013[,c(56)]/Echo_TaysideTo2013[,c(46)],2)

ESC <- rbind(ESC2014,ESC2013)

# add clinical EKO
ESC_Selected <- merge(ESC, BaselineCohort_08[,c("PROCHI","NextEchoDate")], by =c("PROCHI") )
ESC_Selected$Time_Diff <- as.Date(as.character(ESC_Selected$Date),format="%d/%m/%Y")-as.Date(as.character(ESC_Selected$NextEchoDate),format="%Y-%m-%d")
ESC_Selected <- ESC_Selected[ESC_Selected$Time_Diff==0,]
BaselineCohort_08 <- merge(BaselineCohort_08, ESC_Selected[,c(1,4:(dim(ESC_Selected)[2]-2))], by = c("PROCHI"), all.x = TRUE)
BaselineCohort_08 <- unique(BaselineCohort_08)

BaselineCohort_08$duplicate <- duplicated(BaselineCohort_08$PROCHI)
t <- BaselineCohort_08[BaselineCohort_08$duplicate=='TRUE',]
t2 <- BaselineCohort_08[BaselineCohort_08$PROCHI %in% t$PROCHI,]
BaselineCohort_08 <- BaselineCohort_08[-c(20,220,416,472,538),]

# combine US2AI
# read the EKO AI image analysis results update
setwd("P:/Project 3634/Dundee (2)")
getwd()
EKOImageData_update <- read.csv("MeasurementsDundeeJan3.csv")

EKOImageData_update_selected <- EKOImageData_update[,c(1,5,6,10,11,12)]
colnames(EKOImageData_update_selected) <- c("PROCHI","US2AILAVolumn", "US2AILVMass", "US2AIEeRatio", "US2AITRVelecity","US2AIRWT")
BaselineCohort_08 <- merge(BaselineCohort_08, EKOImageData_update_selected, by = "PROCHI", all.x = TRUE)

# calculate BSA
BaselineCohort_08$BSA <- 0.007184*(BaselineCohort_08$patientHeight^0.725)*(BaselineCohort_08$patientWeight^0.425)
# calculate BSA
BaselineCohort_08$BMI <- BaselineCohort_08$patientWeight/((BaselineCohort_08$patientHeight/100)^2)
# index
BaselineCohort_08$US2AILAVolumnIndex <- BaselineCohort_08$US2AILAVolumn/BaselineCohort_08$BSA
# change inf into na
BaselineCohort_08[!is.na(BaselineCohort_08$US2AILAVolumnIndex) & BaselineCohort_08$US2AILAVolumnIndex==Inf,"US2AILAVolumnIndex"] <- NA
# index
BaselineCohort_08$US2AILVMassIndex <- BaselineCohort_08$US2AILVMass/BaselineCohort_08$BSA
# change inf into na
BaselineCohort_08[!is.na(BaselineCohort_08$US2AILVMassIndex) & BaselineCohort_08$US2AILVMassIndex==Inf,"US2AILVMassIndex"] <- NA 

##### compare the fraction empty between Echo report and US2AI ####
{
  BaselineCohort <- BaselineCohort_08[,c("PROCHI","FinalDiagnosis","LVMassIndex","TRVelocity","RWT","US2AIEeRatio","US2AITRVelecity","US2AIRWT","US2AILAVolumnIndex", "US2AILVMassIndex","US2AILAVolumn","US2AILVMass")]
  t <- matrix(nrow=0,ncol = 0)
  compare_Summary <- data.frame(t)
  compare_Summary["LV mass Index","EHR"] <- round((sum(is.na(BaselineCohort$LVMassIndex)==FALSE)/nrow(BaselineCohort))*100,2)
  compare_Summary["LV mass Index","US2AI"] <- round((sum(is.na(BaselineCohort$US2AILVMass )==FALSE)/nrow(BaselineCohort))*100,2)
  compare_Summary["RWT","EHR"] <- round((sum(is.na(BaselineCohort$RWT )==FALSE)/nrow(BaselineCohort))*100,2)
  compare_Summary["RWT","US2AI"] <- round((sum(is.na(BaselineCohort$US2AIRWT )==FALSE)/nrow(BaselineCohort))*100,2)
  compare_Summary["LA volumn index","EHR"] <- round((0/nrow(BaselineCohort))*100,2)
  compare_Summary["LA volumn index","US2AI"] <- round((sum(is.na(BaselineCohort$US2AILAVolumn )==FALSE)/nrow(BaselineCohort))*100,2)
  compare_Summary["E/e","EHR"] <- round((0/nrow(BaselineCohort))*100,2)
  compare_Summary["E/e","US2AI"] <- round((sum(is.na(BaselineCohort$US2AIEeRatio )==FALSE)/nrow(BaselineCohort))*100,2)
  compare_Summary["TRVelocity","EHR"] <- round((sum(is.na(BaselineCohort$TRVelocity )==FALSE)/nrow(BaselineCohort))*100,2)
  compare_Summary["TRVelocity","US2AI"] <- round((sum(is.na(BaselineCohort$US2AITRVelecity )==FALSE)/nrow(BaselineCohort))*100,2)
  
  }


 
 # select columns
 BaselineCohort_08 <- BaselineCohort_08[BaselineCohort_08$FinalDiagnosis=="HFpEF"|BaselineCohort_08$FinalDiagnosis=="Control",
                                        c("PROCHI","FinalDiagnosis","sex",
                                          "LVMassIndex","RWT","TRVelocity",
                                          "US2AILVMassIndex","US2AIRWT","US2AILAVolumnIndex","US2AIEeRatio","US2AITRVelecity")]
#12052023
#### generate a combined ESC #########
 
 # combine LVMassIndex
 BaselineCohort_08 <- BaselineCohort_08 %>% mutate(LVMassIndexCombined = case_when(
   LVMassIndex == 0.0|is.na(LVMassIndex) ~ US2AILVMassIndex,
   TRUE ~ LVMassIndex
 ))
 
 # combine RWT
 BaselineCohort_08 <- BaselineCohort_08 %>% mutate(RWTCombined = case_when(
   is.na(RWT) ~ US2AIRWT,
   TRUE ~ RWT
 ))
 
 # combine LAVolumnIndex
 BaselineCohort_08$LAVolumnIndexCombined <- BaselineCohort_08$US2AILAVolumnIndex
 
 # combine EeRatio
 BaselineCohort_08$EeRatioCombined <- BaselineCohort_08$US2AIEeRatio
 
 # combine TRVelocity
 BaselineCohort_08$TRVelocity <- as.numeric(as.character(BaselineCohort_08$TRVelocity))/100
 BaselineCohort_08 <- BaselineCohort_08 %>% mutate(TRVelocityCombined = case_when(
   TRVelocity == 0.0|is.na(TRVelocity) ~ US2AITRVelecity,
   TRUE ~ TRVelocity
 ))
 
 colnames(BaselineCohort_08)
 # "PROCHI"                "FinalDiagnosis"        "sex"                   "LVMassIndex"           "RWT"                   "TRVelocity"           
 # "US2AILVMassIndex"      "US2AIRWT"              "US2AILAVolumnIndex"    "US2AIEeRatio"          "US2AITRVelecity"      
 # "LVMassIndexCombined"   "RWTUpdate"        "LAVolumnIndexCombined" "EeRatioCombined"       "TRVelocityCombined"   
 #first row is EHR ESC
 #second row is US2AI ESC
 #third row is combined ESC
 
######################################
########## ESC analysis ##############
######################################
BaselineCohort_08[BaselineCohort_08==0] <- NA
 BaselineCohort_08[BaselineCohort_08=="NULL"] <- NA
########## Combined########
# LVmass
BaselineCohort_08 <- BaselineCohort_08 %>% mutate(LVmassbiggerthan95 = case_when(
  LVMassIndexCombined >= 95 & sex== 'F' ~ 1,
  LVMassIndexCombined >= 115 & sex== 'M' ~ 1,
  TRUE ~ 0
))

 # RWT
 BaselineCohort_08 <- BaselineCohort_08 %>% mutate(RWTbiggerthan042 = case_when(
   RWTCombined > 0.42 ~ 1,
   TRUE ~ 0
 ))
 
 # LAVolumnIndex
 BaselineCohort_08 <- BaselineCohort_08 %>% mutate(LAVolumnIndexbiggerthan34 = case_when(
   LAVolumnIndexCombined > 34 ~ 1,
   TRUE ~ 0
 ))
 
 # EeRatio
 BaselineCohort_08 <- BaselineCohort_08 %>% mutate(EeRatiobiggerthan9 = case_when(
   EeRatioCombined > 9 ~ 1,
   TRUE ~ 0
 ))
 
# TRVelocity
BaselineCohort_08 <- BaselineCohort_08 %>% mutate(TRVelocitybiggerthan28 = case_when(
 TRVelocityCombined > 2.8 ~ 1,
 TRUE ~ 0
))
 
BaselineCohort_08$Combined <- BaselineCohort_08$LVmassbiggerthan95 +
  BaselineCohort_08$RWTbiggerthan042 +
  BaselineCohort_08$LAVolumnIndexbiggerthan34 +
  BaselineCohort_08$EeRatiobiggerthan9 +
  BaselineCohort_08$TRVelocitybiggerthan28

########## EHR ########
# LVmass
BaselineCohort_08 <- BaselineCohort_08 %>% mutate(LVmassbiggerthan95 = case_when(
  LVMassIndex >= 95 & sex== 'F' ~ 1,
  LVMassIndex >= 115 & sex== 'M' ~ 1,
  TRUE ~ 0
))

# RWT
BaselineCohort_08 <- BaselineCohort_08 %>% mutate(RWTbiggerthan042 = case_when(
  RWT > 0.42 ~ 1,
  TRUE ~ 0
))

# TRVelocity
BaselineCohort_08 <- BaselineCohort_08 %>% mutate(TRVelocitybiggerthan28 = case_when(
  TRVelocity > 2.8 ~ 1,
  TRUE ~ 0
))

BaselineCohort_08$EHR <- BaselineCohort_08$LVmassbiggerthan95 +
  BaselineCohort_08$RWTbiggerthan042 +
  BaselineCohort_08$TRVelocitybiggerthan28


########## US2AI ########
# LVmass
BaselineCohort_08 <- BaselineCohort_08 %>% mutate(LVmassbiggerthan95 = case_when(
  US2AILVMassIndex >= 95 & sex== 'F' ~ 1,
  US2AILVMassIndex >= 115 & sex== 'M' ~ 1,
  TRUE ~ 0
))

# RWT
BaselineCohort_08 <- BaselineCohort_08 %>% mutate(RWTbiggerthan042 = case_when(
  US2AIRWT > 0.42 ~ 1,
  TRUE ~ 0
))

# LAVolumnIndex
BaselineCohort_08 <- BaselineCohort_08 %>% mutate(LAVolumnIndexbiggerthan34 = case_when(
  US2AILAVolumnIndex > 34 ~ 1,
  TRUE ~ 0
))

# EeRatio
BaselineCohort_08 <- BaselineCohort_08 %>% mutate(EeRatiobiggerthan9 = case_when(
  US2AIEeRatio > 9 ~ 1,
  TRUE ~ 0
))

# TRVelocity
BaselineCohort_08 <- BaselineCohort_08 %>% mutate(TRVelocitybiggerthan28 = case_when(
  US2AITRVelecity > 2.8 ~ 1,
  TRUE ~ 0
))

BaselineCohort_08$US2AI <- BaselineCohort_08$LVmassbiggerthan95 +
  BaselineCohort_08$RWTbiggerthan042 +
  BaselineCohort_08$LAVolumnIndexbiggerthan34 +
  BaselineCohort_08$EeRatiobiggerthan9 +
  BaselineCohort_08$TRVelocitybiggerthan28
#####
# 12052023
#BaselineCohort_08$Combined is the ESC score based on combined 
#BaselineCohort_08$EHR is ESC score based on EHR
#BaselineCohort_08$US2AI is ESC score based on US2AI


# Biomarker BarcodetoProchi mapping
BiomarkerWithProchi <- merge(Biomarker, BarcodetoProchi, by.x=c("Barcode.ID"), by.y = ("Barcode"))
BiomarkerWithProchi_Selected <- subset(BiomarkerWithProchi, select=c(9,11,13,16))
BaselineCohort_08 <- merge(BaselineCohort_08, BiomarkerWithProchi_Selected, by = c("PROCHI"), all.x = TRUE)
# RocheBiomarker BarcodetoProchi mapping
RocheBiomarkerWithProchi <- merge(RocheBiomarker, BarcodetoProchi, by.x=c("Barcode.ID"), by.y = ("Barcode"))
RocheBiomarkerWithProchi_Selected <- subset(RocheBiomarkerWithProchi, select=c(16,8:13))
BaselineCohort_08 <- merge(BaselineCohort_08, RocheBiomarkerWithProchi_Selected, by = c("PROCHI"), all.x = TRUE)

#####################################################
################## table  ###########################
#####################################################
# for combined
#BaselineCohort_08$Value <- BaselineCohort_08$Combined
# for EHR
BaselineCohort_08$Value <- BaselineCohort_08$EHR
# for US2AI
#BaselineCohort_08$Value <- BaselineCohort_08$US2AI

CM <- c(colnames(Biomarker)[c(9,11,13)],colnames(RocheBiomarker)[c(8:13)])
TableResult <- data.frame(c("HFpEF_No",CM))
colnames(TableResult)[1] <- "Biomarker"
TableResult$HFpEFasClinicalDiagnosis <- 0
TableResult$HFpEFanyOneCriteria <- 0
TableResult$HFpEFanyTwoCriteria <- 0
TableResult$HFpEFanyThreeCriteria <- 0
for (i in 0:3) {
  exclude <- BaselineCohort_08[BaselineCohort_08$FinalDiagnosis=="HFpEF" & BaselineCohort_08$Value <i, c("PROCHI")]
  data1 = BaselineCohort_08[!BaselineCohort_08$PROCHI %in% exclude, c("FinalDiagnosis",CM)]
  data1$FinalDiagnosis<-ifelse(data1$FinalDiagnosis=='Control', 0,1)  
  
  # data cleaning
  for (j in 2:dim(data1)[2]) {
    d <- data1[,j]
    if (is.numeric(d)==FALSE) {
      print(j)
      print(colnames(data1)[j])
      data1[,j]<- as.numeric(as.character(data1[,j]))
    }
  }
  
  
  Data <- data1
  l <- dim(Data)[2]
  
  TableResult[1,i+2] <- dim(Data[Data$FinalDiagnosis==1,])[1]
  
  
  # ROC
  
  for (j in 2:l){

    # prepare the data
    d <- Data[,c(1,j)]
    d_full <- BaselineCohort_08[,c(colnames(d)[2])]
    Q <- quantile(d_full, probs = c(.25, .75),na.rm = TRUE)
    iqr <- IQR(d_full, na.rm=TRUE)
    up <- Q[2]+10*iqr
    low <- Q[1]-10*iqr
    elimated <- subset(d, d[,2]>low & d[,2]<up)
    
    # plot 
    #setwd("P:/Project 3634/R/EKOAI/BiomarkerOutlierCut")
    #tiff(paste0(i,"_Criteria",colnames(d)[2],"_outlier.jpg"), width = 700, height = 350)  # open a file to save plot
    #t <- hist(elimated[,2],breaks = 50, main = paste0(colnames(d)[2]," after outlier cut"), xlab = paste0(colnames(d)[2]," value"))
    #text(x=(max(elimated[,2])-min(elimated[,2]))*0.7,y=max(t$counts)*0.7, labels = paste0("up = ",up,"; low=",low))
    #dev.off()
    
    elimated$logValue <- log10(elimated[,2])
    Diagnose <- elimated[,1]
    p <- as.double(as.character(elimated[,3]))
    
    ROC_curve <- pROC::roc(factor(Diagnose),p)
    A <- round(pROC::auc(ROC_curve),2)
    CI <- ci.auc(ROC_curve)
    cp <- cutpointr(elimated[,3],elimated$FinalDiagnosis,method = maximize_metric, metric = youden, na.rm = TRUE)
    ppvnpv <- coords(ROC_curve,cp$optimal_cutpoint,"threshold",ret = c("ppv","npv"))
    cp_optimal=round(10^cp$optimal_cutpoint,2)
    TableResult[j,i+2] <- paste0(A,"; (",round(CI[1],2),"-",round(CI[3],2),"); ",cp_optimal,"; ",round(ppvnpv$ppv,2),"; ",round(ppvnpv$npv,2))
    
    elimated <- elimated %>% mutate(predict = case_when(
      elimated[,2] >= cp_optimal  ~ 1,
      TRUE ~ 0
    ))
    
    confusionMatrix(as.factor(elimated$FinalDiagnosis), as.factor(elimated$predict))
  }
}
write.csv(TableResult, "P:/Project 3634/R/EKOAI/12052023_ESC_with_biomarker_EHR.csv")
