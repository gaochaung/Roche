## clear the workspace
rm(list = ls())

install.packages("tidyverse")
install.packages("vioplot")
install.packages("rstatix")
install.packages("ggpubr")
install.packages("openxlsx")
install.packages("eeptools")
library(eeptools)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(readxl)
library(dplyr)
library(reshape2)
library(vioplot)
library(openxlsx)
library(data.table)


######################################
############# LOAD data ##############
######################################

#### set the working directory and load excel #######
# read the EKO AI image analysis results
setwd("D:/Project-3634/3634_link-6869_08072021_MeasurementWithDates")
getwd()
EKOImageData_update <- read.csv("_MeasurementsDundee_v2.csv")

# read the EKO AI image analysis results update
#setwd("D:/Project-3634/Chuang/Dundee")
#getwd()
#EKOImageData_update <- read.csv("MeasurementsDundeeJan3-v2.csv")

# read the biomarkers
setwd("D:/Project-3634/Biomarkers/Link-6859/Data hic formatted")
getwd()
Biomarker <- read.csv("tt_3634_IVD_Measurement_alldata_hic.csv")

# read the ROCHE biomarkers
setwd("D:/Project-3634/Biomarkers/Link-6859/Data hic formatted")
getwd()
RocheBiomarker <- read.csv("tt_3634_Results_data2.csv")
#change the following column from factor to numeric
RocheBiomarker[,10] <-  as.double(as.character(RocheBiomarker[,10]))  
RocheBiomarker[,11] <-  as.double(as.character(RocheBiomarker[,11]))  
RocheBiomarker[,12] <-  as.double(as.character(RocheBiomarker[,12]))  
RocheBiomarker[,13] <-  as.double(as.character(RocheBiomarker[,13]))  
# read barcodetoprochi
setwd("D:/Project-3634/Biomarkers/tt_3634_HIC_CombinedList")
getwd()
BarcodetoProchi <- read.csv("3634_Barcode_mapping.csv")
# read the clinical data
setwd("D:/Project-3634/Extr_23159-Samples Clinical Data/Biochemistry (Lab Data Restructured)")
getwd()
Lab <- read.csv("Labs_Biochem.csv")

# read the death 
setwd("D:/Project-3634/Extr_23159-Samples Clinical Data/Deaths - CHI")
getwd()
Death <- read.csv("Deaths_CHI.csv")

# read the death NRS
setwd("D:/Project-3634/Extr_23159-Samples Clinical Data/Deaths GRO - NRS")
getwd()
DeathCAUSE <- read.csv("Deaths_NRS.csv")

# read the demography
setwd("D:/Project-3634/Extr_23159-Samples Clinical Data/Demography")
getwd()
Demography <- read.csv("Demography_Current.csv")

# read the diabetes BMI
setwd("D:/Project-3634/Extr_23159-Samples Clinical Data/Diabetes - BMI")
getwd()
DiabetesBMI <- read.csv("Diabetes_BMI.csv")

# read the diabetes BP
setwd("D:/Project-3634/Extr_23159-Samples Clinical Data/Diabetes - BP")
getwd()
DiabetesBP <- read.csv("Diabetes_BP.csv")

# read the Diabetes - Extended Summary
setwd("D:/Project-3634/Extr_23159-Samples Clinical Data/Diabetes - Extended Summary")
getwd()
Diabetes_Ext_Sum <- read.csv("Diabetes_Ext_Sum.csv")


# read the Echo - Tayside From Nov 2014
setwd("D:/Project-3634/Extr_23159-Samples Clinical Data/Echo - Tayside From Nov 2014")
getwd()
Echo_TaysideFrom2014 <- read.csv("Echo_TaysideFrom2014.csv")

# read the Echo - Tayside From Nov 2013
setwd("D:/Project-3634/Extr_23159-Samples Clinical Data/Echo - Tayside to March 2013")
getwd()
Echo_TaysideTo2013 <- read.csv("Echo_TaysideTo2013.csv")

# read the haematology data
setwd("D:/Project-3634/Extr_23159-Samples Clinical Data/Haematology (Lab Data Restructured)")
getwd()
haematologydata <- read.csv("HaematologyRestructured.csv")

# read the prescribing
setwd("D:/Project-3634/Extr_23159-Samples Clinical Data/Prescribing")
getwd()
prescribing<- read.csv("prescribing.csv")
setwd("D:/Project-3634/Extr_23159-Samples Clinical Data/Prescribing/SupportingDataTables")
prescription_items<- read.csv("prescription items.csv")
# need to morge this two to have the fill records

# read the hospital admission
setwd("D:/Project-3634/Extr_23159-Samples Clinical Data/SMR01 - Hospital Admissions")
getwd()
hospital_admission<- read.csv("SMR01_Admissions.csv")

# read the Cancer Register
setwd("D:/Project-3634/Extr_23159-Samples Clinical Data/SMR06 - Cancer Register")
getwd()
Cancer_Register<- read.csv("SMR06_CancerRegister.csv")

# read the HF_latest_Diagnosis
setwd("D:/Project-3634/Extr_23159-Samples Clinical Data/tt_3684_HIC_HF_latest_Diagnosis_3634")
getwd()
HF_latest_Diagnosis<- read.csv("HF_latest_Diagnosis.csv")

# read the Final Cohort_Summary
#setwd("D:/Project-3634/Release 5_Extr_23098/Project3634_Final Cohort_Summary")
#getwd()
#FinalCohort_Summary<- read.csv("Project3634_Final Cohort_Summary_v1.csv")

# read the Final_Cohort_Detail
#setwd("D:/Project-3634/Release 5_Extr_23098/Project3634_Final_Cohort_Detail")
#getwd()
#Final_Cohort_Detail<- read.csv("Project3634_Final_Cohort_Detail.csv")

# read the new Final_Cohort_Detail
setwd("D:/Project-3634/Extr_23175")
getwd()
Final_Cohort_Summary<- read.csv("Cohort_Summary_Revised.csv")

# read the mapping biomarker barcode to prochi to sampledate
setwd("D:/Project-3634/Extr_23175/tt_3634_HIC_CombinedList_with_Date")
getwd()
Mapping_Plasma_Sample<- read.csv("Mapping_Plasma_Sample.csv")

# read the Project3634_Mon_Echo_screening
setwd("D:/Project-3634/Extr_23175/Project3634_ALL_Echo_screening")
getwd()
Project3634_Mon_Echo_screening<- read.csv("Project3634_Mon_Echo_screening.csv")

######################################
############### all EKO ##############
######################################
Echo_TaysideFrom2014$File <- "From2014"
Echo_TaysideTo2013$File <- "To2013"

EKOPart1 <- subset(Echo_TaysideTo2013, select = c(1,9, 61))
EKOPart2 <- subset(Echo_TaysideFrom2014, select = c(1,5, 195))
colnames(EKOPart1)[2] <- "savetime"
EKO <- rbind(EKOPart1, EKOPart2)
EKO <- unique(EKO)      # all the EKO records
colnames(EKO)[2] <- "EchoDate"
EKO$EchoDate <- as.Date(as.character(EKO$EchoDate),format="%d/%m/%Y")

######################################
########### all EKO with Result ######
######################################
# doing some cleaning of Mon's results
Project3634_Mon_Echo_screening_WithoutDuplicates <- subset(Project3634_Mon_Echo_screening, Diagnosis!="")
Project3634_Mon_Echo_screening_WithoutDuplicates[Project3634_Mon_Echo_screening_WithoutDuplicates == "HfpEF"] <- "HFpEF"
Project3634_Mon_Echo_screening_WithoutDuplicates[Project3634_Mon_Echo_screening_WithoutDuplicates == "no HfpEF"] <- "HFrEF"
Project3634_Mon_Echo_screening_WithoutDuplicates[Project3634_Mon_Echo_screening_WithoutDuplicates == "No HFpEF"] <- "HFrEF"
colnames(Project3634_Mon_Echo_screening_WithoutDuplicates)[2] <- "EchoDate"
EchoResult1 <- Project3634_Mon_Echo_screening_WithoutDuplicates
EchoResult1$EchoDate <- as.Date(as.character(EchoResult1$EchoDate),format="%d/%m/%Y")
EchoResult1 <- unique(EchoResult1)

# loading the second 
setwd("D:/Project-3634/R/EKOAI/For Mon and Chim")
getwd()
EchoResult2014<- read_excel("Echo_TaysideFrom2014_Selected-Mon.xlsx")
EchoResult2013<- read_excel("Echo_TaysideTo2013_Selected-Mon.xlsx")
EchoResult2013_Selected <- subset(EchoResult2013, select = c(2:8))
colnames(EchoResult2013_Selected)[2] <- "EchoDate"
colnames(EchoResult2013_Selected)[6] <- "Diagnosis1"
colnames(EchoResult2013_Selected)[7] <- "Diagnosis2"
EchoResult2013_Selected <- EchoResult2013_Selected %>% mutate(Diagnosis = case_when(
  is.na(Diagnosis2)==TRUE ~ Diagnosis1,
  TRUE ~ Diagnosis2
))   # use Diagnosis2 as diagnosis is default, in case na use diagnosis 1

EchoResult2014_Selected <- subset(EchoResult2014, select = c(2:8))
colnames(EchoResult2014_Selected)[2] <- "EchoDate"
colnames(EchoResult2014_Selected)[6] <- "Diagnosis1"
colnames(EchoResult2014_Selected)[7] <- "Diagnosis2"
EchoResult2014_Selected <- EchoResult2014_Selected %>% mutate(Diagnosis = case_when(
  is.na(Diagnosis2)==TRUE ~ Diagnosis1,
  TRUE ~ Diagnosis2
))

EchoResult2 <- rbind(EchoResult2013_Selected, EchoResult2014_Selected)
EchoResult2 <- unique(EchoResult2)
EchoResult2 <- EchoResult2[,c(1,2,8)]

# inlcude Mon's diagnosis on 230 NA results
# sample with both AI and Echo in the range but not able to produce a diagnosis
setwd("D:/Project-3634/R/EKOAI/ProblemCase1")
t <- read.csv("PatientList230.csv")
# build a dataframe to save the results
p <- t[1,"PROCHI"]
print(p)
setwd("D:/Project-3634/R/EKOAI/ProblemCase1")
EKO2013 <- read_excel(paste0(p,".xlsx"),"EKO2013")
EKO2014 <- read_excel(paste0(p,".xlsx"),"EKO2014")
EKO2013 <- EKO2013[,c(1,2,dim(EKO2013)[2])]
EKO2014 <- EKO2014[,c(1,2,dim(EKO2014)[2])]
colnames(EKO2014)[3] <- "Diagnosis"
colnames(EKO2013) <- colnames(EKO2014)
e <- rbind(EKO2013, EKO2014)
e$savetime <- as.Date(as.character(e$savetime),format="%d/%m/%Y")
E <- unique(e)
for (p in t$PROCHI) {
  #print(p)
  setwd("D:/Project-3634/R/EKOAI/ProblemCase1")
  EKO2013 <- read_excel(paste0(p,".xlsx"),"EKO2013")
  EKO2014 <- read_excel(paste0(p,".xlsx"),"EKO2014")
  EKO2013 <- EKO2013[,c(1,2,dim(EKO2013)[2])]
  EKO2014 <- EKO2014[,c(1,2,dim(EKO2014)[2])]
  colnames(EKO2014)[3] <- "Diagnosis"
  colnames(EKO2013) <- colnames(EKO2014)
  e <- rbind(EKO2013, EKO2014)
  e$savetime <- as.Date(as.character(e$savetime),format="%d/%m/%Y")
  e <- unique(e)
  E <- rbind(E,e)
}
EchoResult3 <- unique(E)
colnames(EchoResult3)[2] <- "EchoDate"

# merge ECHO result
colnames(EchoResult1)[3] <- "Diagnosis1"
colnames(EchoResult2)[3] <- "Diagnosis2"
colnames(EchoResult3)[3] <- "Diagnosis3"
EchoResult <- merge(x=EchoResult1, y=EchoResult2, by=c("PROCHI","EchoDate"),all.x = TRUE)
EchoResult <- EchoResult %>% mutate(Diagnosis = case_when(
  is.na(Diagnosis2)==TRUE ~ as.character(Diagnosis1),
  TRUE ~ Diagnosis2
))
EchoResult <- EchoResult[,c(1,2,5)]
EchoResult <- unique(EchoResult)
colnames(EchoResult)[3] <-"Diagnosisinter"

EchoResult <- merge(x=EchoResult, y=EchoResult3, by=c("PROCHI","EchoDate"),all.x = TRUE)
EchoResult <- EchoResult %>% mutate(Diagnosis = case_when(
  is.na(Diagnosis3)==TRUE ~ as.character(Diagnosisinter),
  TRUE ~ Diagnosis3
))
EchoResult <- EchoResult[,c(1,2,5)]
EchoResult <- unique(EchoResult)

# 17/09/2021
# include MOn's ToBeProvided 
setwd("D:/Project-3634/R/EKOAI/EKOToBeProvided")
t <- read_excel("EKO2013ToBeProvided- Mon.xlsx")
EKO2013ToBeProvided <- t[,c(2,3,dim(t)[2])]
t <- read_excel("EKO2014ToBeProvided-Mon.xlsx")
EKO2014ToBeProvided <- t[,c(2,3,dim(t)[2])]
colnames(EKO2013ToBeProvided)[2] <- "EchoDate"
colnames(EKO2013ToBeProvided)[3] <- "Diagnosis"
colnames(EKO2014ToBeProvided)[2] <- "EchoDate"
colnames(EKO2014ToBeProvided)[3] <- "Diagnosis"

EchoResult <- rbind(EchoResult,EKO2013ToBeProvided)
EchoResult <- rbind(EchoResult,EKO2014ToBeProvided)
EchoResult <- unique(EchoResult)

# 22/09/2021
# include MOn's ToBeProvided v2 
setwd("D:/Project-3634/R/EKOAI/EKOToBeProvided")
t <- read_excel("EKO2014ToBeProvided21092021-Mon.xlsx")
EKO2014ToBeProvided <- t[,c(2,3,dim(t)[2])]
colnames(EKO2014ToBeProvided)[2] <- "EchoDate"
colnames(EKO2014ToBeProvided)[3] <- "Diagnosis"
EKO2014ToBeProvided[is.na(EKO2014ToBeProvided$Diagnosis),"Diagnosis"] <- "NA" 

EchoResult <- rbind(EchoResult,EKO2014ToBeProvided)
EchoResult <- unique(EchoResult)

# 22/09/2021
# include MOn's ToBeProvided v3
setwd("D:/Project-3634/R/EKOAI/EKOToBeProvided")
t <- read_excel("EKO2013ToBeProvided22092021-MON.xlsx")
EKO2013ToBeProvided <- t[,c(2,3,dim(t)[2])]
t <- read_excel("EKO2014ToBeProvided22092021-MON.xlsx")
EKO2014ToBeProvided <- t[,c(2,3,dim(t)[2])]

colnames(EKO2013ToBeProvided)[2] <- "EchoDate"
colnames(EKO2013ToBeProvided)[3] <- "Diagnosis"
colnames(EKO2014ToBeProvided)[2] <- "EchoDate"
colnames(EKO2014ToBeProvided)[3] <- "Diagnosis"

EchoResult <- rbind(EchoResult,EKO2013ToBeProvided)
EchoResult <- rbind(EchoResult,EKO2014ToBeProvided)
EchoResult <- unique(EchoResult)

# 22/09/2021
# include MOn's ToBeProvided v4 
setwd("D:/Project-3634/R/EKOAI/EKOToBeProvided")
t <- read_excel("EKO2014ToBeProvided23092021-Mon.xlsx")
EKO2014ToBeProvided <- t[,c(2,3,dim(t)[2])]
colnames(EKO2014ToBeProvided)[2] <- "EchoDate"
colnames(EKO2014ToBeProvided)[3] <- "Diagnosis"

EchoResult <- rbind(EchoResult,EKO2014ToBeProvided)
EchoResult <- unique(EchoResult)

# merge EKO and EchoResult
EKO$Diagnosis1 <- "ToBeProvided"
colnames(EchoResult)[3] <- "Diagnosis2"
EchoResult <- merge(x=EKO, y=EchoResult, by=c("PROCHI","EchoDate"),all.x = TRUE)
EchoResult <- EchoResult %>% mutate(Diagnosis = case_when(
  is.na(Diagnosis2)==TRUE ~ as.character(Diagnosis1),
  TRUE ~ Diagnosis2
))
EchoResult <- EchoResult[,c(1,2,6)]
EchoResult <- unique(EchoResult)

# correct all the Control Patients into COntrol
ControlPatient <- Final_Cohort_Summary[Final_Cohort_Summary$Diagnosis=="Control",]
EchoResult[EchoResult$PROCHI %in% ControlPatient$PROCHI, "Diagnosis"] <- "Control"
EchoResult <- unique(EchoResult)

# 17/02/2022 update on 75 patients
setwd("D:/Project-3634/R/EKOAI/For Mon")
getwd()
t <- read_xlsx("BaselineCohortforLVEFcheck10022022-Mon-LVEF4-Mon-50.xlsx")

for (i in 1:dim(t)[1]){
  tt <- t[i,]
  if (EchoResult[EchoResult$PROCHI==as.character(tt$PROCHI) & EchoResult$EchoDate==as.Date(as.character(tt$EchoDate),format="%Y-%m-%d"),"Diagnosis"] != "Control"){
    EchoResult[EchoResult$PROCHI==as.character(tt$PROCHI) & EchoResult$EchoDate==as.Date(as.character(tt$EchoDate),format="%Y-%m-%d"),"Diagnosis"] <- 'NA'
  }
}

# 17/02/2022 further updates
setwd("D:/Project-3634/R/EKOAI/For Mon")
getwd()
t <- read.csv("BaselineCohortforLVEFcheck10022022-Mon-LVEF4-Mon-50-v2.csv")

for (i in 1:dim(t)[1]){
  tt <- t[i,]
  if (EchoResult[EchoResult$PROCHI==as.character(tt$PROCHI) & EchoResult$EchoDate==as.Date(as.character(tt$NextEchoDate),format="%Y-%m-%d"),"Diagnosis"] != "Control"){
    EchoResult[EchoResult$PROCHI==as.character(tt$PROCHI) & EchoResult$EchoDate==as.Date(as.character(tt$NextEchoDate),format="%Y-%m-%d"),"Diagnosis"] <- 'NA'
  }
}


###
### EchoResult contains all the EKO and its up to date results provided by Mon
###
load("D:/Project-3634/R/EKOAI/EchoResult.RData")

# 25/08/2021 
######################################
##### all US2.AI with results ########
######################################
EKOImageData_WithDiagnosis <- merge(EKOImageData_update[,c(1,2,4)],Final_Cohort_Summary[,c(1,4)], by.x = c("PatientID"), by.y = c("PROCHI"))
colnames(EKOImageData_WithDiagnosis)[3] <- "US2AILVEF"
EKOImageData_WithDiagnosis$US2AILVEF <- round(EKOImageData_WithDiagnosis$US2AILVEF)
EKOImageData_WithDiagnosis <- EKOImageData_WithDiagnosis %>% mutate(US2AIResult = case_when(
  Diagnosis == "Control" ~ "Control",
  US2AILVEF >= 40 & Diagnosis != "Control" ~ "HFpEF",
  US2AILVEF < 40 & Diagnosis != "Control" ~ "HFrEF"
))
US2AIResult <- unique(EKOImageData_WithDiagnosis[,c(1,2,5)])
colnames(US2AIResult)[1] <- "PROCHI"
colnames(US2AIResult)[2] <- "US2AIDate"
US2AIResult$US2AIDate <- as.Date(as.character(US2AIResult$US2AIDate),format="%d/%m/%Y")

# 06/09/2021
######################################
############## CHart #################
######################################
L1 <- unique(EchoResult$PROCHI)
L2 <- unique(US2AIResult$PROCHI)
L3 <- unique(Mapping_Plasma_Sample$PROCHI)
BiomarkerWithProchi <- merge(Biomarker, BarcodetoProchi, by.x=c("Barcode.ID"), by.y = ("Barcode"))
L4 <- unique(BiomarkerWithProchi$PROCHI)
L <- Reduce(intersect, list(L1, L2, L3))  # 1107 patients

Cohort <- merge(Mapping_Plasma_Sample, Final_Cohort_Summary[,c(1,4)], by =c("PROCHI"))
colnames(Cohort)[4] <- "AssignedDiagnosis"
Cohort$SampleDate <- as.Date(as.character(Cohort$SampleDate),format="%d/%m/%Y")
# records with us2ai
L2 <- as.data.frame(L2)
Cohort <- merge(Cohort, L2, by.x =c("PROCHI"), by.y = c("L2"))
# records with biomarker
L4 <- as.data.frame(L4)
Cohort <- merge(Cohort, L4, by.x =c("PROCHI"), by.y = c("L4"))
Cohort$IsPaired <- ""
Cohort$PairNumber <- 0
Cohort$Case <- 0
Cohort$SelectedPairNumber <- 0
Cohort$EchoDate <- as.Date(as.character("3333-03-30"),format="%Y-%m-%d")
Cohort$EchoDiagnosis <- ""
Cohort$US2AIDate <- as.Date(as.character("3333-03-30"),format="%Y-%m-%d")

for (i in 1:dim(Cohort)[1]){
  p= Cohort[i,"PROCHI"]
  
  ECO_p <- subset(EchoResult, PROCHI==as.character(p))
  US2AI_p <- subset(US2AIResult, PROCHI==as.character(p))
  M <- merge(ECO_p,US2AI_p,by = c("PROCHI"))
  M$TimeDiff <- abs(as.Date(as.character(M$EchoDate),format="%Y-%m-%d")-as.Date(as.character(M$US2AIDate),format="%Y-%m-%d"))
  M_Selected <- subset(M, TimeDiff <= 60)
  if (dim(M_Selected)[1]>0){
    Cohort[i,"IsPaired"] <- "YES"
    Cohort[i,"PairNumber"] <- dim(M_Selected)[1]
    SD <- Cohort[i,"SampleDate"]
    M_Selected$TimeDiff2 <- as.Date(as.character(M_Selected$EchoDate),format="%Y-%m-%d")-as.Date(as.character(SD),format="%Y-%m-%d")
    M_Selected_DateFirst <- subset(M_Selected, TimeDiff2<=0) 
    M_Selected_DateLater1 <- subset(M_Selected, TimeDiff2>0 & TimeDiff2<=180)
    M_Selected_DateLater2 <- subset(M_Selected, TimeDiff2>180)
    #print(p)
    if (dim(M_Selected_DateFirst)[1]>0) {
      #print("case 1")
      Cohort[i,"Case"] <- 1
      Cohort[i,"SelectedPairNumber"] <- dim(M_Selected_DateFirst)[1]
      M_Selected_DateFirst <- M_Selected_DateFirst %>% arrange(-TimeDiff2)
      M_records <- M_Selected_DateFirst[1,]
      Cohort[i,"EchoDate"] <- M_records[1,"EchoDate"]
      Cohort[i,"EchoDiagnosis"] <- M_records[1,"Diagnosis"]
      Cohort[i,"US2AIDate"] <- as.Date(as.character(M_records[1,"US2AIDate"]),format="%Y-%m-%d")
      
    } else if (dim(M_Selected_DateFirst)[1]==0 & dim(M_Selected_DateLater1)[1]>0) {
      #print("case 2")
      Cohort[i,"Case"] <- 2
      Cohort[i,"SelectedPairNumber"] <- dim(M_Selected_DateLater1)[1]
      M_Selected_DateLater1 <- M_Selected_DateLater1 %>% arrange(TimeDiff2)
      M_records <- M_Selected_DateLater1[1,]
      Cohort[i,"EchoDate"] <- M_records[1,"EchoDate"]
      Cohort[i,"EchoDiagnosis"] <- M_records[1,"Diagnosis"]
      Cohort[i,"US2AIDate"] <- as.Date(as.character(M_records[1,"US2AIDate"]),format="%Y-%m-%d")
    } else {
      #print("case 3")
      Cohort[i,"Case"] <- 3
    }
    # # rank from big to samll
    #  # rank from small to big
    #M_Selected <- rbind(M_Selected_DateFirst, M_Selected_DateLater)
    #M_Selected$ID <- 1:nrow(M_Selected)
    #M_Selected <- M_Selected %>% mutate(Within = case_when(
    #  TimeDiff2 <= 180 ~ "YES",
    #  TRUE ~ "NO"
    #))
  } else {
    Cohort[i,"IsPaired"] <- "NO"
  }
    
  
  # SD <- EKO_And_US2AI[EKO_And_US2AI$PROCHI==p,"SampleDate"]
  # M_Selected$TimeDiff2 <- abs(as.Date(as.character(M_Selected$EchoDate),format="%Y-%m-%d")-as.Date(as.character(SD),format="%Y-%m-%d"))
  # M_Selected <- M_Selected %>% top_n(-1, TimeDiff2) #group by the column patientID and select the top 1 of each group odered by Time_Diff
  
}
t <- subset(Cohort, IsPaired == "YES")
t <- subset(Cohort, Case == 3)
table(t$AssignedDiagnosis)
PairedWithInTime <- subset(Cohort, SelectedPairNumber > 0)
t <- PairedWithInTime[PairedWithInTime$Case==2,]
table(t$AssignedDiagnosis)

######################################
####### Deal with NA #################
######################################
NApatient <- subset(PairedWithInTime, PairedWithInTime$EchoDiagnosis == "NA")
NApatient$CompareWithOtherECHO <- "CompareWithOtherECHO"
NApatient$NextEchoDate <- as.Date(as.character("3333-03-30"),format="%Y-%m-%d")
NApatient$ALLECHO <- "NA"
for (p in NApatient$PROCHI) {
  assigned <- NApatient[NApatient$PROCHI==as.character(p),"AssignedDiagnosis"]
  selectedechodate <- NApatient[NApatient$PROCHI==as.character(p),"EchoDate"]
  echo <- EchoResult[EchoResult$PROCHI==as.character(p),]
  echo <-  echo[!echo$EchoDate==selectedechodate,]
  
    echo$Time_Diff <- abs(as.Date(as.character(echo$EchoDate),format="%Y-%m-%d")-as.Date(as.character(selectedechodate),format="%Y-%m-%d"))
    echo <- echo %>% arrange(echo$Time_Diff)  # rank from small to big
    NApatient[NApatient$PROCHI==as.character(p),"ALLECHO"] <- paste0(echo$EchoDate,echo$Diagnosis,sep=", ", collapse="")
    s <- echo[1,]
  
  if (s[1,"Diagnosis"]==assigned) {
    NApatient[NApatient$PROCHI==as.character(p),"CompareWithOtherECHO"] <- "Agree"
    NApatient[NApatient$PROCHI==as.character(p),"NextEchoDate"] <- s[1,"EchoDate"]
  } else if (s[1,"Diagnosis"]=="ToBeProvided"){
    NApatient[NApatient$PROCHI==as.character(p),"CompareWithOtherECHO"] <- "ToBeProvided"
    NApatient[NApatient$PROCHI==as.character(p),"NextEchoDate"] <- s[1,"EchoDate"]
  } else {
    NApatient[NApatient$PROCHI==as.character(p),"CompareWithOtherECHO"] <- "Disagree"
    NApatient[NApatient$PROCHI==as.character(p),"NextEchoDate"] <- s[1,"EchoDate"]
    
  }
}

NApatientAgree <- NApatient[NApatient$CompareWithOtherECHO=="Agree",]
NApatientAgree$FinalDiagnosis <- NApatientAgree$AssignedDiagnosis

PairedWithInTime$FinalDiagnosis <- PairedWithInTime$EchoDiagnosis
PairedWithInTimeSelected <- PairedWithInTime[PairedWithInTime$EchoDiagnosis=="Control"|PairedWithInTime$EchoDiagnosis=="HFpEF"|PairedWithInTime$EchoDiagnosis=="HFrEF",]
PairedWithInTimeSelected$NextEchoDate <- PairedWithInTimeSelected$EchoDate

######################################
#### Final Diagnosis update ##########
######################################
# Final Diagnosis update
setwd("D:/Project-3634/R/EKOAI/For Mon")
getwd()
t <- read_xlsx("BaselineCohortforLVEFcheck10022022-Mon-Chim.xlsx")
t <- t[,c(2:14)]
colnames(t)[12] <- "Group"
t <- t[t$Group!=4,]
t <- t[,c(1:11,13)]
setwd("D:/Project-3634/R/EKOAI/For Mon")
getwd()
tt <- read_xlsx("BaselineCohortforLVEFcheck10022022-Mon-LVEF4-Mon-Chim.xlsx")
tt <- tt[,c(2:9,11:13,19)]
colnames(tt)[10] <- "LVEF simpson"
final_check <- rbind(t,tt)

final_check <- final_check %>% mutate(simpson_available = case_when(
  `LVEF simpson` == 0 | `LVEF simpson` == "NA" |  `LVEF simpson` == "NULL" ~ "No",
  is.na(`LVEF simpson`) ~ "No",
  TRUE ~ "Yes"
))

final_check <- final_check %>% mutate(LVEF = case_when(
  `LVEF simpson` == 0 | `LVEF simpson` == "NA" |  `LVEF simpson` == "NULL" ~ as.double(`Chim LVEF`),
  is.na(`LVEF simpson`) ~ as.double(`Chim LVEF`),
  TRUE ~ as.double(`LVEF simpson`)
))

final_check <- final_check %>% mutate(ClinicalDiagnosis = case_when(
  FinalDiagnosis=="Control" ~ "Control",
  LVEF>0 & LVEF <40 & FinalDiagnosis!="Control"~ "HFrEF",
  LVEF>=40 & LVEF <50 & FinalDiagnosis!="Control"~ "HFmrEF",
  LVEF>=50 & FinalDiagnosis!="Control" ~ "HFpEF"
))

final_check <- final_check %>% mutate(simpsonDiagnosis = case_when(
  FinalDiagnosis=="Control" ~ "Control",
  simpson_available=="Yes" & LVEF>0 & LVEF <40 & FinalDiagnosis!="Control"~ "HFrEF",
  simpson_available=="Yes" & LVEF>=40 & LVEF <50 & FinalDiagnosis!="Control"~ "HFmrEF",
  simpson_available=="Yes" & LVEF>=50 & FinalDiagnosis!="Control" ~ "HFpEF",
  TRUE ~ "NA"
))

final_check <- final_check %>% mutate(chimLVEFDiagnosis = case_when(
  FinalDiagnosis=="Control" ~ "Control",
  as.double(`Chim LVEF`)>0 & as.double(`Chim LVEF`) <40 & FinalDiagnosis!="Control"~ "HFrEF",
  as.double(`Chim LVEF`)>=40 & as.double(`Chim LVEF`) <50 & FinalDiagnosis!="Control"~ "HFmrEF",
  as.double(`Chim LVEF`)>=50 & FinalDiagnosis!="Control" ~ "HFpEF",
  TRUE ~ "NA"
))

# 18/02/2022
####################################################
#### update BaselineCohort Final Diagnosis##########
####################################################
#####################################################
BaselineCohort <- rbind(NApatientAgree[,c(1,2,3,9,11,13,15)],PairedWithInTimeSelected[,c(1,2,3,9,11,13,12)])
BaselineCohort <- merge(BaselineCohort,final_check[,c("PROCHI","LVEF")],by="PROCHI", all.x=TRUE)
BaselineCohort$FinalDiagnosis <- as.character(BaselineCohort$FinalDiagnosis)
BaselineCohort[BaselineCohort$PROCHI %in% final_check$PROCHI, "FinalDiagnosis"] <- final_check$ClinicalDiagnosis
BaselineCohort$FinalDiagnosis <- as.factor(BaselineCohort$FinalDiagnosis)


