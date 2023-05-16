
# 17/09/2021
#####################################################
##################### demograph  ####################
#####################################################
BaselineCohort <- rbind(NApatientAgree[,c(1,2,3,9,11,13,15)],PairedWithInTimeSelected[,c(1,2,3,9,11,13,12)])
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
BaselineCohort[BaselineCohort$FinalDiagnosis=="HFmrEF","FinalDiagnosis"] <- "HFrEF"
######################################
########## Demo Summary ##############
######################################
# process the baseline information for prescribing
n <- c("Control", "HFrEF", "HFpEF")
CN <- c("Age","Sex(female)")
Demo_Summary <- data.frame(CN)
colnames(Demo_Summary)[1] <- ""
Demo_Summary$Control <- 0
Demo_Summary$HFrEF <- 0
Demo_Summary$HFpEF <- 0
Demo_Summary$p <- 0
for (i in 1:3) {
  t <- subset(BaselineCohort, BaselineCohort$FinalDiagnosis == n[i])
  a <- paste0(median(t$Calculated_Age),"±",round(sd(t$Calculated_Age),0)) 
  a <- c(a,paste0(round(table(t$sex)[1]/dim(t)[1],2)*100,"%"))
  a <- data.frame(a)
  Demo_Summary[,i+1] <- a
}

t <- aov(BaselineCohort$Calculated_Age ~ BaselineCohort$FinalDiagnosis, data = BaselineCohort)
Demo_Summary[1,"p"] <- summary(t)[[1]][["Pr(>F)"]][1]

t <- chisq.test(BaselineCohort$sex, BaselineCohort$FinalDiagnosis, correct = FALSE)
Demo_Summary[2,"p"] <- t$p.value

write.csv(Demo_Summary, "P:/Project 3634/R/EKOAI/Demo_Summary_11052023.csv")

#####################################################
########## biochemistry Selected ####################
#####################################################
# include MOn's ToBeProvided 
setwd("P:/Project 3634/R/EKOAI/For Mon")
RCList <- read_excel("EKOAI.test_fre3-Mon-CG.xlsx")
RC <- na.omit(unique(RCList$Var1))
RC <- c(RC,"46N3.")
# all urine proteins does not have read code with it

LabClean <- Lab
# add 44CB. to Transferrin
LabClean[LabClean$ReadCodeDescription=="Transferrin","ReadCodeValue"] <- "44CB."
# add 46N3. to "URINE PROTEINS g/L" and "URINE PROTEINS mg/L"
LabClean[LabClean$ReadCodeDescription=="URINE PROTEINS g/L","ReadCodeValue"] <- "46N3."
LabClean[LabClean$ReadCodeDescription=="URINE PROTEINS g/L","QuantityUnit"] <- "g/L"

LabClean[LabClean$ReadCodeDescription=="URINE PROTEINS mg/L","ReadCodeValue"] <- "46N3."
LabClean[LabClean$ReadCodeDescription=="URINE PROTEINS mg/L","QuantityUnit"] <- "mg/L"


Lab_WithSampleDate <-merge(LabClean[,c(1,9, 10:17, 20, 21)], Mapping_Plasma_Sample, by=c("PROCHI"))
Lab_Selected <- Lab_WithSampleDate[Lab_WithSampleDate$ReadCodeValue %in% RC, ]
Lab_Selected$Time_Diff <- as.Date(as.character(Lab_Selected$DateTimeSampled),format="%d/%m/%Y")-as.Date(as.character(Lab_Selected$SampleDate),format="%d/%m/%Y")
Lab_Selected <- Lab_Selected[Lab_Selected$Time_Diff<=180,]
Lab_Selected <- Lab_Selected[Lab_Selected$PROCHI %in% BaselineCohort$PROCHI,]
Lab_Selected$Time_Diff2 <- abs(Lab_Selected$Time_Diff)

Lab_Selected2 <- Lab_Selected %>% group_by(PROCHI, ReadCodeValue) %>% top_n(-1, Time_Diff2) #group by the column patientID and select the top 1 of each group odered by Time_Diff
Lab_Selected2 <- merge(Lab_Selected2, BaselineCohort, by = c("PROCHI"))

Lab_Sanity <- data.frame(RC)
colnames(Lab_Sanity)[1] <- "ReadCodeValue"
Lab_Sanity$UnitNO <- 0
Lab_Sanity$Unit <- "Read"

for (rc in RC) {
  labs <- Lab_Selected2[Lab_Selected2$ReadCodeValue==rc,] 
  Lab_Sanity[Lab_Sanity$ReadCodeValue==rc,"UnitNO"]= length(unique(labs$QuantityUnit))
  Lab_Sanity[Lab_Sanity$ReadCodeValue==rc,"Unit"]= paste0(unique(labs$QuantityUnit),sep=", ", collapse="")

}
# 42R4. has ug/L and ng/mL I think they are the same
# 46N3. has g/L and mg/L should change that into mg/L
Lab_Selected3 <- Lab_Selected2
Lab_Selected3[Lab_Selected3$QuantityUnit=="g/L"&Lab_Selected3$ReadCodeValue=="46N3.","QuantityValue"] <- Lab_Selected3[Lab_Selected3$QuantityUnit=="g/L"&Lab_Selected3$ReadCodeValue=="46N3.","QuantityValue"]*1000
Lab_Selected3[Lab_Selected3$ReadCodeValue=="46N3.","ReadCodeDescription"] <- "URINE PROTEINS mg/L"
Lab_Selected3[Lab_Selected3$ReadCodeValue=="46N3.","QuantityUnit"] <- "mg/L"


# 22/09/2021
######################################
########### Lab Summary ##############
######################################
Lab_Summary <- data.frame(RC)
colnames(Lab_Summary)[1] <- "ReadCodeValue"
Lab_Summary$ReadCodeDescription <- "ReadCodeDescription"
Lab_Summary$Unit <- "Unit"
Lab_Summary$CohortCoverage <- 0
Lab_Summary$Control <- 0
Lab_Summary$HFrEF <- 0
Lab_Summary$HFpEF <- 0
Lab_Summary$p <- 0

for (rc in RC) {
  labs <- Lab_Selected3[Lab_Selected3$ReadCodeValue==rc,] 
  labs$ISDuplicate <- duplicated(labs$PROCHI)
  Lab_Summary[Lab_Summary$ReadCodeValue==rc,"ReadCodeDescription"]= as.character(unique(labs$ReadCodeDescription)[1])
  Lab_Summary[Lab_Summary$ReadCodeValue==rc,"Unit"]= paste0(as.character(unique(labs$QuantityUnit)),collapse = ",")
  Lab_Summary[Lab_Summary$ReadCodeValue==rc,"CohortCoverage"]= paste0(round(length(unique(labs$PROCHI))/dim(BaselineCohort)[1],2)*100,"%")
  Lab_Summary[Lab_Summary$ReadCodeValue==rc,"Control"]= paste0(round(median(labs[labs$FinalDiagnosis=="Control","QuantityValue"],na.rm = TRUE),1),"±",round(sd(labs[labs$FinalDiagnosis=="Control","QuantityValue"],na.rm = TRUE),1)) 
  Lab_Summary[Lab_Summary$ReadCodeValue==rc,"HFpEF"]= paste0(round(median(labs[labs$FinalDiagnosis=="HFpEF","QuantityValue"],na.rm = TRUE),1),"±",round(sd(labs[labs$FinalDiagnosis=="HFpEF","QuantityValue"],na.rm = TRUE),1)) 
  Lab_Summary[Lab_Summary$ReadCodeValue==rc,"HFrEF"]= paste0(round(median(labs[labs$FinalDiagnosis=="HFrEF","QuantityValue"],na.rm = TRUE),1),"±",round(sd(labs[labs$FinalDiagnosis=="HFrEF","QuantityValue"],na.rm = TRUE),1)) 
   
  t <- aov(labs$QuantityValue ~ labs$FinalDiagnosis, data = labs)
  Lab_Summary[Lab_Summary$ReadCodeValue==rc,"p"] <- summary(t)[[1]][["Pr(>F)"]][1]
}
#11052023
################################################
######### creatinine calcualte eGFR ############
################################################
rc="44J3."
labs <- Lab_Selected2[Lab_Selected2$ReadCodeValue==rc,] 
labs$AgeatTest <- floor(age_calc(as.Date(as.character(labs$anon_date_of_birth),format="%d/%m/%Y"),as.Date(as.character(labs$DateTimeSampled),format="%d/%m/%Y"),units="years"))
# eGFR
labs <- labs %>% mutate(eGFR = case_when(
  sex == "F" ~ 175*(0.0113*QuantityValue)^(-1.154)*AgeatTest^(-0.203)*0.742,
  sex == "M" ~ 175*(0.0113*QuantityValue)^(-1.154)*AgeatTest^(-0.203)
))

rc="451E."
Lab_Summary[Lab_Summary$ReadCodeValue==rc,"Control"]= paste0(round(median(labs[labs$FinalDiagnosis=="Control","eGFR"],na.rm = TRUE),1),"±",round(sd(labs[labs$FinalDiagnosis=="Control","eGFR"],na.rm = TRUE),1)) 
Lab_Summary[Lab_Summary$ReadCodeValue==rc,"HFpEF"]= paste0(round(median(labs[labs$FinalDiagnosis=="HFpEF","eGFR"],na.rm = TRUE),1),"±",round(sd(labs[labs$FinalDiagnosis=="HFpEF","eGFR"],na.rm = TRUE),1)) 
Lab_Summary[Lab_Summary$ReadCodeValue==rc,"HFrEF"]= paste0(round(median(labs[labs$FinalDiagnosis=="HFrEF","eGFR"],na.rm = TRUE),1),"±",round(sd(labs[labs$FinalDiagnosis=="HFrEF","eGFR"],na.rm = TRUE),1)) 



write.csv(Lab_Summary, "P:/Project 3634/R/EKOAI/Lab_Summary_11052023.csv")


#####################################################
########### precription Selected ####################
#####################################################

# 22/09/2021
######################################
########### FUROSEMIDE ###############
######################################

FUROSEMIDE <- read.csv("P:/Project 3634/R/EKOAI/precribing/furosemide.csv")
DrugList <- FUROSEMIDE 
prescribing_WithSampleDate <-merge(prescribing[,c(1,2,8)], Mapping_Plasma_Sample, by=c("PROCHI"))
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$PROCHI %in% BaselineCohort$PROCHI,]
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$res_seqno %in% DrugList$res_seqno, ]
prescribing_WithSampleDate$Time_Diff <- as.Date(as.character(prescribing_WithSampleDate$prescribed_date),format="%d/%m/%Y")-as.Date(as.character(prescribing_WithSampleDate$SampleDate),format="%d/%m/%Y")
prescribing_WithSampleDate <- prescribing_WithSampleDate[!is.na(prescribing_WithSampleDate$Time_Diff),]
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$Time_Diff<=180,]
DrugPatientList <- data.frame(unique(prescribing_WithSampleDate$PROCHI)) 
colnames(DrugPatientList)[1] <- "PROCHI"

DrugPatientList$FUROSEMIDE <- 1
BaselineCohort <- merge(BaselineCohort, DrugPatientList, by=c("PROCHI"),all.x = TRUE)
BaselineCohort[is.na(BaselineCohort$FUROSEMIDE),"FUROSEMIDE"] <- 0

# 22/09/2021
######################################
########### BUMETANIDE ###############
######################################

BUMETANIDE <- prescription_items[prescription_items$BNF_Description %like% "BUMETANIDE",]

DrugList <- BUMETANIDE 
prescribing_WithSampleDate <-merge(prescribing[,c(1,2,8)], Mapping_Plasma_Sample, by=c("PROCHI"))
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$PROCHI %in% BaselineCohort$PROCHI,]
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$res_seqno %in% DrugList$res_seqno, ]
prescribing_WithSampleDate$Time_Diff <- as.Date(as.character(prescribing_WithSampleDate$prescribed_date),format="%d/%m/%Y")-as.Date(as.character(prescribing_WithSampleDate$SampleDate),format="%d/%m/%Y")
prescribing_WithSampleDate <- prescribing_WithSampleDate[!is.na(prescribing_WithSampleDate$Time_Diff),]
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$Time_Diff<=180,]
DrugPatientList <- data.frame(unique(prescribing_WithSampleDate$PROCHI)) 
colnames(DrugPatientList)[1] <- "PROCHI"

DrugPatientList$BUMETANIDE <- 1
BaselineCohort <- merge(BaselineCohort, DrugPatientList, by=c("PROCHI"),all.x = TRUE)
BaselineCohort[is.na(BaselineCohort$BUMETANIDE),"BUMETANIDE"] <- 0

######################################
######## SPIRONOLACTONE ##############
######################################
# process the prescribing SPIRONOLACTONE 
SPIRONOLACTONE <- prescription_items[prescription_items$BNF_Description %like% "SPIRONOL",]

DrugList <- SPIRONOLACTONE 
prescribing_WithSampleDate <-merge(prescribing[,c(1,2,8)], Mapping_Plasma_Sample, by=c("PROCHI"))
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$PROCHI %in% BaselineCohort$PROCHI,]
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$res_seqno %in% DrugList$res_seqno, ]
prescribing_WithSampleDate$Time_Diff <- as.Date(as.character(prescribing_WithSampleDate$prescribed_date),format="%d/%m/%Y")-as.Date(as.character(prescribing_WithSampleDate$SampleDate),format="%d/%m/%Y")
prescribing_WithSampleDate <- prescribing_WithSampleDate[!is.na(prescribing_WithSampleDate$Time_Diff),]
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$Time_Diff<=180,]
DrugPatientList <- data.frame(unique(prescribing_WithSampleDate$PROCHI)) 
colnames(DrugPatientList)[1] <- "PROCHI"

DrugPatientList$SPIRONOLACTONE <- 1
BaselineCohort <- merge(BaselineCohort, DrugPatientList, by=c("PROCHI"),all.x = TRUE)
BaselineCohort[is.na(BaselineCohort$SPIRONOLACTONE),"SPIRONOLACTONE"] <- 0

######################################
########### EPLERENONE ###############
######################################
# process the prescribing EPLERENONE 
EPLERENONE <- prescription_items[prescription_items$BNF_Description %like% "EPLERENONE",]

DrugList <- EPLERENONE 
prescribing_WithSampleDate <-merge(prescribing[,c(1,2,8)], Mapping_Plasma_Sample, by=c("PROCHI"))
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$PROCHI %in% BaselineCohort$PROCHI,]
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$res_seqno %in% DrugList$res_seqno, ]
prescribing_WithSampleDate$Time_Diff <- as.Date(as.character(prescribing_WithSampleDate$prescribed_date),format="%d/%m/%Y")-as.Date(as.character(prescribing_WithSampleDate$SampleDate),format="%d/%m/%Y")
prescribing_WithSampleDate <- prescribing_WithSampleDate[!is.na(prescribing_WithSampleDate$Time_Diff),]
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$Time_Diff<=180,]
DrugPatientList <- data.frame(unique(prescribing_WithSampleDate$PROCHI)) 
colnames(DrugPatientList)[1] <- "PROCHI"

DrugPatientList$EPLERENONE <- 1
BaselineCohort <- merge(BaselineCohort, DrugPatientList, by=c("PROCHI"),all.x = TRUE)
BaselineCohort[is.na(BaselineCohort$EPLERENONE),"EPLERENONE"] <- 0

######################################
############## ASPIRIN ###############
######################################
# process the prescribing ASPIRIN
ASPIRIN <- prescription_items[prescription_items$BNF_Description %like% "ASPIRIN",]

DrugList <- ASPIRIN 
prescribing_WithSampleDate <-merge(prescribing[,c(1,2,8)], Mapping_Plasma_Sample, by=c("PROCHI"))
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$PROCHI %in% BaselineCohort$PROCHI,]
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$res_seqno %in% DrugList$res_seqno, ]
prescribing_WithSampleDate$Time_Diff <- as.Date(as.character(prescribing_WithSampleDate$prescribed_date),format="%d/%m/%Y")-as.Date(as.character(prescribing_WithSampleDate$SampleDate),format="%d/%m/%Y")
prescribing_WithSampleDate <- prescribing_WithSampleDate[!is.na(prescribing_WithSampleDate$Time_Diff),]
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$Time_Diff<=180,]
DrugPatientList <- data.frame(unique(prescribing_WithSampleDate$PROCHI)) 
colnames(DrugPatientList)[1] <- "PROCHI"

DrugPatientList$ASPIRIN <- 1
BaselineCohort <- merge(BaselineCohort, DrugPatientList, by=c("PROCHI"),all.x = TRUE)
BaselineCohort[is.na(BaselineCohort$ASPIRIN),"ASPIRIN"] <- 0

######################################
############## STATIN ################
######################################
# process the prescribing STATIN
STATIN <- prescription_items[prescription_items$BNF_Description %like% "PRAVASTATIN",]
t <- prescription_items[prescription_items$BNF_Description %like% "CERIVASTATIN",]
STATIN <- rbind(STATIN, t)
t <- prescription_items[prescription_items$BNF_Description %like% "SIMVASTATIN",]
STATIN <- rbind(STATIN, t)
t <- prescription_items[prescription_items$BNF_Description %like% "FLUVASTATIN",]
STATIN <- rbind(STATIN, t)
t <- prescription_items[prescription_items$BNF_Description %like% "ATORVASTATIN",]
STATIN <- rbind(STATIN, t)
t <- prescription_items[prescription_items$BNF_Description %like% "ROSUVASTATIN",]
STATIN <- rbind(STATIN, t)

DrugList <- STATIN 
prescribing_WithSampleDate <-merge(prescribing[,c(1,2,8)], Mapping_Plasma_Sample, by=c("PROCHI"))
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$PROCHI %in% BaselineCohort$PROCHI,]
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$res_seqno %in% DrugList$res_seqno, ]
prescribing_WithSampleDate$Time_Diff <- as.Date(as.character(prescribing_WithSampleDate$prescribed_date),format="%d/%m/%Y")-as.Date(as.character(prescribing_WithSampleDate$SampleDate),format="%d/%m/%Y")
prescribing_WithSampleDate <- prescribing_WithSampleDate[!is.na(prescribing_WithSampleDate$Time_Diff),]
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$Time_Diff<=180,]
DrugPatientList <- data.frame(unique(prescribing_WithSampleDate$PROCHI)) 
colnames(DrugPatientList)[1] <- "PROCHI"

DrugPatientList$STATIN <- 1
BaselineCohort <- merge(BaselineCohort, DrugPatientList, by=c("PROCHI"),all.x = TRUE)
BaselineCohort[is.na(BaselineCohort$STATIN),"STATIN"] <- 0

######################################
############# METFORMIN ###############
######################################
# process the prescribing METFORMIN
METFORMIN <-prescription_items[prescription_items$BNF_Description %like% "METFORMIN",]

DrugList <- METFORMIN 
prescribing_WithSampleDate <-merge(prescribing[,c(1,2,8)], Mapping_Plasma_Sample, by=c("PROCHI"))
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$PROCHI %in% BaselineCohort$PROCHI,]
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$res_seqno %in% DrugList$res_seqno, ]
prescribing_WithSampleDate$Time_Diff <- as.Date(as.character(prescribing_WithSampleDate$prescribed_date),format="%d/%m/%Y")-as.Date(as.character(prescribing_WithSampleDate$SampleDate),format="%d/%m/%Y")
prescribing_WithSampleDate <- prescribing_WithSampleDate[!is.na(prescribing_WithSampleDate$Time_Diff),]
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$Time_Diff<=180,]
DrugPatientList <- data.frame(unique(prescribing_WithSampleDate$PROCHI)) 
colnames(DrugPatientList)[1] <- "PROCHI"

DrugPatientList$METFORMIN <- 1
BaselineCohort <- merge(BaselineCohort, DrugPatientList, by=c("PROCHI"),all.x = TRUE)
BaselineCohort[is.na(BaselineCohort$METFORMIN),"METFORMIN"] <- 0

######################################
############## DPP4I #################
######################################
# process the prescribing DPP4I 
DPP4I <-prescription_items[prescription_items$BNF_Description %like% "VILDAGLIPTIN",]
t <- prescription_items[prescription_items$BNF_Description %like% "SAXAGLIPTIN",]
DPP4I <- rbind(DPP4I, t)
t <- prescription_items[prescription_items$BNF_Description %like% "LINAGLIPTIN",]
DPP4I <- rbind(DPP4I, t)
t <- prescription_items[prescription_items$BNF_Description %like% "ALOGLIPTIN",]
DPP4I <- rbind(DPP4I, t)

DrugList <- DPP4I 
prescribing_WithSampleDate <-merge(prescribing[,c(1,2,8)], Mapping_Plasma_Sample, by=c("PROCHI"))
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$PROCHI %in% BaselineCohort$PROCHI,]
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$res_seqno %in% DrugList$res_seqno, ]
prescribing_WithSampleDate$Time_Diff <- as.Date(as.character(prescribing_WithSampleDate$prescribed_date),format="%d/%m/%Y")-as.Date(as.character(prescribing_WithSampleDate$SampleDate),format="%d/%m/%Y")
prescribing_WithSampleDate <- prescribing_WithSampleDate[!is.na(prescribing_WithSampleDate$Time_Diff),]
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$Time_Diff<=180,]
DrugPatientList <- data.frame(unique(prescribing_WithSampleDate$PROCHI)) 
colnames(DrugPatientList)[1] <- "PROCHI"

DrugPatientList$DPP4I <- 1
BaselineCohort <- merge(BaselineCohort, DrugPatientList, by=c("PROCHI"),all.x = TRUE)
BaselineCohort[is.na(BaselineCohort$DPP4I),"DPP4I"] <- 0

######################################
############## SGLT2I ################
######################################
# process the prescribing SGLT2I
SGLT2I <-prescription_items[prescription_items$BNF_Description %like% "DAPAGLIFLOZIN",]
t <- prescription_items[prescription_items$BNF_Description %like% "CANAGLIFLOZIN",]
SGLT2I <- rbind(SGLT2I, t)
t <- prescription_items[prescription_items$BNF_Description %like% "EMPAGLIFLOZIN",]
SGLT2I <- rbind(SGLT2I, t)

DrugList <- SGLT2I 
prescribing_WithSampleDate <-merge(prescribing[,c(1,2,8)], Mapping_Plasma_Sample, by=c("PROCHI"))
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$PROCHI %in% BaselineCohort$PROCHI,]
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$res_seqno %in% DrugList$res_seqno, ]
prescribing_WithSampleDate$Time_Diff <- as.Date(as.character(prescribing_WithSampleDate$prescribed_date),format="%d/%m/%Y")-as.Date(as.character(prescribing_WithSampleDate$SampleDate),format="%d/%m/%Y")
prescribing_WithSampleDate <- prescribing_WithSampleDate[!is.na(prescribing_WithSampleDate$Time_Diff),]
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$Time_Diff<=180,]
DrugPatientList <- data.frame(unique(prescribing_WithSampleDate$PROCHI)) 
colnames(DrugPatientList)[1] <- "PROCHI"

DrugPatientList$SGLT2I <- 1
BaselineCohort <- merge(BaselineCohort, DrugPatientList, by=c("PROCHI"),all.x = TRUE)
BaselineCohort[is.na(BaselineCohort$SGLT2I),"SGLT2I"] <- 0

######################################
############# insulin ################
######################################
# process the prescribing SGLT2I
insulin <-prescription_items[prescription_items$BNF_Description %like% "HUMALOG",]
t <- prescription_items[prescription_items$Approved_Name %like% "INSULIN",]
insulin <- rbind(insulin, t)


DrugList <- insulin 
prescribing_WithSampleDate <-merge(prescribing[,c(1,2,8)], Mapping_Plasma_Sample, by=c("PROCHI"))
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$PROCHI %in% BaselineCohort$PROCHI,]
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$res_seqno %in% DrugList$res_seqno, ]
prescribing_WithSampleDate$Time_Diff <- as.Date(as.character(prescribing_WithSampleDate$prescribed_date),format="%d/%m/%Y")-as.Date(as.character(prescribing_WithSampleDate$SampleDate),format="%d/%m/%Y")
prescribing_WithSampleDate <- prescribing_WithSampleDate[!is.na(prescribing_WithSampleDate$Time_Diff),]
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$Time_Diff<=180,]
DrugPatientList <- data.frame(unique(prescribing_WithSampleDate$PROCHI)) 
colnames(DrugPatientList)[1] <- "PROCHI"

DrugPatientList$insulin <- 1
BaselineCohort <- merge(BaselineCohort, DrugPatientList, by=c("PROCHI"),all.x = TRUE)
BaselineCohort[is.na(BaselineCohort$insulin),"insulin"] <- 0


######################################
############## Beta blocker #################
######################################
# process the prescribing Beta blocker
 
Betablocker<-prescription_items[prescription_items$BNF_Description %like%"SOTALOL",]
t <- prescription_items[prescription_items$BNF_Description %like% "PROPRANOLOL",]
Betablocker <- rbind(Betablocker, t)
t <- prescription_items[prescription_items$BNF_Description %like% "CARVEDILOL",]
Betablocker <- rbind(Betablocker, t)
t <- prescription_items[prescription_items$BNF_Description %like% "NEBIVOLOL",]
Betablocker <- rbind(Betablocker, t)
t <- prescription_items[prescription_items$BNF_Description %like% "LABETALOL",]
Betablocker <- rbind(Betablocker, t)
t <- prescription_items[prescription_items$BNF_Description %like% "OXPRENOLOL",]
Betablocker <- rbind(Betablocker, t)
t <- prescription_items[prescription_items$BNF_Description %like% "CELIPROLOL",]
Betablocker <- rbind(Betablocker, t)
t <- prescription_items[prescription_items$BNF_Description %like% "ATENOLOL",]
Betablocker <- rbind(Betablocker, t)
t <- prescription_items[prescription_items$BNF_Description %like% "METOPROLOL",]
Betablocker <- rbind(Betablocker, t)
t <- prescription_items[prescription_items$BNF_Description %like% "BISOPROLOL",]
Betablocker <- rbind(Betablocker, t)
t <- prescription_items[prescription_items$BNF_Description %like% "PINDOLOL",]
Betablocker <- rbind(Betablocker, t)
t <- prescription_items[prescription_items$BNF_Description %like% "NEBIVOLOL",]
Betablocker <- rbind(Betablocker, t)


DrugList <- Betablocker
prescribing_WithSampleDate <-merge(prescribing[,c(1,2,8)], Mapping_Plasma_Sample, by=c("PROCHI"))
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$PROCHI %in% BaselineCohort$PROCHI,]
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$res_seqno %in% DrugList$res_seqno, ]
prescribing_WithSampleDate$Time_Diff <- as.Date(as.character(prescribing_WithSampleDate$prescribed_date),format="%d/%m/%Y")-as.Date(as.character(prescribing_WithSampleDate$SampleDate),format="%d/%m/%Y")
prescribing_WithSampleDate <- prescribing_WithSampleDate[!is.na(prescribing_WithSampleDate$Time_Diff),]
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$Time_Diff<=180,]
DrugPatientList <- data.frame(unique(prescribing_WithSampleDate$PROCHI)) 
colnames(DrugPatientList)[1] <- "PROCHI"

DrugPatientList$Betablocker <- 1
BaselineCohort <- merge(BaselineCohort, DrugPatientList, by=c("PROCHI"),all.x = TRUE)
BaselineCohort[is.na(BaselineCohort$Betablocker),"Betablocker"] <- 0


######################################
############## ARB ################
######################################
# process the prescribing ARB
ARB <-prescription_items[prescription_items$BNF_Description %like% "IRBESARTAN",]
t <- prescription_items[prescription_items$BNF_Description %like% "CANDESARTAN",]
ARB <- rbind(ARB, t)
t <- prescription_items[prescription_items$BNF_Description %like% "TELMISARTAN",]
ARB <- rbind(ARB, t)
t <- prescription_items[prescription_items$BNF_Description %like% "EPROSARTAN",]
ARB <- rbind(ARB, t)
t <- prescription_items[prescription_items$BNF_Description %like% "LOSARTAN",]
ARB <- rbind(ARB, t)
t <- prescription_items[prescription_items$BNF_Description %like% "OLMESARTAN",]
ARB <- rbind(ARB, t)
t <- prescription_items[prescription_items$BNF_Description %like% "VALSARTAN",]
ARB <- rbind(ARB, t)
t <- prescription_items[prescription_items$BNF_Description %like% "AZILSARTAN",]
ARB <- rbind(ARB, t)
t <- prescription_items[prescription_items$BNF_Description %like% "IRBESARTAN",]
ARB <- rbind(ARB, t)
t <- prescription_items[prescription_items$BNF_Description %like% "EPROSARTAN",]



DrugList <- ARB
prescribing_WithSampleDate <-merge(prescribing[,c(1,2,8)], Mapping_Plasma_Sample, by=c("PROCHI"))
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$PROCHI %in% BaselineCohort$PROCHI,]
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$res_seqno %in% DrugList$res_seqno, ]
prescribing_WithSampleDate$Time_Diff <- as.Date(as.character(prescribing_WithSampleDate$prescribed_date),format="%d/%m/%Y")-as.Date(as.character(prescribing_WithSampleDate$SampleDate),format="%d/%m/%Y")
prescribing_WithSampleDate <- prescribing_WithSampleDate[!is.na(prescribing_WithSampleDate$Time_Diff),]
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$Time_Diff<=180,]
DrugPatientList <- data.frame(unique(prescribing_WithSampleDate$PROCHI)) 
colnames(DrugPatientList)[1] <- "PROCHI"

DrugPatientList$ARB <- 1
BaselineCohort <- merge(BaselineCohort, DrugPatientList, by=c("PROCHI"),all.x = TRUE)
BaselineCohort[is.na(BaselineCohort$ARB),"ARB"] <- 0


######################################
############## ACEI ################
######################################
# process the prescribing ACEI

ACEI<-prescription_items[prescription_items$BNF_Description %like% "CAPTOPRIL",]
t <- prescription_items[prescription_items$BNF_Description %like% "ENALAPRIL",]
ACEI<- rbind(ACEI, t)
t <- prescription_items[prescription_items$BNF_Description %like% "RAMIPRIL",]
ACEI<- rbind(ACEI, t)
t <- prescription_items[prescription_items$BNF_Description %like% "IMIDAPRIL",]
ACEI<- rbind(ACEI, t)
t <- prescription_items[prescription_items$BNF_Description %like% "LISINOPRIL",]
ACEI<- rbind(ACEI, t)
t <- prescription_items[prescription_items$BNF_Description %like% "PERINDOPRIL",]
ACEI<- rbind(ACEI, t)
t <- prescription_items[prescription_items$BNF_Description %like% "TRANDOLAPRIL",]
ACEI<- rbind(ACEI, t)
t <- prescription_items[prescription_items$BNF_Description %like% "QUINAPRIL",]
ACEI<- rbind(ACEI, t)
t <- prescription_items[prescription_items$BNF_Description %like% "FOSINOPRIL",]
ACEI<- rbind(ACEI, t)
t <- prescription_items[prescription_items$BNF_Description %like% "CILAZAPRIL",]
ACEI<- rbind(ACEI, t)
t <- prescription_items[prescription_items$BNF_Description %like% "MOEXIPRIL",]
ACEI<- rbind(ACEI, t)

DrugList <- ACEI
prescribing_WithSampleDate <-merge(prescribing[,c(1,2,8)], Mapping_Plasma_Sample, by=c("PROCHI"))
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$PROCHI %in% BaselineCohort$PROCHI,]
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$res_seqno %in% DrugList$res_seqno, ]
prescribing_WithSampleDate$Time_Diff <- as.Date(as.character(prescribing_WithSampleDate$prescribed_date),format="%d/%m/%Y")-as.Date(as.character(prescribing_WithSampleDate$SampleDate),format="%d/%m/%Y")
prescribing_WithSampleDate <- prescribing_WithSampleDate[!is.na(prescribing_WithSampleDate$Time_Diff),]
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$Time_Diff<=180,]
DrugPatientList <- data.frame(unique(prescribing_WithSampleDate$PROCHI)) 
colnames(DrugPatientList)[1] <- "PROCHI"

DrugPatientList$ACEI <- 1
BaselineCohort <- merge(BaselineCohort, DrugPatientList, by=c("PROCHI"),all.x = TRUE)
BaselineCohort[is.na(BaselineCohort$ACEI),"ACEI"] <- 0


######################################
######### Drug Summary ###############
######################################
# process the baseline information for prescribing
n <- c("Control", "HFrEF", "HFpEF")
PI <- c("FUROSEMIDE","BUMETANIDE","SPIRONOLACTONE","EPLERENONE","ASPIRIN","STATIN","METFORMIN","DPP4I","SGLT2I","insulin","Betablocker","ARB","ACEI")
Precribing_Summary <- data.frame(PI)
colnames(Precribing_Summary)[1] <- "PrecribingItem"
Precribing_Summary$Control <- 0
Precribing_Summary$HFrEF <- 0
Precribing_Summary$HFpEF <- 0
Precribing_Summary$p <- 0
for (i in 1:3) {
t <- subset(BaselineCohort, BaselineCohort$FinalDiagnosis == n[i])
a <- paste0(round(table(t$FUROSEMIDE)[2]/dim(t)[1],2)*100,"%") 
a <- c(a,paste0(round(table(t$BUMETANIDE)[2]/dim(t)[1],2)*100,"%"))
a <- c(a,paste0(round(table(t$SPIRONOLACTONE)[2]/dim(t)[1],2)*100,"%"))
a <- c(a,paste0(round(table(t$EPLERENONE)[2]/dim(t)[1],2)*100,"%"))
a <- c(a,paste0(round(table(t$ASPIRIN)[2]/dim(t)[1],2)*100,"%"))
a <- c(a,paste0(round(table(t$STATIN)[2]/dim(t)[1],2)*100,"%"))
a <- c(a,paste0(round(table(t$METFORMIN)[2]/dim(t)[1],2)*100,"%"))
a <- c(a,paste0(round(table(t$DPP4I)[2]/dim(t)[1],2)*100,"%"))
a <- c(a,paste0(round(table(t$SGLT2I)[2]/dim(t)[1],2)*100,"%"))
a <- c(a,paste0(round(table(t$insulin)[2]/dim(t)[1],2)*100,"%"))
a <- c(a,paste0(round(table(t$Betablocker)[2]/dim(t)[1],2)*100,"%"))
a <- c(a,paste0(round(table(t$ARB)[2]/dim(t)[1],2)*100,"%"))
a <- c(a,paste0(round(table(t$ACEI)[2]/dim(t)[1],2)*100,"%"))
a <- data.frame(a)
Precribing_Summary[,i+1] <- a
}
t <- chisq.test(BaselineCohort$FUROSEMIDE, BaselineCohort$FinalDiagnosis, correct = FALSE)
Precribing_Summary[1,"p"] <- t$p.value
t <- chisq.test(BaselineCohort$BUMETANIDE, BaselineCohort$FinalDiagnosis, correct = FALSE)
Precribing_Summary[2,"p"] <- t$p.value
t <- chisq.test(BaselineCohort$SPIRONOLACTONE, BaselineCohort$FinalDiagnosis, correct = FALSE)
Precribing_Summary[3,"p"] <- t$p.value
t <- chisq.test(BaselineCohort$EPLERENONE, BaselineCohort$FinalDiagnosis, correct = FALSE)
Precribing_Summary[4,"p"] <- t$p.value
t <- chisq.test(BaselineCohort$ASPIRIN, BaselineCohort$FinalDiagnosis, correct = FALSE)
Precribing_Summary[5,"p"] <- t$p.value
t <- chisq.test(BaselineCohort$STATIN, BaselineCohort$FinalDiagnosis, correct = FALSE)
Precribing_Summary[6,"p"] <- t$p.value
t <- chisq.test(BaselineCohort$METFORMIN, BaselineCohort$FinalDiagnosis, correct = FALSE)
Precribing_Summary[7,"p"] <- t$p.value
t <- chisq.test(BaselineCohort$DPP4I, BaselineCohort$FinalDiagnosis, correct = FALSE)
Precribing_Summary[8,"p"] <- t$p.value
t <- chisq.test(BaselineCohort$SGLT2I, BaselineCohort$FinalDiagnosis, correct = FALSE)
Precribing_Summary[9,"p"] <- t$p.value
t <- chisq.test(BaselineCohort$insulin, BaselineCohort$FinalDiagnosis, correct = FALSE)
Precribing_Summary[10,"p"] <- t$p.value
t <- chisq.test(BaselineCohort$Betablocker, BaselineCohort$FinalDiagnosis, correct = FALSE)
Precribing_Summary[11,"p"] <- t$p.value
t <- chisq.test(BaselineCohort$ARB, BaselineCohort$FinalDiagnosis, correct = FALSE)
Precribing_Summary[12,"p"] <- t$p.value
t <- chisq.test(BaselineCohort$ACEI, BaselineCohort$FinalDiagnosis, correct = FALSE)
Precribing_Summary[13,"p"] <- t$p.value


write.csv(Precribing_Summary, "P:/Project 3634/R/EKOAI/Precribing_Summary_11052023.csv")

#####################################################
################ Comorbidities ######################
#####################################################
#08/09/2021
######################################
############# Diabetes ###############
######################################
Diabetes_Ext_Sum_Selected <- subset(Diabetes_Ext_Sum, select = c(1,3, 4))
Diabetes_Ext_Sum_Selected_WithSampleDate <-merge(Diabetes_Ext_Sum_Selected, Mapping_Plasma_Sample, by = ("PROCHI"))
Diabetes_Ext_Sum_Selected_WithSampleDate$Time_Diff <- as.Date(as.character(Diabetes_Ext_Sum_Selected_WithSampleDate$DateOfDiagnosisDiabetes_Date),format="%d/%m/%Y")-as.Date(as.character(Diabetes_Ext_Sum_Selected_WithSampleDate$SampleDate),format="%d/%m/%Y")
Diabetes_List <- Diabetes_Ext_Sum_Selected_WithSampleDate[Diabetes_Ext_Sum_Selected_WithSampleDate$Time_Diff<=180,]
#Diabetes_Ext_Sum_Selected2 <- subset(Diabetes_Ext_Sum_Selected_WithSampleDate, Time_Diff <=0, select = c(1,3))
#Diabetes_Ext_Sum_Selected2 <- Diabetes_Ext_Sum_Selected2 %>% mutate(Diabetes = case_when(
#  Diabetes_Ext_Sum_Selected2$DiabetesMellitusType_Mapped == "Impaired Glucose Metabolism and Other" ~ 1,
#  Diabetes_Ext_Sum_Selected2$DiabetesMellitusType_Mapped == "Diabetes in Remission" ~ 2,
#  Diabetes_Ext_Sum_Selected2$DiabetesMellitusType_Mapped == "Type 1 Diabetes Mellitus" ~ 3,
#  Diabetes_Ext_Sum_Selected2$DiabetesMellitusType_Mapped == "Type 2 Diabetes Mellitus" ~ 4,
#  Diabetes_Ext_Sum_Selected2$DiabetesMellitusType_Mapped == "Other" ~ 5
#))
#Diabetes_Ext_Sum_Selected3 <- subset(Diabetes_Ext_Sum_Selected2, select = c(1,3))
Diabetes_List$Diabetes <- 1
BaselineCohort <- merge(BaselineCohort, Diabetes_List[,c(1,7)], by=c("PROCHI"),all.x = TRUE)
BaselineCohort[is.na(BaselineCohort$Diabetes),"Diabetes"] <- 0

# BaselineCohort is the 671 cohort that we will focus on
t <- subset(BaselineCohort, BaselineCohort$FinalDiagnosis == "HFrEF")
print("Diabete percentage")
round(table(t$Diabetes)[2]/dim(t)[1],2)


#08/09/2021
######################################
################ AF ##################
######################################
# process the prescribing STATIN 

# find patients precribed one of the following drags
AF_DrugList <- prescription_items[prescription_items$BNF_Description %like% "EDOXABAN",]   # 3 ITEMS
t <- prescription_items[prescription_items$BNF_Description %like% "APIXABAN",]   # 2 ITEMS
AF_DrugList <- rbind(AF_DrugList, t)
t <- prescription_items[prescription_items$BNF_Description %like% "RIVAROXABAN",]   # 4 ITEMS
AF_DrugList <- rbind(AF_DrugList, t)
t <- prescription_items[prescription_items$BNF_Description %like% "DABIGATRAN",]   # 3 ITEMS
AF_DrugList <- rbind(AF_DrugList, t)
t <- prescription_items[prescription_items$BNF_Description %like% "wARFARIN",]   # 31 ITEMS
AF_DrugList <- rbind(AF_DrugList, t)

DrugList <- AF_DrugList
prescribing_WithSampleDate <-merge(prescribing[,c(1,2,8)], BaselineCohort, by=c("PROCHI"))
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$res_seqno %in% DrugList$res_seqno, ]
prescribing_WithSampleDate$Time_Diff <- as.Date(as.character(prescribing_WithSampleDate$prescribed_date),format="%d/%m/%Y")-as.Date(as.character(prescribing_WithSampleDate$EchoDate),format="%Y-%m-%d")
prescribing_WithSampleDate <- prescribing_WithSampleDate[!is.na(prescribing_WithSampleDate$Time_Diff),]
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$Time_Diff<=0,]
DrugPatientList <- unique(prescribing_WithSampleDate$PROCHI) 

AF_ICD <- c("I480", "I481", "I482", "I483", "I484", "I489", "I48X")
DVT_ICD <- c("I800", "I801", "I802", "I803", "I808", "I809", "I820", "I821", "I822", "I823", "I824", "I825", "I826", "I827", "I828", "I829")
PE_ICD <- c("I260", "I269", "I278", "Z867", "T790", "T791", "T800", "T817")
ProstheticValve_ICD <- c("Z952")

hospital_admission_s_WithSampleDate <-merge(hospital_admission[,c(1,8,17:22)], BaselineCohort, by=c("PROCHI"))
hospital_admission_s_WithSampleDate$Time_Diff <- as.Date(as.character(hospital_admission_s_WithSampleDate$ADMISSION_DATE),format="%d/%m/%Y")-as.Date(as.character(hospital_admission_s_WithSampleDate$EchoDate),format="%Y-%m-%d")
hospital_admission_s_WithSampleDate <- hospital_admission_s_WithSampleDate[!is.na(hospital_admission_s_WithSampleDate$Time_Diff),]
hospital_admission_s_WithSampleDate <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$Time_Diff<=0,]

# AF patients
IcdCode <- AF_ICD

hospital_admission_patient <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$MAIN_CONDITION %in% IcdCode,]
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_1 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_2 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_3 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_4 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_5 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)

hospital_admission_s_WithSampleDate_AF <- hospital_admission_patient
# DVT patients
IcdCode <- DVT_ICD

hospital_admission_patient <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$MAIN_CONDITION %in% IcdCode,]
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_1 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_2 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_3 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_4 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_5 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)

hospital_admission_s_WithSampleDate_DVT <- hospital_admission_patient

# PE patients
IcdCode <- PE_ICD

hospital_admission_patient <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$MAIN_CONDITION %in% IcdCode,]
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_1 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_2 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_3 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_4 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_5 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)

hospital_admission_s_WithSampleDate_PE <- hospital_admission_patient

# ProstheticValve patients
IcdCode <- ProstheticValve_ICD

hospital_admission_patient <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$MAIN_CONDITION %in% IcdCode,]
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_1 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_2 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_3 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_4 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_5 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)

hospital_admission_s_WithSampleDate_ProstheticValve <- hospital_admission_patient


AFPatientList <- unique(hospital_admission_s_WithSampleDate_AF[,1])
DVTPatientList <- unique(hospital_admission_s_WithSampleDate_DVT[,1])
PEPatientList <- unique(hospital_admission_s_WithSampleDate_PE[,1])
ProstheticValvePatientList <- unique(hospital_admission_s_WithSampleDate_ProstheticValve[,1])

C <- c(as.character(AFPatientList),as.character(DrugPatientList))
C <- unique(C)
ToBeExclude <- c(as.character(DVTPatientList),as.character(PEPatientList),as.character(ProstheticValvePatientList))
ToBeExclude <- unique(ToBeExclude)

Method2 <- C[!C %in% ToBeExclude]
Method1 <- as.character(AFPatientList)

AF <- c(Method1,Method2)
AF <- data.frame(unique(AF))
colnames(AF)[1] <- "PROCHI"
AF$AF <- 1
BaselineCohort <- merge(BaselineCohort, AF, by=c("PROCHI"),all.x = TRUE)
BaselineCohort[is.na(BaselineCohort$AF),"AF"] <- 0

# BaselineCohort is the 671 cohort that we will focus on
t <- subset(BaselineCohort, BaselineCohort$FinalDiagnosis == "HFrEF")
print("AF percentage")
round(table(t$AF)[2]/dim(t)[1],2)

#05/10/2021
######################################
############### CAD ##################
######################################
# read ICD Codes 
setwd("P:/Project 3634")
getwd()
CAD_ICD_Codes <- read_excel("ICD code for CAD.xlsx")
IcdCode <- CAD_ICD_Codes$I20.0 

hospital_admission_s_WithSampleDate <-merge(hospital_admission[,c(1,8,17:22)], BaselineCohort, by=c("PROCHI"))
hospital_admission_s_WithSampleDate$Time_Diff <- as.Date(as.character(hospital_admission_s_WithSampleDate$ADMISSION_DATE),format="%d/%m/%Y")-as.Date(as.character(hospital_admission_s_WithSampleDate$EchoDate),format="%Y-%m-%d")
hospital_admission_s_WithSampleDate <- hospital_admission_s_WithSampleDate[!is.na(hospital_admission_s_WithSampleDate$Time_Diff),]
hospital_admission_s_WithSampleDate <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$Time_Diff<=180,]


hospital_admission_patient <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$MAIN_CONDITION %in% IcdCode,]
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_1 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_2 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_3 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_4 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_5 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)

hospital_admission_s_WithSampleDate_CAD <- hospital_admission_patient

CADPatientList <- unique(hospital_admission_s_WithSampleDate_CAD[,1])

CAD <- data.frame(unique(CADPatientList))
colnames(CAD)[1] <- "PROCHI"
CAD$CAD <- 1
BaselineCohort <- merge(BaselineCohort, CAD, by=c("PROCHI"),all.x = TRUE)
BaselineCohort[is.na(BaselineCohort$CAD),"CAD"] <- 0

# BaselineCohort is the 671 cohort that we will focus on
t <- subset(BaselineCohort, BaselineCohort$FinalDiagnosis == "HFrEF")
print("CAD percentage")
round(table(t$CAD)[2]/dim(t)[1],2)

#12/05/2022
######################################
################ CKD #################
######################################
# include MOn's ToBeProvided 
setwd("P:/Project 3634/R/EKOAI/For Mon")
RCList <- read_excel("EKOAI.test_fre3-Mon-CG.xlsx")
RC <- na.omit(unique(RCList$Var1))
RC <- c(RC,"46N3.")
# all urine proteins does not have read code with it

LabClean <- Lab
# add 44CB. to Transferrin
LabClean[LabClean$ReadCodeDescription=="Transferrin","ReadCodeValue"] <- "44CB."
# add 46N3. to "URINE PROTEINS g/L" and "URINE PROTEINS mg/L"
LabClean[LabClean$ReadCodeDescription=="URINE PROTEINS g/L","ReadCodeValue"] <- "46N3."
LabClean[LabClean$ReadCodeDescription=="URINE PROTEINS g/L","QuantityUnit"] <- "g/L"

LabClean[LabClean$ReadCodeDescription=="URINE PROTEINS mg/L","ReadCodeValue"] <- "46N3."
LabClean[LabClean$ReadCodeDescription=="URINE PROTEINS mg/L","QuantityUnit"] <- "mg/L"


Lab_WithSampleDate <-merge(LabClean[,c(1,9, 10:17, 20, 21)], Mapping_Plasma_Sample, by=c("PROCHI"))
Lab_Selected <- Lab_WithSampleDate[Lab_WithSampleDate$ReadCodeValue %in% RC, ]
Lab_Selected$Time_Diff <- as.Date(as.character(Lab_Selected$DateTimeSampled),format="%d/%m/%Y")-as.Date(as.character(Lab_Selected$SampleDate),format="%d/%m/%Y")
Lab_Selected <- Lab_Selected[Lab_Selected$Time_Diff<=180,]
Lab_Selected <- Lab_Selected[Lab_Selected$PROCHI %in% BaselineCohort$PROCHI,]
Lab_Selected$Time_Diff2 <- abs(Lab_Selected$Time_Diff)

Lab_Selected2 <- Lab_Selected %>% group_by(PROCHI, ReadCodeValue) %>% top_n(-1, Time_Diff2) #group by the column patientID and select the top 1 of each group odered by Time_Diff
Lab_Selected2 <- merge(Lab_Selected2, BaselineCohort, by = c("PROCHI"))


# creatinine
rc="44J3."
labs <- Lab_Selected2[Lab_Selected2$ReadCodeValue==rc,] 
labs$AgeatTest <- floor(age_calc(as.Date(as.character(labs$anon_date_of_birth),format="%d/%m/%Y"),as.Date(as.character(labs$DateTimeSampled),format="%d/%m/%Y"),units="years"))
# eGFR
labs <- labs %>% mutate(eGFR = case_when(
  sex == "F" ~ 175*(0.0113*QuantityValue)^(-1.154)*AgeatTest^(-0.203)*0.742,
  sex == "M" ~ 175*(0.0113*QuantityValue)^(-1.154)*AgeatTest^(-0.203)
))

labs$duplicated <- duplicated(labs$PROCHI)
t <- labs[labs$PROCHI %in% unique(labs[labs$duplicated==TRUE,"PROCHI"]),]  #find all duplicate records # find eGFR lower than 60, all those are consider to be CKD

labs_v2 <- labs %>% group_by(PROCHI) %>% summarise_at(vars(eGFR),list(name=mean))
colnames(labs_v2)[2] <- "eGFR"
CKD_List <- labs_v2 %>% mutate(CKD = case_when(
  labs_v2$eGFR <= 59 ~ 1,
  TRUE ~ 0
))

BaselineCohort <- merge(BaselineCohort, CKD_List[,c(1,3)], by=c("PROCHI"),all.x = TRUE)
BaselineCohort[is.na(BaselineCohort$CKD),"CKD"] <- 0
t <- BaselineCohort[BaselineCohort$CKD==1,]
table(t$FinalDiagnosis)
t <- chisq.test(BaselineCohort$CKD, BaselineCohort$FinalDiagnosis, correct = FALSE)
t$p.value


#23/05/2022
######################################
############### COPD #################
######################################
#23/05/2022
################ COPD ##################

# find patients precribed one of the following drags
COPD_DrugList <- prescription_items[prescription_items$BNF_Description %like% "FLUTICASONE/SALMETEROL",]   #fluticasone/salmeterol
t <- prescription_items[prescription_items$BNF_Description %like% "SALBUTAMOL",]   #salbutamol
COPD_DrugList <- rbind(COPD_DrugList, t)
t <- prescription_items[prescription_items$BNF_Description %like% "PREDNISOLONE_TAB",]  #prednisolone_tab
COPD_DrugList <- rbind(COPD_DrugList, t)

DrugList <- COPD_DrugList
prescribing_WithSampleDate <-merge(prescribing[,c(1,2,8)], BaselineCohort, by=c("PROCHI"))
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$res_seqno %in% DrugList$res_seqno, ]
prescribing_WithSampleDate$Time_Diff <- as.Date(as.character(prescribing_WithSampleDate$prescribed_date),format="%d/%m/%Y")-as.Date(as.character(prescribing_WithSampleDate$EchoDate),format="%Y-%m-%d")
prescribing_WithSampleDate <- prescribing_WithSampleDate[!is.na(prescribing_WithSampleDate$Time_Diff),]
prescribing_WithSampleDate <- prescribing_WithSampleDate[prescribing_WithSampleDate$Time_Diff<=0,]
DrugPatientList <- unique(prescribing_WithSampleDate$PROCHI) 

COPD_ICD <- c("J432", "J441", "J440", "J449", "J43", "J439", "J430", "J44", "J438", "J448", "J431")

hospital_admission_s_WithSampleDate <-merge(hospital_admission[,c(1,8,17:22)], BaselineCohort, by=c("PROCHI"))
hospital_admission_s_WithSampleDate$Time_Diff <- as.Date(as.character(hospital_admission_s_WithSampleDate$ADMISSION_DATE),format="%d/%m/%Y")-as.Date(as.character(hospital_admission_s_WithSampleDate$EchoDate),format="%Y-%m-%d")
hospital_admission_s_WithSampleDate <- hospital_admission_s_WithSampleDate[!is.na(hospital_admission_s_WithSampleDate$Time_Diff),]
hospital_admission_s_WithSampleDate <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$Time_Diff<=0,]

# COPD patients
IcdCode <- COPD_ICD

hospital_admission_patient <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$MAIN_CONDITION %in% IcdCode,]
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_1 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_2 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_3 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_4 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
t <- hospital_admission_s_WithSampleDate[hospital_admission_s_WithSampleDate$OTHER_CONDITION_5 %in% IcdCode,]
hospital_admission_patient <- rbind(hospital_admission_patient,t)

hospital_admission_s_WithSampleDate_COPD <- hospital_admission_patient

COPDPatientList <- unique(hospital_admission_s_WithSampleDate_COPD[,1])

COPD <- union(DrugPatientList, COPDPatientList)
COPD <- data.frame(unique(COPD))
colnames(COPD)[1] <- "PROCHI"
COPD$COPD <- 1
BaselineCohort <- merge(BaselineCohort, COPD, by=c("PROCHI"),all.x = TRUE)
BaselineCohort[is.na(BaselineCohort$COPD),"COPD"] <- 0

######################################
#### Comorbidity Summary #############
######################################
# process the baseline information for prescribing
n <- c("Control", "HFrEF", "HFpEF")
CM <- c("Diabetes","AF","CAD", "COPD", "CKD")
Comorbidity_Summary <- data.frame(CM)
colnames(Comorbidity_Summary)[1] <- "Comorbidity"
Comorbidity_Summary$Control <- 0
Comorbidity_Summary$HFrEF <- 0
Comorbidity_Summary$HFpEF <- 0
Comorbidity_Summary$p <- 0

for (i in 1:3) {
  t <- subset(BaselineCohort, BaselineCohort$FinalDiagnosis == n[i])
  a <- paste0(round(table(t$Diabetes)[2]/dim(t)[1],2)*100,"%") 
  a <- c(a,paste0(round(table(t$AF)[2]/dim(t)[1],2)*100,"%"))
  a <- c(a,paste0(round(table(t$CAD)[2]/dim(t)[1],2)*100,"%"))
  a <- c(a,paste0(round(table(t$COPD)[2]/dim(t)[1],2)*100,"%"))
  a <- c(a,paste0(round(table(t$CKD)[2]/dim(t)[1],2)*100,"%"))
  a <- data.frame(a)
  Comorbidity_Summary[,i+1] <- a
}


t <- chisq.test(BaselineCohort$Diabetes, BaselineCohort$FinalDiagnosis, correct = FALSE)
Comorbidity_Summary[1,"p"] <- t$p.value
t <- chisq.test(BaselineCohort$AF, BaselineCohort$FinalDiagnosis, correct = FALSE)
Comorbidity_Summary[2,"p"] <- t$p.value
t <- chisq.test(BaselineCohort$CAD, BaselineCohort$FinalDiagnosis, correct = FALSE)
Comorbidity_Summary[3,"p"] <- t$p.value
t <- chisq.test(BaselineCohort$COPD, BaselineCohort$FinalDiagnosis, correct = FALSE)
Comorbidity_Summary[4,"p"] <- t$p.value
t <- chisq.test(BaselineCohort$CKD, BaselineCohort$FinalDiagnosis, correct = FALSE)
Comorbidity_Summary[5,"p"] <- t$p.value
write.csv(Comorbidity_Summary, "P:/Project 3634/R/EKOAI/Comorbidity_Summary_11052023.csv")


######################################
######### US2.AI Summary #############
######################################
setwd("P:/Project 3634/Dundee (2)")
getwd()
EKOImageData <- read.csv("MeasurementsDundeeJan3-v3.csv")
### 03032023 data update ##########
setwd("P:/Project 3634/USAI-50")
getwd()
EKOImageData_03032023 <- read.csv("Dundee-2023-02-15-v3 - strain.csv")
EKOImageData <- EKOImageData[,!colnames(EKOImageData) %in% colnames(EKOImageData_03032023)[7:10]]
EKOImageData <- merge(EKOImageData, EKOImageData_03032023[,c(1,6:10)], by=c("PatientID","Visit"))
###################################
BaselineCohort <- merge(BaselineCohort, EKOImageData[,c(3,7:dim(EKOImageData)[2])], by = c("PROCHI"), all.x = TRUE)
# process the baseline information for prescribing
n <- c("Control", "HFpEF", "HFrEF")
CM <- colnames(EKOImageData)[7:dim(EKOImageData)[2]]
US2AIFeatures_Summary <- data.frame(CM)
colnames(US2AIFeatures_Summary)[1] <- "US2AIFeatures"
US2AIFeatures_Summary$Coverage <- 0
US2AIFeatures_Summary$Control <- 0
US2AIFeatures_Summary$HFpEF <- 0
US2AIFeatures_Summary$HFrEF <- 0
US2AIFeatures_Summary$p <- 0
US2AIFeatures_Summary$p_ControlHFpEF <- 0
US2AIFeatures_Summary$p_ControlHFrEF <- 0
US2AIFeatures_Summary$p_HFrEFHFpEF <- 0

for (i in 1:3) {
  t <- subset(BaselineCohort, BaselineCohort$FinalDiagnosis == n[i])
  for (j in 1:length(CM)){
    p <- round(length(na.omit(BaselineCohort[,CM[j]]))/length(BaselineCohort[,CM[j]]),2)
    m <- round(median(as.double(as.character(t[,CM[j]])),na.rm=TRUE),1)
    s <- round(sd(as.double(as.character(t[,CM[j]])),na.rm=TRUE),1)
    US2AIFeatures_Summary[j,"Coverage"] <- p
    US2AIFeatures_Summary[j,i+2] <- paste0(m,"±",s)
    #US2AIFeatures_Summary[j,i+1] <- p
    
  }
}

for (j in c(1:70,72:length(CM))){
  j
  d <- BaselineCohort[,c("FinalDiagnosis", CM[j])]
  d <- d[is.na(d[,2])==FALSE,]
  d[,2] <- as.double(as.character(d[,2])) 
  t <- aov(d[,2] ~ d[,1], data = d)
  US2AIFeatures_Summary[j,"p"] <- summary(t)[[1]][["Pr(>F)"]][1]
}

for (j in c(1:70,72:length(CM))){
  j
  d <- BaselineCohort[,c("FinalDiagnosis", CM[j])]
  d <- d[d$FinalDiagnosis!="HFrEF",]
  d <- d[is.na(d[,2])==FALSE,]
  d[,2] <- as.double(as.character(d[,2])) 
  t <- aov(d[,2] ~ d[,1], data = d)
  US2AIFeatures_Summary[j,"p_ControlHFpEF"] <- summary(t)[[1]][["Pr(>F)"]][1]
}


for (j in c(1:70,72:length(CM))){
  j
  d <- BaselineCohort[,c("FinalDiagnosis", CM[j])]
  d <- d[d$FinalDiagnosis!="HFpEF",]
  d <- d[is.na(d[,2])==FALSE,]
  d[,2] <- as.double(as.character(d[,2])) 
  t <- aov(d[,2] ~ d[,1], data = d)
  US2AIFeatures_Summary[j,"p_ControlHFrEF"] <- summary(t)[[1]][["Pr(>F)"]][1]
}

for (j in c(1:70,72:length(CM))){
  j
  d <- BaselineCohort[,c("FinalDiagnosis", CM[j])]
  d <- d[d$FinalDiagnosis!="Control",]
  d <- d[is.na(d[,2])==FALSE,]
  d[,2] <- as.double(as.character(d[,2])) 
  t <- aov(d[,2] ~ d[,1], data = d)
  US2AIFeatures_Summary[j,"p_HFrEFHFpEF"] <- summary(t)[[1]][["Pr(>F)"]][1]
}

write.csv(US2AIFeatures_Summary, "P:/Project 3634/R/EKOAI/US2AIFeatures_Summary_11052023.csv")
