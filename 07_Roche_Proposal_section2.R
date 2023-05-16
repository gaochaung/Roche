# 18/02/2022
####################################################
#### update BaselineCohort Final Diagnosis##########
####################################################
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

# Biomarker BarcodetoProchi mapping
BiomarkerWithProchi <- merge(Biomarker, BarcodetoProchi, by.x=c("Barcode.ID"), by.y = ("Barcode"))
BiomarkerWithProchi_Selected <- subset(BiomarkerWithProchi, select=c(9,11,13,16))
BaselineCohort <- merge(BaselineCohort, BiomarkerWithProchi_Selected, by = c("PROCHI"), all.x = TRUE)
# RocheBiomarker BarcodetoProchi mapping
RocheBiomarkerWithProchi <- merge(RocheBiomarker, BarcodetoProchi, by.x=c("Barcode.ID"), by.y = ("Barcode"))
RocheBiomarkerWithProchi_Selected <- subset(RocheBiomarkerWithProchi, select=c(16,8:13))
BaselineCohort <- merge(BaselineCohort, RocheBiomarkerWithProchi_Selected, by = c("PROCHI"), all.x = TRUE)

# demography
BaselineCohort<- merge(Demography[,c(1,3,5,7)], BaselineCohort, by = c("PROCHI"))
BaselineCohort[BaselineCohort$FinalDiagnosis=="HFmrEF","FinalDiagnosis"] <- "HFrEF"

#####################################################
##################### table 2.1  ####################
#####################################################
#death  ####
Death_Selected <- merge(Death, DeathCAUSE, by=c("PROCHI"), all = TRUE)
Death_Selected<- merge(Death_Selected[,c(1,2,31:34)], BaselineCohort[,c(1,7)], by = c("PROCHI"))
#Death_Selected$Time_Diff <- as.Date(as.character(Death_Selected$date_of_death),format="%d/%m/%Y")-as.Date(as.character(Death_Selected$EchoDate),format="%Y-%m-%d")
#Death_Selected <- subset(Death_Selected, Death_Selected$Time_Diff<=3650)
Death_Selected$Death <- 1
BaselineCohort <- merge(BaselineCohort, Death_Selected[,c("PROCHI","Death","date_of_death")], by=c("PROCHI"),all.x = TRUE)
BaselineCohort[is.na(BaselineCohort$Death),"Death"] <- 0
#hosptial admission #####
hospital_admission_Selected <-merge(hospital_admission[,c(1,8,17:22)], BaselineCohort[,c("PROCHI","EchoDate")], by=c("PROCHI"))
hospital_admission_Selected$Time_Diff <- as.Date(as.character(hospital_admission_Selected$ADMISSION_DATE),format="%d/%m/%Y")-as.Date(as.character(hospital_admission_Selected$EchoDate),format="%Y-%m-%d")
hospital_admission_Selected <- hospital_admission_Selected[!is.na(hospital_admission_Selected$Time_Diff),]
hospital_admission_Selected <- hospital_admission_Selected[hospital_admission_Selected$Time_Diff>=0,]
hospital_admission_patient <- hospital_admission_Selected[hospital_admission_Selected$MAIN_CONDITION %like% "^I50",]
t <- hospital_admission_Selected[hospital_admission_Selected$OTHER_CONDITION_1 %like% "^I50",]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
t <- hospital_admission_Selected[hospital_admission_Selected$OTHER_CONDITION_2 %like% "^I50",]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
t <- hospital_admission_Selected[hospital_admission_Selected$OTHER_CONDITION_3 %like% "^I50",]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
t <- hospital_admission_Selected[hospital_admission_Selected$OTHER_CONDITION_4 %like% "^I50",]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
t <- hospital_admission_Selected[hospital_admission_Selected$OTHER_CONDITION_5 %like% "^I50",]
hospital_admission_patient <- rbind(hospital_admission_patient,t)
hospital_admission_patient <- hospital_admission_patient %>% group_by(PROCHI) %>% top_n(-1, Time_Diff)
hospital_admission_patient <- unique(hospital_admission_patient[,c(1,2)])
hospital_admission_patient$Hospital_Admission <- 1
BaselineCohort <- merge(BaselineCohort, hospital_admission_patient, by=c("PROCHI"),all.x = TRUE)
BaselineCohort[is.na(BaselineCohort$Hospital_Admission),"Hospital_Admission"] <- 0

n <- c("Control", "HFrEF", "HFpEF")
TableResult <- data.frame("Event")
colnames(TableResult)[1] <- ""
TableResult$Control <- 0
TableResult$HFrEF <- 0
TableResult$HFpEF <- 0
for (i in 1:3) {
  t <- subset(BaselineCohort, BaselineCohort$FinalDiagnosis == n[i])
  t_s <- subset(t,t$Death==1|t$Hospital_Admission==1)
  TableResult[1,i+1] <- paste0(dim(t_s)[1],"/",dim(t)[1],", ",round(dim(t_s)[1]/dim(t)[1],2)*100,"%")
  
}
write.csv(TableResult, "P:/Project 3634/R/EKOAI/12052023_table2.1.csv")
#####################################################
##################### table 2.2  ####################
#####################################################
BaselineCohort <- BaselineCohort %>% mutate(Survivor = case_when(
  BaselineCohort$Death == 0 & BaselineCohort$Hospital_Admission==0 ~ 1,
  TRUE ~ 0
))
# Bimarker Summary #######
CM <- c(colnames(Biomarker)[c(9,11,13)],colnames(RocheBiomarker)[c(8:13)])
TableResult <- data.frame(CM)
colnames(TableResult)[1] <- "Biomarker"
TableResult$Survivors <- 0
TableResult$None_Survivors <- 0
TableResult$AUC <- 0
TableResult$p_value <- 0
for (j in 1:length(CM)){
  t <- BaselineCohort[BaselineCohort$FinalDiagnosis!="Control",c(CM[j],"Survivor")]
  t_survivor <- t[t$Survivor==1,]
  t_nonesurvivor <- t[t$Survivor==0,]
  q <- quantile(as.double(as.character(t_survivor[,CM[j]])), probs = c(.25,.5,.75), na.rm = TRUE)
  c <- dim(t_survivor)[1]
  TableResult[j,2] <- paste0(round(q[2],1)," [",round(q[1],1),"-",round(q[3],1),"], ",c)
  q <- quantile(as.double(as.character(t_nonesurvivor[,CM[j]])), probs = c(.25,.5,.75), na.rm = TRUE)
  c <- dim(t_nonesurvivor)[1]
  TableResult[j,3] <- paste0(round(q[2],1)," [",round(q[1],1),"-",round(q[3],1),"], ",c)
  p <- as.double(as.character(t[,CM[j]]))
  ROC_curve <- pROC::roc(factor(t$Survivor),log10(p))
  A <- round(pROC::auc(ROC_curve),2)
  CI <- ci.auc(ROC_curve)
  cp <- cutpointr(log10(p),t$Survivor,method = maximize_metric, metric = sum_sens_spec, na.rm = TRUE)
  ppvnpv <- coords(ROC_curve,cp$optimal_cutpoint,"threshold",ret = c("ppv","npv"))
  cp_optimal=round(10^cp$optimal_cutpoint,2)
  TableResult[j,4] <- paste0(A,"; (",round(CI[1],2),"-",round(CI[3],2),"); ",cp_optimal,"; ",round(ppvnpv$ppv,2),"; ",round(ppvnpv$npv,2))
  TableResult[j,5] <- wilcox.test(as.double(as.character(t_survivor[,1])),as.double(as.character(t_nonesurvivor[,1])))$p.value
}
write.csv(TableResult, "P:/Project 3634/R/EKOAI/12052023_table2.2.csv")


write.csv(TableResult, "P:/Project 3634/R/EKOAI/table2.2.1_update.csv")

# 2.2.2  #######
CM <- c(colnames(Biomarker)[c(9,11,13)],colnames(RocheBiomarker)[c(8:13)])
TableResult <- data.frame(CM)
colnames(TableResult)[1] <- "Biomarker"
TableResult$Survivors <- 0
TableResult$None_Survivors <- 0
TableResult$AUC <- 0
TableResult$p_value <- 0
for (j in 1:length(CM)){
  t <- BaselineCohort[BaselineCohort$FinalDiagnosis=="HFpEF",c(CM[j],"Survivor")]
  t_survivor <- t[t$Survivor==1,]
  t_nonesurvivor <- t[t$Survivor==0,]
  q <- quantile(as.double(as.character(t_survivor[,CM[j]])), probs = c(.25,.5,.75), na.rm = TRUE)
  c <- dim(t_survivor)[1]
  TableResult[j,2] <- paste0(round(q[2],1)," [",round(q[1],1),"-",round(q[3],1),"], ",c)
  q <- quantile(as.double(as.character(t_nonesurvivor[,CM[j]])), probs = c(.25,.5,.75), na.rm = TRUE)
  c <- dim(t_nonesurvivor)[1]
  TableResult[j,3] <- paste0(round(q[2],1)," [",round(q[1],1),"-",round(q[3],1),"], ",c)
  p <- as.double(as.character(t[,CM[j]]))
  ROC_curve <- pROC::roc(factor(t$Survivor),log10(p))
  A <- round(pROC::auc(ROC_curve),2)
  CI <- ci.auc(ROC_curve)
  cp <- cutpointr(log10(p),t$Survivor,method = maximize_metric, metric = sum_sens_spec, na.rm = TRUE)
  ppvnpv <- coords(ROC_curve,cp$optimal_cutpoint,"threshold",ret = c("ppv","npv"))
  cp_optimal=round(10^cp$optimal_cutpoint,2)
  TableResult[j,4] <- paste0(A,"; (",round(CI[1],2),"-",round(CI[3],2),"); ",cp_optimal,"; ",round(ppvnpv$ppv,2),"; ",round(ppvnpv$npv,2))
  TableResult[j,5] <- wilcox.test(as.double(as.character(t_survivor[,1])),as.double(as.character(t_nonesurvivor[,1])))$p.value
}
write.csv(TableResult, "P:/Project 3634/R/EKOAI/table2.2.2_update.csv")


#####################################################
##################### table 2.3  ####################
#####################################################
BaselineCohort <- BaselineCohort %>% mutate(SurvivalTime = case_when(
  BaselineCohort$Hospital_Admission==1 ~ as.Date(as.character(BaselineCohort$ADMISSION_DATE),format="%d/%m/%Y")-as.Date(as.character(BaselineCohort$EchoDate),format="%Y-%m-%d"),
  BaselineCohort$Hospital_Admission==0 & BaselineCohort$Death==1  ~ as.Date(as.character(BaselineCohort$date_of_death),format="%d/%m/%Y")-as.Date(as.character(BaselineCohort$EchoDate),format="%Y-%m-%d")-1,
  BaselineCohort$Hospital_Admission==0 & BaselineCohort$Death==0 ~ as.Date(as.character("15/07/2021"),format="%d/%m/%Y")-as.Date(as.character(BaselineCohort$EchoDate),format="%Y-%m-%d")
))
BaselineCohort <- BaselineCohort %>% mutate(Censor = case_when(
  BaselineCohort$Survivor==1 ~ 0,
  BaselineCohort$Survivor==0 ~ 1
))

BaselineCohort$surv_object <- Surv(time=BaselineCohort$SurvivalTime, event=BaselineCohort$Censor)
###### 2.3 #######
CM <- c(colnames(Biomarker)[c(9,11,13)],colnames(RocheBiomarker)[c(8:13)])
TableResult <- data.frame(CM)
colnames(TableResult)[1] <- "Biomarker"
TableResult$HR <- 0
for (j in 1:length(CM)){
  t <- BaselineCohort[BaselineCohort$FinalDiagnosis!="Control",c(CM[j],"surv_object")]
  t[,1] <- as.double(as.character(t[,1]))
  t$logbiomarker <- log2(t[,1])
  fit.coxph <- coxph(t$surv_object ~ t$logbiomarker, data = t)
  fit.coxphsummary <- summary(fit.coxph)
  HR <- round(fit.coxphsummary$coef[,2],2)
  CI <- fit.coxphsummary$conf.int[,c(3:4)]
  TableResult[j,2] <- paste0(HR," [",round(CI[1],2),"-",round(CI[2],2),"] ")
}
write.csv(TableResult, "P:/Project 3634/R/EKOAI/table2.3.csv")
###### 2.3.1 ######
CM <- c(colnames(Biomarker)[c(9,11,13)],colnames(RocheBiomarker)[c(8:13)])
TableResult <- data.frame(CM)
colnames(TableResult)[1] <- "Biomarker"
TableResult$HR <- 0
for (j in 1:length(CM)){
  t <- BaselineCohort[BaselineCohort$FinalDiagnosis=="HFrEF",c(CM[j],"surv_object")]
  t[,1] <- as.double(as.character(t[,1]))
  t$logbiomarker <- log2(t[,1])
  fit.coxph <- coxph(t$surv_object ~ t$logbiomarker, data = t)
  fit.coxphsummary <- summary(fit.coxph)
  HR <- round(fit.coxphsummary$coef[,2],2)
  CI <- fit.coxphsummary$conf.int[,c(3:4)]
  TableResult[j,2] <- paste0(HR," [",round(CI[1],2),"-",round(CI[2],2),"] ")
}
write.csv(TableResult, "P:/Project 3634/R/EKOAI/table2.3.1.csv")
###### 2.3.2 #######
CM <- c(colnames(Biomarker)[c(9,11,13)],colnames(RocheBiomarker)[c(8:13)])
TableResult <- data.frame(CM)
colnames(TableResult)[1] <- "Biomarker"
TableResult$HR <- 0
for (j in 1:length(CM)){
  t <- BaselineCohort[BaselineCohort$FinalDiagnosis=="HFpEF",c(CM[j],"surv_object")]
  t[,1] <- as.double(as.character(t[,1]))
  t$logbiomarker <- log2(t[,1])
  fit.coxph <- coxph(t$surv_object ~ t$logbiomarker, data = t)
  fit.coxphsummary <- summary(fit.coxph)
  HR <- round(fit.coxphsummary$coef[,2],2)
  CI <- fit.coxphsummary$conf.int[,c(3:4)]
  TableResult[j,2] <- paste0(HR," [",round(CI[1],2),"-",round(CI[2],2),"] ")
}
write.csv(TableResult, "P:/Project 3634/R/EKOAI/table2.3.2.csv")