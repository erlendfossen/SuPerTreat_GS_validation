# Script for analysis GS1 (172-GS)


# libraries ----
library(readr)
library(dplyr)
library(tidyselect)
library(tidyverse)
library(hacksig)
library(survival)
library(survminer)
library(lubridate)
library(sm)
library(tidycmprsk)
library(ggsurvfit)
library(contsurvplot)
library(gtsummary)
library(flextable)
library(officer)
library(forestplot)
library(survRM2)



# load data ----
path  <- "M:/data_directory/"

# Overall survival data
super_os_name <- "GS1_GS2_all_OS.csv"
super_os<- readr::read_delim(paste0(path, super_os_name), delim = ";", escape_double = FALSE, col_names = T, trim_ws = TRUE)

# Disease-free survival data
super_dfs_name <- "GS1_GS2_all_DFS.csv"
super_dfs<- readr::read_delim(paste0(path, super_dfs_name), delim = ";", escape_double = FALSE, col_names = T, trim_ws = TRUE)

# calculate GS scores and z-scale the score ----
## Hacksig package requires a normalized gene expression matrix with symbols as row names and samples as columns

# Subset to only include omics data, then transpose and make into a matrix
super_os_df <- data.frame(super_os) # needed to give rownames
rownames(super_os_df) <- super_os_df$supertreat_id
super_os_omics_matrix <- super_os_df[,-c(1:20)] %>% data.matrix() %>% t() #1:20 are the non-omics columns

super_dfs_df <- data.frame(super_dfs) # needed to give rownames
rownames(super_dfs_df) <- super_dfs_df$supertreat_id
super_dfs_omics_matrix <- super_dfs_df[,-c(1:20)] %>% data.matrix() %>% t() #1:20 are the non-omics columns

# check how many genes are present in GS1
check_sig(super_os_omics_matrix, signatures="dececco2014_int172") #76 % present
check_sig(super_dfs_omics_matrix, signatures="dececco2014_int172") #80 % present

# obtain scores per sample
GS1_scores_os<- hack_sig(super_os_omics_matrix, signatures="dececco2014_int172", method="original") 
GS1_scores_dfs<- hack_sig(super_dfs_omics_matrix, signatures="dececco2014_int172", method="original") 

# scale scores: z-score with mean 0 and SD 1
GS1_scores_os$GS1_score_scaled <- scale(GS1_scores_os$dececco2014_int172)
GS1_scores_dfs$GS1_score_scaled <- scale(GS1_scores_dfs$dececco2014_int172)

# Analysis dataset: clinical variables and GS scores ----
# rename id column and signature & remove redundant column
GS1_scores_os <- GS1_scores_os %>% rename(supertreat_id = sample_id, GS1_score = GS1_score_scaled) %>%
  select(-dececco2014_int172)
GS1_scores_dfs <- GS1_scores_dfs %>% rename(supertreat_id = sample_id, GS1_score = GS1_score_scaled) %>%
  select(-dececco2014_int172)

# select only clinical variables
super_os_analysis <- super_os %>% select(c(1:20)) # clinical variables
super_dfs_analysis <- super_dfs %>% select(c(1:20))

# merge GS1_score with super_os_analysis /super_dfs_analysis
super_os_analysis <- merge(super_os_analysis, GS1_scores_os, by="supertreat_id") 
super_dfs_analysis <- merge(super_dfs_analysis, GS1_scores_dfs, by="supertreat_id") 

# new variables needed ----
# os DATA: 0 1 instead of alive dead
super_os_analysis <- super_os_analysis %>% mutate(CensoredStatus = ifelse(
  follow_status_of_patient=="Alive",yes=0,no=1))

# median follow-up time
## Using reverse kaplan-meier method (set events to 0 and censored to 1)
fit_followup <- survfit(Surv(osprimarydays, CensoredStatus==0) ~ 1, data = super_os_analysis)
surv_median(fit_followup) # 1325 days
surv_median(fit_followup)$median/365 # 3.63 years


# censor OS at 2 years and 5 years
super_os_analysis$osdays_2y<-ifelse(super_os_analysis$osprimarydays>2*365.25, yes=2*365.25, no=super_os_analysis$osprimarydays)
super_os_analysis$osstatus_2y<-ifelse(super_os_analysis$osprimarydays>2*365.25, yes=0, no=super_os_analysis$CensoredStatus)
super_os_analysis$osdays_5y<-ifelse(super_os_analysis$osprimarydays>5*365.25, yes=5*365.25, no=super_os_analysis$osprimarydays)
super_os_analysis$osstatus_5y<-ifelse(super_os_analysis$osprimarydays>5*365.25, yes=0, no=super_os_analysis$CensoredStatus)

# make new variable for HPV status and remake tumor region to not split by hpv
super_os_analysis <- super_os_analysis %>% mutate(HPV_status = if_else(
  tumor_region_hpv=="oropharynx_HPVpositive",true="Positive",false="Negative")) %>% 
  mutate(HPV_status = relevel(factor(HPV_status), "Positive"))# pos vs neg hpv status
super_os_analysis <- super_os_analysis %>% mutate(tumor_region = if_else(
  grepl("oropharynx_HPV",tumor_region_hpv),true="oropharynx",false=tumor_region_hpv)) %>% 
  mutate(tumor_region = relevel(factor(tumor_region), "oropharynx")) # dont separate oropharynx 

#DFS: 
# censor dfs at 2 years and 5 years
super_dfs_analysis$dfsdays_2y<-ifelse(super_dfs_analysis$dfsdays>2*365.25, yes=2*365.25, no=super_dfs_analysis$dfsdays)
super_dfs_analysis$dfs_event_2y<-ifelse(super_dfs_analysis$dfsdays>2*365.25, yes=0, no=super_dfs_analysis$dfs_event)
super_dfs_analysis$dfsdays_5y<-ifelse(super_dfs_analysis$dfsdays>5*365.25, yes=5*365.25, no=super_dfs_analysis$dfsdays)
super_dfs_analysis$dfs_event_5y<-ifelse(super_dfs_analysis$dfsdays>5*365.25, yes=0, no=super_dfs_analysis$dfs_event)

# make new variable for HPV status and remake tumor region to not split by hpv
super_dfs_analysis <- super_dfs_analysis %>% mutate(HPV_status = if_else(
  tumor_region_hpv=="oropharynx_HPVpositive",true="Positive",false="Negative")) %>% 
  mutate(HPV_status = relevel(factor(HPV_status), "Positive"))# pos vs neg hpv status
super_dfs_analysis <- super_dfs_analysis %>% mutate(tumor_region = if_else(
  grepl("oropharynx_HPV",tumor_region_hpv),true="oropharynx",false=tumor_region_hpv)) %>% 
  mutate(tumor_region = relevel(factor(tumor_region), "oropharynx")) # dont separate oropharynx 


# Classes per variable ----
## Define default levels etc for categories

# clinical_sex: majority class as first level 
super_os_analysis <- super_os_analysis %>% mutate(clinical_sex = relevel(factor(clinical_sex), names(which.max((table(clinical_sex)))))) 
super_dfs_analysis <- super_dfs_analysis %>% mutate(clinical_sex = relevel(factor(clinical_sex), names(which.max((table(clinical_sex))))))

# ctn_disease_extension_diagnosis
super_os_analysis <- super_os_analysis %>% mutate(ctn_disease_extension_diagnosis = relevel(factor(ctn_disease_extension_diagnosis), "Early disease"))
super_dfs_analysis <- super_dfs_analysis %>% mutate(ctn_disease_extension_diagnosis = relevel(factor(ctn_disease_extension_diagnosis), "Early disease"))

# smoking_category
super_os_analysis <- super_os_analysis %>% mutate(smoking_category = relevel(factor(smoking_category), "Never"))
super_dfs_analysis <- super_dfs_analysis %>% mutate(smoking_category = relevel(factor(smoking_category), "Never"))

# tumor_region: already set to oropharynx

# HPV_status: already set to positive

# surge_undergone_cancer_surgery
super_os_analysis <- super_os_analysis %>% mutate(surge_undergone_cancer_surgery = relevel(factor(surge_undergone_cancer_surgery), "No"))
super_dfs_analysis <- super_dfs_analysis %>% mutate(surge_undergone_cancer_surgery = relevel(factor(surge_undergone_cancer_surgery), "No"))

# radio_radiotherapy_treatment
super_os_analysis <- super_os_analysis %>% mutate(radio_radiotherapy_treatment = relevel(factor(radio_radiotherapy_treatment), "No"))
super_dfs_analysis <- super_dfs_analysis %>% mutate(radio_radiotherapy_treatment = relevel(factor(radio_radiotherapy_treatment), "No"))

# chemo_chemotherapy_treatment
super_os_analysis <- super_os_analysis %>% mutate(chemo_chemotherapy_treatment = relevel(factor(chemo_chemotherapy_treatment), "No"))
super_dfs_analysis <- super_dfs_analysis %>% mutate(chemo_chemotherapy_treatment = relevel(factor(chemo_chemotherapy_treatment), "No"))

# Analysis: gs1_os_2y ----

# full model
gs1_os_2y <- coxph(Surv(osdays_2y,osstatus_2y)~ dataset+ GS1_score*HPV_status+tumor_region+clinical_age_at_diagnosis+
                                        clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                                        surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
                                        chemo_chemotherapy_treatment, data=super_os_analysis, x=TRUE, model=T) 
summary(gs1_os_2y)
survMisc::rsq(gs1_os_2y) # r^2 

# Obtain marginal effect using relevel
super_os_analysis2 <- super_os_analysis %>% mutate(HPV_status = fct_relevel(HPV_status, "Negative"))
gs1_os_2y_relevel <- coxph(Surv(osdays_2y,osstatus_2y)~ dataset+ GS1_score*HPV_status+tumor_region+clinical_age_at_diagnosis+
                     clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                     surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
                     chemo_chemotherapy_treatment, data=super_os_analysis2, x=TRUE, model=T) 
summary(gs1_os_2y_relevel)

# equivalent clinical base model: removing GS1 to obtain new and directly comparable C-index
summary(coxph(Surv(osdays_2y,osstatus_2y)~ dataset+ HPV_status+tumor_region+clinical_age_at_diagnosis+
                clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
                chemo_chemotherapy_treatment, data=super_os_analysis, x=TRUE, model=T) )
## Concordance
survMisc::rsq(coxph(Surv(osdays_2y,osstatus_2y)~ dataset+ HPV_status+tumor_region+clinical_age_at_diagnosis+
                      clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                      surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
                      chemo_chemotherapy_treatment, data=super_os_analysis, x=TRUE, model=T))
# R^2 


# main plot:  ----
# make main plot
## start with figuring out labels
facet_names_numbers <- c(
  "Positive" = paste0(" Positive (N = ", nrow(super_os_analysis[super_os_analysis$HPV_status=="Positive",]), ")" )
  , "Negative" = paste0(" Negative (N = ", nrow(super_os_analysis[super_os_analysis$HPV_status=="Negative",]), ")" )
)
facet_names_numbers
# " Positive (N = 156)" 
# " Negative (N = 941)" 

## make new labeller, connecting factor levels with full labels. 
new_labeller <- 
  as_labeller( c(
    "Positive" = "bold((A))~' Positive (N = 156)'"
    , "Negative" = "bold((B))~ 'Negative (N = 941)'"
  ), label_parsed)

## make plot
gs1_os_2y_plot_v2 <- plot_surv_area(time="osdays_2y", status="osstatus_2y", variable="GS1_score", title = "HPV status", 
               xlab="Time since diagnosis (days)", ylab= "Overall survival probability",legend.title= "172-GS score",
               group="HPV_status", data=super_os_analysis,model=gs1_os_2y,
               facet_args = list(labeller = new_labeller)) + ylim(0,1)


# Add 2 histograms

gg_hist_pos<- super_os_analysis %>% filter(HPV_status=="Positive") %>% ggplot(aes(GS1_score))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.tag = element_text()) +
  geom_histogram( color="#e9ecef", alpha=0.7, position = 'identity')+
  #scale_fill_manual(values=c("#69b3a2", "#404080")) +
  labs(title="HPV-positive",  tag = "(C)", y = "Count")+ xlim(c(-5, 3.7)) + 
  xlab("172-GS score") 

gg_hist_neg<- super_os_analysis %>% filter(HPV_status=="Negative") %>% ggplot(aes(GS1_score))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.tag = element_text()) +
  geom_histogram( color="#e9ecef", alpha=0.7, position = 'identity')+
  #scale_fill_manual(values=c("#69b3a2", "#404080")) +
  labs(title="HPV-negative", tag = "(D)", y = "Count")+ xlim(c(-5, 3.7)) + 
  xlab("172-GS score") 

# combine into one plot

gridExtra::grid.arrange(gs1_os_2y_plot_v2, gg_hist_pos, gg_hist_neg, 
                        layout_matrix = rbind(c(1, 1), c(2,3)), heights=c(4,1.5))


# make tiff file
file_path <- "M:/data_directory/plot_GS1.tiff"
tiff(file_path, units="mm", width=120, height=150, res=300)

gridExtra::grid.arrange(gs1_os_2y_plot_v2, gg_hist_pos, gg_hist_neg, 
                        layout_matrix = rbind(c(1, 1), c(2,3)), heights=c(4,1.5))

dev.off()



# Analysis: gs1_os_5y ----
# full model
gs1_os_5y <- coxph(Surv(osdays_5y,osstatus_5y)~ dataset+ GS1_score*HPV_status+tumor_region+clinical_age_at_diagnosis+
                     clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                     surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
                     chemo_chemotherapy_treatment, data=super_os_analysis, x=TRUE, model=T) 
summary(gs1_os_5y)


# Obtain marginal effect using relevel
super_os_analysis2 <- super_os_analysis %>% mutate(HPV_status = fct_relevel(HPV_status, "Negative"))
gs1_os_5y_relevel <- coxph(Surv(osdays_5y,osstatus_5y)~ dataset+ GS1_score*HPV_status+tumor_region+clinical_age_at_diagnosis+
                             clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                             surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
                             chemo_chemotherapy_treatment, data=super_os_analysis2, x=TRUE, model=T) 
summary(gs1_os_5y_relevel)

# Analysis: gs1_dfs_2y ----

# full model
gs1_dfs_2y <- coxph(Surv(dfsdays_2y,dfs_event_2y)~ dataset+ GS1_score*HPV_status+tumor_region+clinical_age_at_diagnosis+
                     clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                     surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
                     chemo_chemotherapy_treatment, data=super_dfs_analysis, x=TRUE, model=T) 
summary(gs1_dfs_2y)
survMisc::rsq(gs1_dfs_2y) # R^2 

# Obtain marginal effect using relevel
super_dfs_analysis2 <- super_dfs_analysis %>% mutate(HPV_status = fct_relevel(HPV_status, "Negative"))
gs1_dfs_2y_relevel <- coxph(Surv(dfsdays_2y,dfs_event_2y)~ dataset+ GS1_score*HPV_status+tumor_region+clinical_age_at_diagnosis+
                              clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                              surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
                              chemo_chemotherapy_treatment, data=super_dfs_analysis2, x=TRUE, model=T)
summary(gs1_dfs_2y_relevel)

# equivalent clinical base model: removing GS1 to obtain new and directly comparable C-index
summary(coxph(Surv(dfsdays_2y,dfs_event_2y)~ dataset+ HPV_status+tumor_region+clinical_age_at_diagnosis+
                clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
                chemo_chemotherapy_treatment, data=super_dfs_analysis, x=TRUE, model=T)  )
## Concordance
survMisc::rsq(coxph(Surv(dfsdays_2y,dfs_event_2y)~ dataset+ HPV_status+tumor_region+clinical_age_at_diagnosis+
                      clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                      surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
                      chemo_chemotherapy_treatment, data=super_dfs_analysis, x=TRUE, model=T) )
# R^2 

# Analysis: gs1_dfs_5y ----

# full model
gs1_dfs_5y <- coxph(Surv(dfsdays_5y,dfs_event_5y)~ dataset+ GS1_score*HPV_status+tumor_region+clinical_age_at_diagnosis+
                      clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                      surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
                      chemo_chemotherapy_treatment, data=super_dfs_analysis, x=TRUE, model=T) 
summary(gs1_dfs_5y)


# Obtain marginal effect using relevel
super_dfs_analysis2 <- super_dfs_analysis %>% mutate(HPV_status = fct_relevel(HPV_status, "Negative"))
gs1_dfs_5y_relevel <- coxph(Surv(dfsdays_5y,dfs_event_5y)~ dataset+ GS1_score*HPV_status+tumor_region+clinical_age_at_diagnosis+
                              clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                              surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
                              chemo_chemotherapy_treatment, data=super_dfs_analysis2, x=TRUE, model=T)
summary(gs1_dfs_5y_relevel)

