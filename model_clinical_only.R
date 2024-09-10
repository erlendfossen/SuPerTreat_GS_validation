# Script for analysis clinical data using all relevant proprietary data and some key public datasets
# 


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
super_os_name <- "clinical_base_all_os_aug2023.csv"
super_os<- readr::read_delim(paste0(path, super_os_name), delim = ";", escape_double = FALSE, col_names = T, trim_ws = TRUE)

# Disease-free survival data
super_dfs_name <- "clinical_base_all_dfs_aug2023.csv"
super_dfs<- readr::read_delim(paste0(path, super_dfs_name), delim = ";", escape_double = FALSE, col_names = T, trim_ws = TRUE)


# new variables needed ----
# os DATA: 0 1 instead of alive dead
super_os <- super_os %>% mutate(CensoredStatus = ifelse(
  follow_status_of_patient=="Alive",yes=0,no=1))

# censor OS at 2 years and 5 years
super_os$osdays_2y<-ifelse(super_os$osprimarydays>2*365.25, yes=2*365.25, no=super_os$osprimarydays)
super_os$osstatus_2y<-ifelse(super_os$osprimarydays>2*365.25, yes=0, no=super_os$CensoredStatus)
super_os$osdays_5y<-ifelse(super_os$osprimarydays>5*365.25, yes=5*365.25, no=super_os$osprimarydays)
super_os$osstatus_5y<-ifelse(super_os$osprimarydays>5*365.25, yes=0, no=super_os$CensoredStatus)

#DFS: 
# censor dfs at 2 years and 5 years
super_dfs$dfsdays_2y<-ifelse(super_dfs$dfsdays>2*365.25, yes=2*365.25, no=super_dfs$dfsdays)
super_dfs$dfs_event_2y<-ifelse(super_dfs$dfsdays>2*365.25, yes=0, no=super_dfs$dfs_event)
super_dfs$dfsdays_5y<-ifelse(super_dfs$dfsdays>5*365.25, yes=5*365.25, no=super_dfs$dfsdays)
super_dfs$dfs_event_5y<-ifelse(super_dfs$dfsdays>5*365.25, yes=0, no=super_dfs$dfs_event)

# Classes per variable ----
## Define default levels etc for categories

# clinical_sex: majority class (male) as first level 
super_os <- super_os %>% mutate(clinical_sex = relevel(factor(clinical_sex), names(which.max((table(clinical_sex)))))) 
super_dfs <- super_dfs %>% mutate(clinical_sex = relevel(factor(clinical_sex), names(which.max((table(clinical_sex))))))

# ctn_disease_extension_diagnosis
super_os <- super_os %>% mutate(ctn_disease_extension_diagnosis = relevel(factor(ctn_disease_extension_diagnosis), "Early disease"))
super_dfs <- super_dfs %>% mutate(ctn_disease_extension_diagnosis = relevel(factor(ctn_disease_extension_diagnosis), "Early disease"))

# smoking_category
super_os <- super_os %>% mutate(smoking_category = relevel(factor(smoking_category), "Never"))
super_dfs <- super_dfs %>% mutate(smoking_category = relevel(factor(smoking_category), "Never"))

# tumor_region_hpv: make HPV+ group the base level since it is most different from the rest
super_os <- super_os %>% mutate(tumor_region_hpv = relevel(factor(tumor_region_hpv), "oropharynx_HPVpositive")) 
super_dfs <- super_dfs %>% mutate(tumor_region_hpv = relevel(factor(tumor_region_hpv), "oropharynx_HPVpositive")) 

# surge_undergone_cancer_surgery
super_os <- super_os %>% mutate(surge_undergone_cancer_surgery = relevel(factor(surge_undergone_cancer_surgery), "No"))
super_dfs <- super_dfs %>% mutate(surge_undergone_cancer_surgery = relevel(factor(surge_undergone_cancer_surgery), "No"))

# radio_radiotherapy_treatment
super_os <- super_os %>% mutate(radio_radiotherapy_treatment = relevel(factor(radio_radiotherapy_treatment), "No"))
super_dfs <- super_dfs %>% mutate(radio_radiotherapy_treatment = relevel(factor(radio_radiotherapy_treatment), "No"))

# chemo_chemotherapy_treatment
super_os <- super_os %>% mutate(chemo_chemotherapy_treatment = relevel(factor(chemo_chemotherapy_treatment), "No"))
super_dfs <- super_dfs %>% mutate(chemo_chemotherapy_treatment = relevel(factor(chemo_chemotherapy_treatment), "No"))

# Analysis: clinical_base_os_2years ----
# Base model with dataset, age (continuous), sex, disease_extension, smoking, tumor_region_hpv, surgery, chemo, radio

# full model
clinical_base_os_2years <- coxph(Surv(osdays_2y,osstatus_2y)~ dataset+tumor_region_hpv+clinical_age_at_diagnosis+
                              clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                              surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
                              chemo_chemotherapy_treatment, data=super_os, x=FALSE, model=TRUE) 
summary(clinical_base_os_2years)

# Analysis: clinical_base_os_5years ----
# Base model with dataset, age (continuous), sex, disease_extension, smoking, tumor_region_hpv, surgery, chemo, radio

# full model
clinical_base_os_5years <- coxph(Surv(osdays_5y,osstatus_5y)~ dataset+tumor_region_hpv+clinical_age_at_diagnosis+
                                   clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                                   surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
                                   chemo_chemotherapy_treatment, data=super_os, x=FALSE, model=TRUE) 
summary(clinical_base_os_5years)

# Analysis: clinical_base_dfs_2years ----
# Base model with dataset, age (continuous), sex, disease_extension, smoking, tumor_region_hpv, surgery, chemo, radio

# full model
clinical_base_dfs_2years <- coxph(Surv(dfsdays_2y,dfs_event_2y)~ dataset+tumor_region_hpv+clinical_age_at_diagnosis+
                                   clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                                   surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
                                   chemo_chemotherapy_treatment, data=super_dfs, x=FALSE, model=TRUE) 
summary(clinical_base_dfs_2years)

# Analysis: clinical_base_dfs_2years ----
# Base model with dataset, age (continuous), sex, disease_extension, smoking, tumor_region_hpv, surgery, chemo, radio

# full model
clinical_base_dfs_5years <- coxph(Surv(dfsdays_5y,dfs_event_5y)~ dataset+tumor_region_hpv+clinical_age_at_diagnosis+
                                    clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                                    surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
                                    chemo_chemotherapy_treatment, data=super_dfs, x=FALSE, model=TRUE) 
summary(clinical_base_dfs_5years)
