# Script for analysis GS3 RSI, prediction of radiosensitivity
# endpoint: overall survival (os) and disease-free survival (DFS)
# HNSCC patients receiving RT or not

# libraries ----
library(readr)
library(dplyr)
library(tidyselect)
library(tidyverse)
library(hacksig)
library(survival)
library(ggplot2)
library(ggpubr)
library(survminer)
library(lubridate)
library(sm)
library(tidycmprsk)
library(ggsurvfit)
library(contsurvplot)
library(survRM2)
library(forestplot)

#
# load data ----
path  <- "M:/data_directory/"

# Overall survival data
super_os_name <- "GS3_os_v3.csv"
super_os<- readr::read_delim(paste0(path, super_os_name), delim = ";", escape_double = FALSE, col_names = T, trim_ws = TRUE)

# Disease-free survival data
super_dfs_name <- "GS3_dfs_v3.csv"
super_dfs<- readr::read_delim(paste0(path, super_dfs_name), delim = ";", escape_double = FALSE, col_names = T, trim_ws = TRUE)

# Functions to simplify things: ----

#Function for calculating the gene signature score
RSI_function <- function(data){
    # data is the dataset you want the score for
    
    ## based on Dai et al 2021 npj Genomic Medicine & eschrich et al 2009
    ## Score is calculated by ranking the 10 genes and giving the highest expressed gene a score of 10 and the lowest a score of 1
    ## Unclear from the text how to rank it, but to get different scores per patient, need to rank within patients, not across patients.
    ## Then using the ranks in the equation
    
    # make a function for calculating the score based on coefficients and ranked gene expression
    RSI_score_function <- function(AR, JUN, STAT1, PRKCB, RELA, ABL1, SUMO1, CDK1, HDAC1, IRF1){
        RSI_score <- (-0.0098009 * AR) + (0.0128283 * JUN) +
            (0.0254552 * STAT1) + (-0.0017589 * PRKCB) + (-0.0038171 * RELA) +
            (0.1070213 * ABL1) + (-0.0002509 * SUMO1) +
            (-0.0092431 * CDK1) +  (-0.0204469 * HDAC1) + (-0.0441683 * IRF1)
        return(RSI_score)
    }
    
    # Select the 10 genes and keep only these per patient
    RSI_df <- data %>% select(supertreat_id, AR, JUN, STAT1, PRKCB, RELA,
                                  ABL1, SUMO1, CDK1, HDAC1, IRF1)
    
    # For each patient, sort genes by rank, assign 10 to highest rank, and calculate the RSI score
    RSI_df$RSI_score <- NA
    for (i in 1:length(RSI_df$supertreat_id)){
        # sort genes from low to high
        rsi_sorted <- sort(as.vector(RSI_df[i, 2:11]), decreasing=F)
        rsi_sorted_rank_value <- data.frame(1,2,3,4,5,6,7,8,9,10)
        colnames(rsi_sorted_rank_value) <- colnames(rsi_sorted)
        RSI_df$RSI_score[i] <- RSI_score_function(AR= rsi_sorted_rank_value$AR, JUN=rsi_sorted_rank_value$JUN, 
                                            STAT1=rsi_sorted_rank_value$STAT1,
                                            PRKCB = rsi_sorted_rank_value$PRKCB, 
                                            RELA=rsi_sorted_rank_value$RELA, ABL1 = rsi_sorted_rank_value$ABL1,
                                            SUMO1 = rsi_sorted_rank_value$SUMO1, 
                                            CDK1 = rsi_sorted_rank_value$CDK1, HDAC1 = rsi_sorted_rank_value$HDAC1,
                                            IRF1 = rsi_sorted_rank_value$IRF1)
    }
    return(RSI_df)
}

# calculate GS scores and scale ----

# Scores based on custom RSI function
RSI_scores_os <- RSI_function(super_os)
RSI_scores_dfs <- RSI_function(super_dfs)

#scale RSI scores
RSI_scores_os$RSI_score_scaled <- scale(RSI_scores_os$RSI_score)
RSI_scores_dfs$RSI_score_scaled <- scale(RSI_scores_dfs$RSI_score)

# Analysis dataset: clinical variables and RSI scores ----
super_os_analysis <- super_os %>% select(c(1:19)) # clinical variables
super_os_analysis <- super_os_analysis %>% mutate(RSI_score = RSI_scores_os$RSI_score_scaled) #add score

super_dfs_analysis <- super_dfs %>% select(c(1:20))
super_dfs_analysis <- super_dfs_analysis %>% mutate(RSI_score = RSI_scores_dfs$RSI_score_scaled)

# new variables needed ----
# os DATA: 0 1 instead of alive dead
super_os_analysis <- super_os_analysis %>% mutate(CensoredStatus = ifelse(
    follow_status_of_patient=="Alive",yes=0,no=1))

# censor OS at 2 years and 5 years
super_os_analysis$osdays_2y<-ifelse(super_os_analysis$osprimarydays>2*365.25, yes=2*365.25, no=super_os_analysis$osprimarydays)
super_os_analysis$osstatus_2y<-ifelse(super_os_analysis$osprimarydays>2*365.25, yes=0, no=super_os_analysis$CensoredStatus)
super_os_analysis$osdays_5y<-ifelse(super_os_analysis$osprimarydays>5*365.25, yes=5*365.25, no=super_os_analysis$osprimarydays)
super_os_analysis$osstatus_5y<-ifelse(super_os_analysis$osprimarydays>5*365.25, yes=0, no=super_os_analysis$CensoredStatus)

#DFS: 
# censor dfs at 2 years and 5 years
super_dfs_analysis$dfsdays_2y<-ifelse(super_dfs_analysis$dfsdays>2*365.25, yes=2*365.25, no=super_dfs_analysis$dfsdays)
super_dfs_analysis$dfs_event_2y<-ifelse(super_dfs_analysis$dfsdays>2*365.25, yes=0, no=super_dfs_analysis$dfs_event)
super_dfs_analysis$dfsdays_5y<-ifelse(super_dfs_analysis$dfsdays>5*365.25, yes=5*365.25, no=super_dfs_analysis$dfsdays)
super_dfs_analysis$dfs_event_5y<-ifelse(super_dfs_analysis$dfsdays>5*365.25, yes=0, no=super_dfs_analysis$dfs_event)

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

# tumor_region_hpv: make HPV+ group the base level since it is most different from the rest
super_os_analysis <- super_os_analysis %>% mutate(tumor_region_hpv = relevel(factor(tumor_region_hpv), "oropharynx_HPVpositive")) 
super_dfs_analysis <- super_dfs_analysis %>% mutate(tumor_region_hpv = relevel(factor(tumor_region_hpv), "oropharynx_HPVpositive")) 

# surge_undergone_cancer_surgery
super_os_analysis <- super_os_analysis %>% mutate(surge_undergone_cancer_surgery = relevel(factor(surge_undergone_cancer_surgery), "No"))
super_dfs_analysis <- super_dfs_analysis %>% mutate(surge_undergone_cancer_surgery = relevel(factor(surge_undergone_cancer_surgery), "No"))

# radio_radiotherapy_treatment
super_os_analysis <- super_os_analysis %>% mutate(radio_radiotherapy_treatment = 
                                                      factor(radio_radiotherapy_treatment, levels= c("No", "Yes", "Missing")))
super_dfs_analysis <- super_dfs_analysis %>% mutate(radio_radiotherapy_treatment = 
                                                        factor(radio_radiotherapy_treatment, levels= c("No", "Yes", "Missing")))

# chemo_chemotherapy_treatment
super_os_analysis <- super_os_analysis %>% mutate(chemo_chemotherapy_treatment = relevel(factor(chemo_chemotherapy_treatment), "No"))
super_dfs_analysis <- super_dfs_analysis %>% mutate(chemo_chemotherapy_treatment = relevel(factor(chemo_chemotherapy_treatment), "No"))

# Analysis: gs3_os_2y ----

# full model
gs3_os_2y <- coxph(Surv(osdays_2y,osstatus_2y)~ dataset+ RSI_score*radio_radiotherapy_treatment 
                   +tumor_region_hpv+clinical_age_at_diagnosis+
                       clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                       surge_undergone_cancer_surgery+
                       chemo_chemotherapy_treatment, data=super_os_analysis, x=TRUE, model=T) 

summary(gs3_os_2y)
survMisc::rsq(gs3_os_2y) # R^2 

# Obtain marginal effect using relevel
super_os_analysis2 <- super_os_analysis %>% mutate(radio_radiotherapy_treatment = fct_relevel(radio_radiotherapy_treatment, "Yes"))
gs3_os_2y_relevel <- coxph(Surv(osdays_2y,osstatus_2y)~ dataset+ RSI_score*radio_radiotherapy_treatment 
                   +tumor_region_hpv+clinical_age_at_diagnosis+
                       clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                       surge_undergone_cancer_surgery+
                       chemo_chemotherapy_treatment, data=super_os_analysis2, x=TRUE, model=T) 
summary(gs3_os_2y_relevel)

# equivalent base model where RSI is removed and can compare c index
summary( coxph(Surv(osdays_2y,osstatus_2y)~ dataset+ radio_radiotherapy_treatment 
      +tumor_region_hpv+clinical_age_at_diagnosis+
        clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
        surge_undergone_cancer_surgery+
        chemo_chemotherapy_treatment, data=super_os_analysis, x=TRUE, model=T) )
# c index: 
survMisc::rsq(coxph(Surv(osdays_2y,osstatus_2y)~ dataset+ radio_radiotherapy_treatment 
                    +tumor_region_hpv+clinical_age_at_diagnosis+
                      clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                      surge_undergone_cancer_surgery+
                      chemo_chemotherapy_treatment, data=super_os_analysis, x=TRUE, model=T))
# R2 

# plot for publication ----
# make main plot
## start with figuring out labels
facet_names_numbers <- c(
  "No" = paste0(" No (N = ", nrow(super_os_analysis[super_os_analysis$radio_radiotherapy_treatment=="No",]), ")" )
  , "Yes" = paste0(" Yes (N = ", nrow(super_os_analysis[super_os_analysis$radio_radiotherapy_treatment=="Yes",]), ")" )
  , "Missing" = paste0(" Missing (N = ", nrow(super_os_analysis[super_os_analysis$radio_radiotherapy_treatment=="Missing",]), ")" )
)

facet_names_numbers
# " No (N = 208)" 
# " Yes (N = 752)" 
# " Missing (N = 137)"


## make new labeller, connecting factor levels with full labels. 
new_labeller <- 
  as_labeller( c(
    "No" = "bold((A))~' No (N = 208)'"
    , "Yes" = "bold((B))~' Yes (N = 752)'"
    , "Missing" = "bold((C))~'Missing (N = 137)'"
  ), label_parsed)

## make plot
gs3_os_2y_plot_v2 <-  plot_surv_area(time="osdays_2y", status="osstatus_2y", variable="RSI_score", title = "Radiotherapy", 
                                  xlab="Time since diagnosis (days)", ylab= "Overall survival probability",legend.title= "RSI score",
                                  group="radio_radiotherapy_treatment", data=super_os_analysis,model=gs3_os_2y,
                                  facet_args = list(labeller = new_labeller)) + ylim(0,1)



# make tiff file
file_path <- "M:/data_directory/GS3_plot.tiff"
tiff(file_path, units="mm", width=140, height=70, res=300)

gs3_os_2y_plot_v2

dev.off()






# Analysis: gs3_os_5y ----

# full model
gs3_os_5y <- coxph(Surv(osdays_5y,osstatus_5y)~ dataset+ RSI_score*radio_radiotherapy_treatment 
                   +tumor_region_hpv+clinical_age_at_diagnosis+
                       clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                       surge_undergone_cancer_surgery+
                       chemo_chemotherapy_treatment, data=super_os_analysis, x=TRUE, model=T) 

summary(gs3_os_5y)


# Obtain marginal effect using relevel
super_os_analysis2 <- super_os_analysis %>% mutate(radio_radiotherapy_treatment = fct_relevel(radio_radiotherapy_treatment, "Yes"))
gs3_os_5y_relevel <- coxph(Surv(osdays_5y,osstatus_5y)~ dataset+ RSI_score*radio_radiotherapy_treatment 
                           +tumor_region_hpv+clinical_age_at_diagnosis+
                               clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                               surge_undergone_cancer_surgery+
                               chemo_chemotherapy_treatment, data=super_os_analysis2, x=TRUE, model=T) 
summary(gs3_os_5y_relevel)


# Analysis: gs3_dfs_2y ----

# full model
gs3_dfs_2y <- coxph(Surv(dfsdays_2y,dfs_event_2y)~ dataset+ RSI_score*radio_radiotherapy_treatment 
                   +tumor_region_hpv+clinical_age_at_diagnosis+
                       clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                       surge_undergone_cancer_surgery+
                       chemo_chemotherapy_treatment, data=super_dfs_analysis, x=TRUE, model=T) 

summary(gs3_dfs_2y)
survMisc::rsq(gs3_dfs_2y) # R2 

# Obtain marginal effect using relevel
super_dfs_analysis2 <- super_dfs_analysis %>% mutate(radio_radiotherapy_treatment = fct_relevel(radio_radiotherapy_treatment, "Yes"))
gs3_dfs_2y_relevel <- coxph(Surv(dfsdays_2y,dfs_event_2y)~ dataset+ RSI_score*radio_radiotherapy_treatment 
                           +tumor_region_hpv+clinical_age_at_diagnosis+
                               clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                               surge_undergone_cancer_surgery+
                               chemo_chemotherapy_treatment, data=super_dfs_analysis2, x=TRUE, model=T) 
summary(gs3_dfs_2y_relevel)

# equivalent base model where RSI is removed and can compare c index
summary( coxph(Surv(dfsdays_2y,dfs_event_2y)~ dataset+ radio_radiotherapy_treatment 
               +tumor_region_hpv+clinical_age_at_diagnosis+
                 clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                 surge_undergone_cancer_surgery+
                 chemo_chemotherapy_treatment, data=super_dfs_analysis, x=TRUE, model=T) )
# c index: 

survMisc::rsq(coxph(Surv(dfsdays_2y,dfs_event_2y)~ dataset+ radio_radiotherapy_treatment 
                    +tumor_region_hpv+clinical_age_at_diagnosis+
                      clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                      surge_undergone_cancer_surgery+
                      chemo_chemotherapy_treatment, data=super_dfs_analysis, x=TRUE, model=T))
# r2 

# Analysis: gs3_dfs_5y ----

# full model
gs3_dfs_5y <- coxph(Surv(dfsdays_5y,dfs_event_5y)~ dataset+ RSI_score*radio_radiotherapy_treatment 
                    +tumor_region_hpv+clinical_age_at_diagnosis+
                        clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                        surge_undergone_cancer_surgery+
                        chemo_chemotherapy_treatment, data=super_dfs_analysis, x=TRUE, model=T) 

summary(gs3_dfs_5y)

# Obtain marginal effect using relevel
super_dfs_analysis2 <- super_dfs_analysis %>% mutate(radio_radiotherapy_treatment = fct_relevel(radio_radiotherapy_treatment, "Yes"))
gs3_dfs_5y_relevel <- coxph(Surv(dfsdays_5y,dfs_event_5y)~ dataset+ RSI_score*radio_radiotherapy_treatment 
                            +tumor_region_hpv+clinical_age_at_diagnosis+
                                clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                                surge_undergone_cancer_surgery+
                                chemo_chemotherapy_treatment, data=super_dfs_analysis2, x=TRUE, model=T) 
summary(gs3_dfs_5y_relevel)

