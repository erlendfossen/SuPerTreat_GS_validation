# Script for analysis of GS4 (pancancer-cisplatin) platin based chemotherapy vs. others vs. no chemo
# endpoint: overall survival (os) and disease-free survival (DFS)


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

# load data ----
path  <- "M:/data_directory/"

# Overall survival data
super_os_name <- "GS4_os_v5.csv"
super_os<- readr::read_delim(paste0(path, super_os_name), delim = ";", escape_double = FALSE, col_names = T, trim_ws = TRUE)

# Disease-free survival data
super_dfs_name <- "GS4_dfs_v5.csv"
super_dfs<- readr::read_delim(paste0(path, super_dfs_name), delim = ";", escape_double = FALSE, col_names = T, trim_ws = TRUE)

# Functions to simplify things: ----

#Function for calculating the gene signature score
gs4_cisplatin <- function(df,
                          unmeasured_genes = c("SLC25A45", "TMC4", "MALAT1")){
    # this function implements a logistic regression model for predicting 
    # drug response for cisplatin based chemotherapy 
    # the model was built using a pan-cancer approach. 
    # data.frame df: gene expression data.frame with genes as columns
    # numeric p_threshold: probability threshold for prediction
    # vector unmeasured_genes: provide a vector of unmeasured genes to be removed
    #                           this is equivalent to imputing 0 for the GE.
    # We miss the following genes in our data: 17 - SLC25A45, 19 - TMC4, 21 - MALAT1
    # change to NULL if all genes are measured
    # reference: 
    # J. D. Wells et al. 2021 DOI:10.1177/11769351211002494
    # "Pan-Cancer Transcriptional Models Predicting Chemosensitivity in Human Tumors"
    
    # betas coefficients from GLM elastic net
    gs4 <- data.frame(Symbol = c('SLFN11', 'TASOR', 'DDX50', 'BCAT1', 
                                 'SALL4', 'PDLIM4', 'METAP2', 
                                 'LRRC8C', 'QKI', 'STOML2', 'SLC35A2', 
                                 'NEBL', 'ABHD11', 'CALCOCO1', 'ZNF91', 
                                 'ZDHHC9', 'SLC25A45', 'MYO5B', 'TMC4', 
                                 'TFF3', 'BSPRY', 'CNPPD1', 'CLDN3', 'BSPRY'),
                      B = c(0.07806286, 0.0584335, 0.055776578, 0.029278636, 
                            0.025196954, 0.021625076, 0.020362816, 0.012095118,
                            0.008967228, 0.001460649, -0.002316639, -0.002369246, 
                            -0.003712791, -0.007393136, -0.009960618, -0.01652152,
                            -0.023694588, -0.034605258, -0.040749427, -0.044766622,
                            -0.056125044, -0.056762982, -0.125547451, -0.146640606))
    
    require(dplyr)
    
    
    # if unmeasured genes are provided as a vector
    # we remove them from the predictor
    # this equivalent to impute 0 to these genes
    if ( is.vector(unmeasured_genes)) {
        for ( gene in unmeasured_genes) {
            gs4 <- gs4[!(gs4["Symbol"]==gene), ]
        }
    }
    
    linear_predictor <- rep(0.0, nrow(df))
    for (model_term in 1:nrow(gs4)) {
        linear_predictor <- linear_predictor +
            df[, gs4[[model_term, "Symbol"]]] * gs4[[model_term, "B"]]
        
    }
    
    results <- linear_predictor[[1]]
    
    return(results)
}


# calculate GS scores & scale ----
# make df with the gene names and coefs
gs4_coefs <- data.frame(Symbol = c('SLFN11', 'TASOR', 'DDX50', 'BCAT1', 
                                          'SALL4', 'PDLIM4', 'METAP2', 
                                          'LRRC8C', 'QKI', 'STOML2', 'SLC35A2', 
                                          'NEBL', 'ABHD11', 'CALCOCO1', 'ZNF91', 
                                          'ZDHHC9', 'SLC25A45', 'MYO5B', 'TMC4', 
                                          'TFF3', 'BSPRY', 'CNPPD1', 'CLDN3', 'BSPRY'),
                               B = c(0.07806286, 0.0584335, 0.055776578, 0.029278636, 
                                     0.025196954, 0.021625076, 0.020362816, 0.012095118,
                                     0.008967228, 0.001460649, -0.002316639, -0.002369246, 
                                     -0.003712791, -0.007393136, -0.009960618, -0.01652152,
                                     -0.023694588, -0.034605258, -0.040749427, -0.044766622,
                                     -0.056125044, -0.056762982, -0.125547451, -0.146640606))

# calculate the signature scores while excluding genes that are not found in the respective datasets
## both are missing 2 genes
GS4_os_scores <- gs4_cisplatin(super_os, unmeasured_genes = gs4_coefs$Symbol[(gs4_coefs$Symbol %in% colnames(super_os))==F])
GS4_dfs_scores <- gs4_cisplatin(super_dfs, unmeasured_genes = gs4_coefs$Symbol[(gs4_coefs$Symbol %in% colnames(super_dfs))==F])

# scale the scores
GS4_os_scores_scaled <- scale(GS4_os_scores)
GS4_dfs_scores_scaled <- scale(GS4_dfs_scores)


# Analysis dataset: clinical variables and GS4 scores ----
super_os_analysis <- super_os %>% select(c(1:19)) # clinical variables
super_os_analysis <- super_os_analysis %>% mutate(GS4_score = GS4_os_scores_scaled ) #add score

super_dfs_analysis <- super_dfs %>% select(c(1:19))
super_dfs_analysis <- super_dfs_analysis %>% mutate(GS4_score = GS4_dfs_scores_scaled)

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

# ctn_stage_7ed_modified 
super_os_analysis <- super_os_analysis %>% mutate(ctn_stage_7ed_modified  = relevel(factor(ctn_stage_7ed_modified), "Stage III"))
super_dfs_analysis <- super_dfs_analysis %>% mutate(ctn_stage_7ed_modified = relevel(factor(ctn_stage_7ed_modified), "Stage III"))

# smoking_category
super_os_analysis <- super_os_analysis %>% mutate(smoking_category = relevel(factor(smoking_category), "Never"))
super_dfs_analysis <- super_dfs_analysis %>% mutate(smoking_category = relevel(factor(smoking_category), "Never"))

# tumor_region_hpv: make HPV+ group the base level since it is most different from the rest
super_os_analysis <- super_os_analysis %>% mutate(tumor_region_hpv = relevel(factor(tumor_region_hpv), "oropharynx_HPVpositive")) 
super_dfs_analysis <- super_dfs_analysis %>% mutate(tumor_region_hpv = relevel(factor(tumor_region_hpv), "oropharynx_HPVpositive")) 

# chemo_platin_agent
super_os_analysis <- super_os_analysis %>% mutate(chemo_platin_agent = 
                                                      factor(chemo_platin_agent, 
                                                             levels = c("Platinum_based", "Non_platinum_based", "No_chemotherapy", "Missing_chemotherapy")))
super_dfs_analysis <- super_dfs_analysis %>% mutate(chemo_platin_agent = 
                                                        factor(chemo_platin_agent, 
                                                               levels = c("Platinum_based", "Non_platinum_based", "No_chemotherapy", "Missing_chemotherapy")))

# chemo in general
super_os_analysis <- super_os_analysis %>% mutate(chemo_chemotherapy_treatment = 
                                                      factor(chemo_chemotherapy_treatment, levels = c("No", "Yes", "Missing")))
super_dfs_analysis <- super_dfs_analysis %>% mutate(chemo_chemotherapy_treatment = 
                                                        factor(chemo_chemotherapy_treatment, levels = c("No", "Yes", "Missing")))

# surgery
super_os_analysis <- super_os_analysis %>% mutate(surge_undergone_cancer_surgery = 
                                                      factor(surge_undergone_cancer_surgery, levels = c("No", "Yes", "Missing")))
super_dfs_analysis <- super_dfs_analysis %>% mutate(surge_undergone_cancer_surgery = 
                                                        factor(surge_undergone_cancer_surgery, levels = c("No", "Yes", "Missing")))


# radiotherapy
super_os_analysis <- super_os_analysis %>% mutate(radio_radiotherapy_treatment = 
                                                      factor(radio_radiotherapy_treatment, levels = c("No", "Yes", "Missing")))
super_dfs_analysis <- super_dfs_analysis %>% mutate(radio_radiotherapy_treatment = 
                                                        factor(radio_radiotherapy_treatment, levels = c("No", "Yes", "Missing")))

# Analysis: gs4_os_2y ----

# full model
gs4_os_2y <- coxph(Surv(osdays_2y,osstatus_2y)~ dataset+ GS4_score*chemo_platin_agent 
                   +tumor_region_hpv+clinical_age_at_diagnosis+
                       surge_undergone_cancer_surgery + radio_radiotherapy_treatment +
                       clinical_sex+ ctn_stage_7ed_modified + smoking_category
                       , data=super_os_analysis, x=TRUE, model=T) 
summary(gs4_os_2y)
survMisc::rsq(gs4_os_2y) #R^2 

# Obtain marginal effects using relevel
## non-platinum
super_os_analysis_nonplatinum <- super_os_analysis %>% mutate(chemo_platin_agent = fct_relevel(chemo_platin_agent, "Non_platinum_based"))
gs4_os_2y_nonplatinum <- coxph(Surv(osdays_2y,osstatus_2y)~ dataset+ GS4_score*chemo_platin_agent 
                               +tumor_region_hpv+clinical_age_at_diagnosis+
                                   surge_undergone_cancer_surgery + radio_radiotherapy_treatment +
                                   clinical_sex+ ctn_stage_7ed_modified + smoking_category
                               , data=super_os_analysis_nonplatinum, x=TRUE, model=T) 
summary(gs4_os_2y_nonplatinum)

## no-chemo
super_os_analysis_nochemo <- super_os_analysis %>% mutate(chemo_platin_agent = fct_relevel(chemo_platin_agent, "No_chemotherapy"))
gs4_os_2y_nochemo <- coxph(Surv(osdays_2y,osstatus_2y)~ dataset+ GS4_score*chemo_platin_agent 
                               +tumor_region_hpv+clinical_age_at_diagnosis+
                               surge_undergone_cancer_surgery + radio_radiotherapy_treatment +
                                   clinical_sex+ ctn_stage_7ed_modified + smoking_category
                               , data=super_os_analysis_nochemo, x=TRUE, model=T) 
summary(gs4_os_2y_nochemo)


# comparable base model without gs4. comparable c index
summary( coxph(Surv(osdays_2y,osstatus_2y)~ dataset+ chemo_platin_agent 
      +tumor_region_hpv+clinical_age_at_diagnosis+
          surge_undergone_cancer_surgery + radio_radiotherapy_treatment +
          clinical_sex+ ctn_stage_7ed_modified + smoking_category
      , data=super_os_analysis, x=TRUE, model=T) )
# c index: 
survMisc::rsq(coxph(Surv(osdays_2y,osstatus_2y)~ dataset+ chemo_platin_agent 
                    +tumor_region_hpv+clinical_age_at_diagnosis+
                      surge_undergone_cancer_surgery + radio_radiotherapy_treatment +
                      clinical_sex+ ctn_stage_7ed_modified + smoking_category
                    , data=super_os_analysis, x=TRUE, model=T))
#r2 

# plot for publication ----
# make main plot
## start with figuring out labels
facet_names_numbers <- c(
  "Platinum_based" = paste0("Platinum-based (N = ", nrow(super_os_analysis[super_os_analysis$chemo_platin_agent=="Platinum_based",]), ")" )
  , "Non_platinum_based" = paste0("Not platinum-based (N = ", nrow(super_os_analysis[super_os_analysis$chemo_platin_agent=="Non_platinum_based",]), ")" )
  , "No_chemotherapy" = paste0("No systemic treatment (N = ", nrow(super_os_analysis[super_os_analysis$chemo_platin_agent=="No_chemotherapy",]), ")" )
  , "Missing_chemotherapy" = paste0("Missing systemic treatment (N = ", nrow(super_os_analysis[super_os_analysis$chemo_platin_agent=="Missing_chemotherapy",]), ")" )
)

facet_names_numbers
# "Platinum-based (N = 362)"       
# "Not platinum-based (N = 51)"      
# "No systemic treatment (N = 155)" 
# "Missing systemic treatment (N = 182)" 


## make new labeller, connecting factor levels with full labels. 
new_labeller <- 
  as_labeller( c(
    "Platinum_based" = "bold((A))~'Platinum-based (N = 362)'"
    , "Non_platinum_based" = "bold((B))~'Not platinum-based (N = 51)'"
    , "No_chemotherapy" = "bold((C))~'No systemic treatment (N = 155)'"
    , "Missing_chemotherapy" = "bold((D))~'Missing systemic treatment (N = 182)'"
  ), label_parsed)

## make plot
gs4_os_2y_plot_v2 <-  plot_surv_area(time="osdays_2y", status="osstatus_2y", variable="GS4_score", title = "Systemic treatment", 
                                  xlab="Time since diagnosis (days)", ylab= "Overall survival probability",legend.title= "Pancancer-cisplatin \nscore",
                                  group="chemo_platin_agent", data=super_os_analysis,model=gs4_os_2y,
                                  facet_args = list(labeller = new_labeller)) + ylim(0,1)


# make tiff file
file_path <- "M:/data_directory/GS4_plot.tiff"
tiff(file_path, units="mm", width=177, height=100, res=300)

gs4_os_2y_plot_v2

dev.off()


# Analysis: gs4_os_5y ----

# full model
gs4_os_5y <- coxph(Surv(osdays_5y,osstatus_5y)~ dataset+ GS4_score*chemo_platin_agent 
                   +tumor_region_hpv+clinical_age_at_diagnosis+
                       surge_undergone_cancer_surgery + radio_radiotherapy_treatment +
                       clinical_sex+ ctn_stage_7ed_modified + smoking_category
                   , data=super_os_analysis, x=TRUE, model=T) 
summary(gs4_os_5y)

# Obtain marginal effects using relevel
## non-platinum
super_os_analysis_nonplatinum <- super_os_analysis %>% mutate(chemo_platin_agent = fct_relevel(chemo_platin_agent, "Non_platinum_based"))
gs4_os_5y_nonplatinum <- coxph(Surv(osdays_5y,osstatus_5y)~ dataset+ GS4_score*chemo_platin_agent 
                               +tumor_region_hpv+clinical_age_at_diagnosis+
                                   surge_undergone_cancer_surgery + radio_radiotherapy_treatment +
                                   clinical_sex+ ctn_stage_7ed_modified + smoking_category
                               , data=super_os_analysis_nonplatinum, x=TRUE, model=T) 
summary(gs4_os_5y_nonplatinum)

## no-chemo
super_os_analysis_nochemo <- super_os_analysis %>% mutate(chemo_platin_agent = fct_relevel(chemo_platin_agent, "No_chemotherapy"))
gs4_os_5y_nochemo <- coxph(Surv(osdays_5y,osstatus_5y)~ dataset+ GS4_score*chemo_platin_agent 
                           +tumor_region_hpv+clinical_age_at_diagnosis+
                               surge_undergone_cancer_surgery + radio_radiotherapy_treatment +
                               clinical_sex+ ctn_stage_7ed_modified + smoking_category
                           , data=super_os_analysis_nochemo, x=TRUE, model=T) 
summary(gs4_os_5y_nochemo)


# Analysis: gs4_dfs_2y ----

# full model
gs4_dfs_2y <- coxph(Surv(dfsdays_2y,dfs_event_2y)~ dataset+ GS4_score*chemo_platin_agent 
                   +tumor_region_hpv+clinical_age_at_diagnosis+
                       surge_undergone_cancer_surgery + radio_radiotherapy_treatment +
                       clinical_sex+ ctn_stage_7ed_modified + smoking_category
                   , data=super_dfs_analysis, x=TRUE, model=T) 
summary(gs4_dfs_2y)
survMisc::rsq(gs4_dfs_2y) # r2 

# Obtain marginal effects using relevel
## non-platinum
super_dfs_analysis_nonplatinum <- super_dfs_analysis %>% mutate(chemo_platin_agent = fct_relevel(chemo_platin_agent, "Non_platinum_based"))
gs4_dfs_2y_nonplatinum <- coxph(Surv(dfsdays_2y,dfs_event_2y)~ dataset+ GS4_score*chemo_platin_agent 
                               +tumor_region_hpv+clinical_age_at_diagnosis+
                                   surge_undergone_cancer_surgery + radio_radiotherapy_treatment +
                                   clinical_sex+ ctn_stage_7ed_modified + smoking_category
                               , data=super_dfs_analysis_nonplatinum, x=TRUE, model=T) 
summary(gs4_dfs_2y_nonplatinum)

## no-chemo
super_dfs_analysis_nochemo <- super_dfs_analysis %>% mutate(chemo_platin_agent = fct_relevel(chemo_platin_agent, "No_chemotherapy"))
gs4_dfs_2y_nochemo <- coxph(Surv(dfsdays_2y,dfs_event_2y)~ dataset+ GS4_score*chemo_platin_agent 
                           +tumor_region_hpv+clinical_age_at_diagnosis+
                               surge_undergone_cancer_surgery + radio_radiotherapy_treatment +
                               clinical_sex+ ctn_stage_7ed_modified + smoking_category
                           , data=super_dfs_analysis_nochemo, x=TRUE, model=T) 
summary(gs4_dfs_2y_nochemo)

# comparable base model without gs4. 
summary( coxph(Surv(dfsdays_2y,dfs_event_2y)~ dataset+ chemo_platin_agent 
               +tumor_region_hpv+clinical_age_at_diagnosis+
                   surge_undergone_cancer_surgery + radio_radiotherapy_treatment +
                   clinical_sex+ ctn_stage_7ed_modified + smoking_category
               , data=super_dfs_analysis, x=TRUE, model=T)  )
# c index: 
survMisc::rsq(coxph(Surv(dfsdays_2y,dfs_event_2y)~ dataset+ chemo_platin_agent 
                    +tumor_region_hpv+clinical_age_at_diagnosis+
                      surge_undergone_cancer_surgery + radio_radiotherapy_treatment +
                      clinical_sex+ ctn_stage_7ed_modified + smoking_category
                    , data=super_dfs_analysis, x=TRUE, model=T))
# R^2 

# Analysis: gs4_dfs_5y ----

# full model
gs4_dfs_5y <- coxph(Surv(dfsdays_5y,dfs_event_5y)~ dataset+ GS4_score*chemo_platin_agent 
                    +tumor_region_hpv+clinical_age_at_diagnosis+
                        surge_undergone_cancer_surgery + radio_radiotherapy_treatment +
                        clinical_sex+ ctn_stage_7ed_modified + smoking_category
                    , data=super_dfs_analysis, x=TRUE, model=T) 
summary(gs4_dfs_5y)


# Obtain marginal effects using relevel
## non-platinum
super_dfs_analysis_nonplatinum <- super_dfs_analysis %>% mutate(chemo_platin_agent = fct_relevel(chemo_platin_agent, "Non_platinum_based"))
gs4_dfs_5y_nonplatinum <- coxph(Surv(dfsdays_5y,dfs_event_5y)~ dataset+ GS4_score*chemo_platin_agent 
                                +tumor_region_hpv+clinical_age_at_diagnosis+
                                    surge_undergone_cancer_surgery + radio_radiotherapy_treatment +
                                    clinical_sex+ ctn_stage_7ed_modified + smoking_category
                                , data=super_dfs_analysis_nonplatinum, x=TRUE, model=T) 
summary(gs4_dfs_5y_nonplatinum)

## no-chemo
super_dfs_analysis_nochemo <- super_dfs_analysis %>% mutate(chemo_platin_agent = fct_relevel(chemo_platin_agent, "No_chemotherapy"))
gs4_dfs_5y_nochemo <- coxph(Surv(dfsdays_5y,dfs_event_5y)~ dataset+ GS4_score*chemo_platin_agent 
                            +tumor_region_hpv+clinical_age_at_diagnosis+
                                surge_undergone_cancer_surgery + radio_radiotherapy_treatment +
                                clinical_sex+ ctn_stage_7ed_modified + smoking_category
                            , data=super_dfs_analysis_nochemo, x=TRUE, model=T) 
summary(gs4_dfs_5y_nochemo)
