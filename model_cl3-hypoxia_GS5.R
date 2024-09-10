# Script for analysis of GS5 (cl3-hypoxia)

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
super_os_name <- "GS5_os_v5.csv"
super_os<- readr::read_delim(paste0(path, super_os_name), delim = ";", escape_double = FALSE, col_names = T, trim_ws = TRUE)

# Disease-free survival data
super_dfs_name <- "GS5_dfs_v5.csv"
super_dfs<- readr::read_delim(paste0(path, super_dfs_name), delim = ";", escape_double = FALSE, col_names = T, trim_ws = TRUE)

# Functions to simplify things: ----

#Function for calculating the gene signature score
gs5_cl3 <- function(df,
                          unmeasured_genes = c( "JAK3",   "CLDN17", "PLPP4",  "PEAR1" , "FIBIN"  )){
  # data.frame df: gene expression data.frame with genes as columns
  # vector unmeasured_genes: provide a vector of unmeasured genes to be removed
  #                           this is equivalent to imputing 0 for the GE.
  # change to NULL if all genes are measured
  # reference: 
  
  # betas coefficients for weighted sum
  gs5 <- data.frame(Symbol = c("CBR3",
                               "COL5A1",
                               "CYP27B1",
                               "DIO2",
                               "DSG1",
                               "ERCC1",
                               "GBE1",
                               "GRIA3",
                               "HNRNPA2B1",
                               "IGFBP5",
                               "JAK3",
                               "KRT2",
                               "PBX1",
                               "PDE7A",
                               "SERPINB8",
                               "SLC22A1",
                               "SLC22A3",
                               "SLPI",
                               "SOX11",
                               "ITPRID2",
                               "TGFB3",
                               "TGM3",
                               "TIMP1",
                               "VEGFC",
                               "NRP2",
                               "SPTLC1",
                               "KLF12",
                               "ANGPTL2",
                               "CLDN17",
                               "RGS17",
                               "HSD17B11",
                               "CYP39A1",
                               "MAP3K20",
                               "RASIP1",
                               "AJAP1",
                               "STK31",
                               "ZC3HAV1",
                               "PGLYRP4",
                               "PELI2",
                               "SLC28A3",
                               "SYNC",
                               "DCUN1D5",
                               "SLC12A8",
                               "ZNF518B",
                               "SLC38A5",
                               "MRGPRX3",
                               "GALNT15",
                               "PLPP4",
                               "PEAR1",
                               "FIBIN"),
                    B = c(0.22613,
                          -0.933,
                          -0.55053,
                          0.62787,
                          -0.06938,
                          0.33007,
                          0.04216,
                          -0.00463,
                          -0.25011,
                          0.03963,
                          0.83392,
                          0.10631,
                          0.3535,
                          -1.03292,
                          0.6833,
                          -0.13804,
                          -0.35736,
                          -0.27494,
                          0.09411,
                          -0.65881,
                          0.10795,
                          0.52079,
                          0.47357,
                          -0.20517,
                          0.77979,
                          0.23361,
                          0.07201,
                          0.67754,
                          0.0211,
                          0.15437,
                          -0.01754,
                          -0.02575,
                          -0.49814,
                          0.15911,
                          -0.03048,
                          -0.28313,
                          -0.00574,
                          -0.53303,
                          0.0073,
                          -0.01494,
                          -0.1753,
                          0.10983,
                          0.01278,
                          -0.00826,
                          -0.00254,
                          0.19623,
                          0.13027,
                          0.15115,
                          -0.30407,
                          0.05492
                    ))
  
  require(dplyr)
 
  # if unmeasured genes are provided as a vector
  # we remove them from the predictor
  # the three missing genes have negative associations with the phenotype
  # this equivalent to impute 0 to these genes
  if ( is.vector(unmeasured_genes)) {
    for ( gene in unmeasured_genes) {
      gs5 <- gs5[!(gs5["Symbol"]==gene), ]
    }
  }
  
  linear_predictor <- rep(0.0, nrow(df))
  for (model_term in 1:nrow(gs5)) {
    linear_predictor <- linear_predictor +
      df[, gs5[[model_term, "Symbol"]]] * gs5[[model_term, "B"]]
    
  } # weighted sum: sum (coef_geneX * geneExpression_geneX)
  
  results <- linear_predictor[[1]]
  
  return(results)
}


# calculate GS scores & scale ----
# make df with the gene names and coefs
gs5_coefs <- data.frame(Symbol = c("CBR3","COL5A1","CYP27B1", "DIO2", "DSG1","ERCC1",
                                          "GBE1", "GRIA3", "HNRNPA2B1","IGFBP5", "JAK3",
                                          "KRT2","PBX1","PDE7A", "SERPINB8","SLC22A1",
                                          "SLC22A3", "SLPI", "SOX11", "ITPRID2",
                                          "TGFB3", "TGM3","TIMP1","VEGFC",
                                          "NRP2",   "SPTLC1",
                                          "KLF12",     "ANGPTL2",
                                          "CLDN17",    "RGS17",
                                          "HSD17B11", "CYP39A1",
                                          "MAP3K20","RASIP1",
                                          "AJAP1","STK31",
                                          "ZC3HAV1",
                                          "PGLYRP4",
                                          "PELI2",
                                          "SLC28A3",
                                          "SYNC",
                                          "DCUN1D5",
                                          "SLC12A8",
                                          "ZNF518B",
                                          "SLC38A5",
                                          "MRGPRX3",
                                          "GALNT15",
                                          "PLPP4",
                                          "PEAR1",
                                          "FIBIN"),
                               B = c(0.22613,
                                     -0.933,
                                     -0.55053,
                                     0.62787,
                                     -0.06938,
                                     0.33007,
                                     0.04216,
                                     -0.00463,
                                     -0.25011,
                                     0.03963,
                                     0.83392,
                                     0.10631,
                                     0.3535,
                                     -1.03292,
                                     0.6833,
                                     -0.13804,
                                     -0.35736,
                                     -0.27494,
                                     0.09411,
                                     -0.65881,
                                     0.10795,
                                     0.52079,
                                     0.47357,
                                     -0.20517,
                                     0.77979,
                                     0.23361,
                                     0.07201,
                                     0.67754,
                                     0.0211,
                                     0.15437,
                                     -0.01754,
                                     -0.02575,
                                     -0.49814,
                                     0.15911,
                                     -0.03048,
                                     -0.28313,
                                     -0.00574,
                                     -0.53303,
                                     0.0073,
                                     -0.01494,
                                     -0.1753,
                                     0.10983,
                                     0.01278,
                                     -0.00826,
                                     -0.00254,
                                     0.19623,
                                     0.13027,
                                     0.15115,
                                     -0.30407,
                                     0.05492
                               ))

# calculate the signature scores while excluding genes that are not found in the respective datasets
## both are missing 3 genes: "PLPP4" "PEAR1" "FIBIN"
GS5_os_scores <- gs5_cl3(super_os, unmeasured_genes = gs5_coefs$Symbol[(gs5_coefs$Symbol %in% colnames(super_os))==F])
GS5_dfs_scores <- gs5_cl3(super_dfs, unmeasured_genes = gs5_coefs$Symbol[(gs5_coefs$Symbol %in% colnames(super_dfs))==F])

# scale the scores
GS5_os_scores_scaled <- scale(GS5_os_scores)
GS5_dfs_scores_scaled <- scale(GS5_dfs_scores)



# Analysis dataset: clinical variables and GS5 scores ----
super_os_analysis <- super_os %>% select(c(1:19)) # clinical variables
super_os_analysis <- super_os_analysis %>% mutate(GS5_score = GS5_os_scores_scaled ) #add score

super_dfs_analysis <- super_dfs %>% select(c(1:19))
super_dfs_analysis <- super_dfs_analysis %>% mutate(GS5_score = GS5_dfs_scores_scaled)

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

# chemo_cetuximab_agent
super_os_analysis <- super_os_analysis %>% mutate(chemo_cetuximab_agent = 
                                                    factor(chemo_cetuximab_agent, 
                                                           levels = c("Cetuximab_based", "Non_cetuximab_based", "No_chemotherapy", "Missing_chemotherapy")))
super_dfs_analysis <- super_dfs_analysis %>% mutate(chemo_cetuximab_agent = 
                                                      factor(chemo_cetuximab_agent, 
                                                             levels = c("Cetuximab_based", "Non_cetuximab_based", "No_chemotherapy", "Missing_chemotherapy")))

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


# analysis: gs5_os_2y ----

# full model
gs5_os_2y <- coxph(Surv(osdays_2y,osstatus_2y)~ dataset+ GS5_score*chemo_cetuximab_agent
                   +tumor_region_hpv+clinical_age_at_diagnosis+
                     radio_radiotherapy_treatment + surge_undergone_cancer_surgery + 
                     clinical_sex+ ctn_stage_7ed_modified + smoking_category
                   , data=super_os_analysis, x=TRUE, model=T) 
summary(gs5_os_2y)
survMisc::rsq(gs5_os_2y) # R2 

# Obtain marginal effects using relevel
## non-cetuximab
super_os_analysis_noncetuximab <- super_os_analysis %>% mutate(chemo_cetuximab_agent = fct_relevel(chemo_cetuximab_agent, "Non_cetuximab_based"))
gs5_os_2y_noncetuximab <- coxph(Surv(osdays_2y,osstatus_2y)~ dataset+ GS5_score*chemo_cetuximab_agent 
                               +tumor_region_hpv+clinical_age_at_diagnosis+
                                 radio_radiotherapy_treatment + surge_undergone_cancer_surgery + 
                                 clinical_sex+ ctn_stage_7ed_modified + smoking_category
                               , data=super_os_analysis_noncetuximab, x=TRUE, model=T) 
summary(gs5_os_2y_noncetuximab)

## no-chemo
super_os_analysis_nochemo <- super_os_analysis %>% mutate(chemo_cetuximab_agent = fct_relevel(chemo_cetuximab_agent, "No_chemotherapy"))
gs5_os_2y_nochemo <- coxph(Surv(osdays_2y,osstatus_2y)~ dataset+ GS5_score*chemo_cetuximab_agent 
                           +tumor_region_hpv+clinical_age_at_diagnosis+
                             radio_radiotherapy_treatment + surge_undergone_cancer_surgery + 
                             clinical_sex+ ctn_stage_7ed_modified + smoking_category
                           , data=super_os_analysis_nochemo, x=TRUE, model=T) 
summary(gs5_os_2y_nochemo)


# comparable c index by removing gs5 and having a base model
summary( coxph(Surv(osdays_2y,osstatus_2y)~ dataset+ chemo_cetuximab_agent
      +tumor_region_hpv+clinical_age_at_diagnosis+
        radio_radiotherapy_treatment + surge_undergone_cancer_surgery + 
        clinical_sex+ ctn_stage_7ed_modified + smoking_category
      , data=super_os_analysis, x=TRUE, model=T) )
# c index: 
survMisc::rsq(coxph(Surv(osdays_2y,osstatus_2y)~ dataset+ chemo_cetuximab_agent
                    +tumor_region_hpv+clinical_age_at_diagnosis+
                      radio_radiotherapy_treatment + surge_undergone_cancer_surgery + 
                      clinical_sex+ ctn_stage_7ed_modified + smoking_category
                    , data=super_os_analysis, x=TRUE, model=T))
# R2 

# plot for publication ----

# make main plot
## start with figuring out labels
facet_names_numbers <- c(
  "Cetuximab_based" = paste0("Cetuximab-based (N = ", nrow(super_os_analysis[super_os_analysis$chemo_cetuximab_agent=="Cetuximab_based",]), ")" )
  , "Non_cetuximab_based" = paste0("Not cetuximab-based (N = ", nrow(super_os_analysis[super_os_analysis$chemo_cetuximab_agent=="Non_cetuximab_based",]), ")" )
  , "No_chemotherapy" = paste0("No systemic treatment (N = ", nrow(super_os_analysis[super_os_analysis$chemo_cetuximab_agent=="No_chemotherapy",]), ")" )
  , "Missing_chemotherapy" = paste0("Missing systemic treatment (N = ", nrow(super_os_analysis[super_os_analysis$chemo_cetuximab_agent=="Missing_chemotherapy",]), ")" )
)

facet_names_numbers
# "Cetuximab-based (N = 35)"        
# "Not cetuximab-based (N = 378)"      
# "No systemic treatment (N = 155)" 
# "Missing systemic treatment (N = 182)"


## make new labeller, connecting factor levels with full labels. 
new_labeller <- 
  as_labeller( c(
    "Cetuximab_based" = "bold((A))~'Cetuximab-based (N = 35)'"
    , "Non_cetuximab_based" = "bold((B))~'Not cetuximab-based (N = 378)'"
    , "No_chemotherapy" = "bold((C))~'No systemic treatment (N = 155)'"
    , "Missing_chemotherapy" = "bold((D))~'Missing systemic treatment (N = 182)'"
  ), label_parsed)

## make plot
gs5_os_2y_plot_v2 <-  
  plot_surv_area(time="osdays_2y", status="osstatus_2y", variable="GS5_score", title = "Systemic treatment", 
                 xlab="Time since diagnosis (days)", ylab= "Overall survival probability",legend.title= "Cl3-hypoxia \nscore",
                 group="chemo_cetuximab_agent", data=super_os_analysis,model=gs5_os_2y,
                 facet_args = list(labeller = new_labeller)) + ylim(0,1)


# make tiff file
file_path <- "M:/data_directory/GS5_plot.tiff"
tiff(file_path, units="mm", width=165, height=100, res=300)

gs5_os_2y_plot_v2

dev.off()



# analysis: gs5_os_5y ----

# full model
gs5_os_5y <- coxph(Surv(osdays_5y,osstatus_5y)~ dataset+ GS5_score*chemo_cetuximab_agent
                   +tumor_region_hpv+clinical_age_at_diagnosis+
                     radio_radiotherapy_treatment + surge_undergone_cancer_surgery +
                     clinical_sex+ ctn_stage_7ed_modified + smoking_category
                   , data=super_os_analysis, x=TRUE, model=T) 
summary(gs5_os_5y)

# Obtain marginal effects using relevel
## non-cetuximab
super_os_analysis_noncetuximab <- super_os_analysis %>% mutate(chemo_cetuximab_agent = fct_relevel(chemo_cetuximab_agent, "Non_cetuximab_based"))
gs5_os_5y_noncetuximab <- coxph(Surv(osdays_5y,osstatus_5y)~ dataset+ GS5_score*chemo_cetuximab_agent 
                                +tumor_region_hpv+clinical_age_at_diagnosis+
                                  radio_radiotherapy_treatment + surge_undergone_cancer_surgery +
                                  clinical_sex+ ctn_stage_7ed_modified + smoking_category
                                , data=super_os_analysis_noncetuximab, x=TRUE, model=T) 
summary(gs5_os_5y_noncetuximab)

## no-chemo
super_os_analysis_nochemo <- super_os_analysis %>% mutate(chemo_cetuximab_agent = fct_relevel(chemo_cetuximab_agent, "No_chemotherapy"))
gs5_os_5y_nochemo <- coxph(Surv(osdays_5y,osstatus_5y)~ dataset+ GS5_score*chemo_cetuximab_agent 
                           +tumor_region_hpv+clinical_age_at_diagnosis+
                             radio_radiotherapy_treatment + surge_undergone_cancer_surgery +
                             clinical_sex+ ctn_stage_7ed_modified + smoking_category
                           , data=super_os_analysis_nochemo, x=TRUE, model=T) 
summary(gs5_os_5y_nochemo)

# analysis: gs5_dfs_2y ----

# full model
gs5_dfs_2y <- coxph(Surv(dfsdays_2y,dfs_event_2y)~ dataset+ GS5_score*chemo_cetuximab_agent
                   +tumor_region_hpv+clinical_age_at_diagnosis+
                     radio_radiotherapy_treatment + surge_undergone_cancer_surgery +
                     clinical_sex+ ctn_stage_7ed_modified + smoking_category
                   , data=super_dfs_analysis, x=TRUE, model=T) 
summary(gs5_dfs_2y)
survMisc::rsq(gs5_dfs_2y) # r2 

# Obtain marginal effects using relevel
## non-cetuximab
super_dfs_analysis_noncetuximab <- super_dfs_analysis %>% mutate(chemo_cetuximab_agent = fct_relevel(chemo_cetuximab_agent, "Non_cetuximab_based"))
gs5_dfs_2y_noncetuximab <- coxph(Surv(dfsdays_2y,dfs_event_2y)~ dataset+ GS5_score*chemo_cetuximab_agent 
                                +tumor_region_hpv+clinical_age_at_diagnosis+
                                  radio_radiotherapy_treatment + surge_undergone_cancer_surgery +
                                  clinical_sex+ ctn_stage_7ed_modified + smoking_category
                                , data=super_dfs_analysis_noncetuximab, x=TRUE, model=T) 
summary(gs5_dfs_2y_noncetuximab)

## no-chemo
super_dfs_analysis_nochemo <- super_dfs_analysis %>% mutate(chemo_cetuximab_agent = fct_relevel(chemo_cetuximab_agent, "No_chemotherapy"))
gs5_dfs_2y_nochemo <- coxph(Surv(dfsdays_2y,dfs_event_2y)~ dataset+ GS5_score*chemo_cetuximab_agent 
                           +tumor_region_hpv+clinical_age_at_diagnosis+
                             radio_radiotherapy_treatment + surge_undergone_cancer_surgery +
                             clinical_sex+ ctn_stage_7ed_modified + smoking_category
                           , data=super_dfs_analysis_nochemo, x=TRUE, model=T) 
summary(gs5_dfs_2y_nochemo)


# comparable c index by removing gs5 and having a base model
summary( coxph(Surv(dfsdays_2y,dfs_event_2y)~ dataset+ chemo_cetuximab_agent
               +tumor_region_hpv+clinical_age_at_diagnosis+
                 radio_radiotherapy_treatment + surge_undergone_cancer_surgery +
                 clinical_sex+ ctn_stage_7ed_modified + smoking_category
               , data=super_dfs_analysis, x=TRUE, model=T) )
# c index: 
survMisc::rsq(coxph(Surv(dfsdays_2y,dfs_event_2y)~ dataset+ chemo_cetuximab_agent
                    +tumor_region_hpv+clinical_age_at_diagnosis+
                      radio_radiotherapy_treatment + surge_undergone_cancer_surgery +
                      clinical_sex+ ctn_stage_7ed_modified + smoking_category
                    , data=super_dfs_analysis, x=TRUE, model=T))
# r2 

# analysis: gs5_dfs_5y ----

# full model
gs5_dfs_5y <- coxph(Surv(dfsdays_5y,dfs_event_5y)~ dataset+ GS5_score*chemo_cetuximab_agent
                    +tumor_region_hpv+clinical_age_at_diagnosis+
                      radio_radiotherapy_treatment + surge_undergone_cancer_surgery +
                      clinical_sex+ ctn_stage_7ed_modified + smoking_category
                    , data=super_dfs_analysis, x=TRUE, model=T) 
summary(gs5_dfs_5y)

# Obtain marginal effects using relevel
## non-cetuximab
super_dfs_analysis_noncetuximab <- super_dfs_analysis %>% mutate(chemo_cetuximab_agent = fct_relevel(chemo_cetuximab_agent, "Non_cetuximab_based"))
gs5_dfs_5y_noncetuximab <- coxph(Surv(dfsdays_5y,dfs_event_5y)~ dataset+ GS5_score*chemo_cetuximab_agent 
                                 +tumor_region_hpv+clinical_age_at_diagnosis+
                                   radio_radiotherapy_treatment + surge_undergone_cancer_surgery +
                                   clinical_sex+ ctn_stage_7ed_modified + smoking_category
                                 , data=super_dfs_analysis_noncetuximab, x=TRUE, model=T) 
summary(gs5_dfs_5y_noncetuximab)

## no-chemo
super_dfs_analysis_nochemo <- super_dfs_analysis %>% mutate(chemo_cetuximab_agent = fct_relevel(chemo_cetuximab_agent, "No_chemotherapy"))
gs5_dfs_5y_nochemo <- coxph(Surv(dfsdays_5y,dfs_event_5y)~ dataset+ GS5_score*chemo_cetuximab_agent 
                            +tumor_region_hpv+clinical_age_at_diagnosis+
                              radio_radiotherapy_treatment + surge_undergone_cancer_surgery +
                              clinical_sex+ ctn_stage_7ed_modified + smoking_category
                            , data=super_dfs_analysis_nochemo, x=TRUE, model=T )
summary(gs5_dfs_5y_nochemo)
