# Script for analysis of GS2 (3-cluster HPV) 


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


# Functions to simplify things: ----

#Function for calculating the gene signature score
gs2_hpv <- function(df, unmeasured_genes = c("FCER2",       "IKBKG" ,      "IGHVII-44-2", "BPIFA3" ,     "ZNF300P1"  ,  "DNAJB8"    ,
                                               "FLCN" ,       "LINC00927" ,  "EPHA1-AS1" ,  "LINC00226"  , "FKBP9P1" ,    "BEND4" ,     
                                               "MIR30E" ,     "TRIM51CP" ,   "GAPLINC"  )){
  # data.frame df: gene expression data.frame with genes as columns
  # vector unmeasured_genes: provide a vector of unmeasured genes to be removed
  #                           this is equivalent to imputing 0 for the GE.
  # We miss the following genes in our data: 17 - SLC25A45, 19 - TMC4, 21 - MALAT1
  # change to NULL if all genes are measured

  # betas coefficients for weighted sum
  

  
  gs2 <- data.frame(Symbol = c("COL6A1","FCER2", "GPR4", "PTK7","S100A1", "SNAI2",
                               "SPARC",  "IKBKG","ANGPTL2", "SEC61G","IGHVII-44-2",    "EFEMP2",
                               "CMC2", "RINT1", "LENG1","SGIP1", "BPIFA3", "ZNF300P1", "DNAJB8", "FLCN",
                               "LINC00927",  "EPHA1-AS1", "LINC00226", "FKBP9P1", "BEND4",  "MIR30E",
                               "TRIM51CP", "GAPLINC"),
                    B = c(-0.2803, 0.04178,  -0.92826,  0.00923, -1.27554,
                          -1.08628,   -0.07866,  -0.08676, -0.40695,
                          -0.16735,  0.02409, -0.37154,-0.13505,
                          0.16931,   -1.50627, -0.87416,   -0.52965,
                          -1.03186,  -0.51708, 0.90364, -0.22078,
                          -0.606, 1.62887, -1.29993, 0.89997,
                          0.31241,  -1.8549,  -0.67773))
  
  require(dplyr)
  # gs2$Symbol %in% colnames(super_os) #  find which genes are missing in our dataset
  # sum(gs2$Symbol %in% colnames(super_os)) # only 13 of 28 genes are found
  # sum(gs2$Symbol %in% colnames(super_dfs)) # only 16 of 28 genes are found
  # gs2$Symbol[(gs2$Symbol %in% colnames(super_os))==F] #gives names of missing genes for OS
  # gs2$Symbol[(gs2$Symbol %in% colnames(super_dfs))==F] #gives names of missing genes for OS
  
  # if unmeasured genes are provided as a vector
  # we remove them from the predictor
  # this equivalent to impute 0 to these genes
  if ( is.vector(unmeasured_genes)) {
    for ( gene in unmeasured_genes) {
      gs2 <- gs2[!(gs2["Symbol"]==gene), ]
    }
  }
  
  linear_predictor <- rep(0.0, nrow(df))
  for (model_term in 1:nrow(gs2)) {
    linear_predictor <- linear_predictor +
      df[, gs2[[model_term, "Symbol"]]] * gs2[[model_term, "B"]]
    
  } # weighted sum: sum (coef_geneX * geneExpression_geneX)
  
  results <- linear_predictor[[1]]
  
  return(results)
}

# calculate GS scores & scale ----
# make df with the gene names and coefs
gs2_coefs <- data.frame(Symbol = c("COL6A1","FCER2", "GPR4", "PTK7","S100A1", "SNAI2",
                             "SPARC",  "IKBKG","ANGPTL2", "SEC61G","IGHVII-44-2",    "EFEMP2",
                             "CMC2", "RINT1", "LENG1","SGIP1", "BPIFA3", "ZNF300P1", "DNAJB8", "FLCN",
                             "LINC00927",  "EPHA1-AS1", "LINC00226", "FKBP9P1", "BEND4",  "MIR30E",
                             "TRIM51CP", "GAPLINC"),
                  B = c(-0.2803, 0.04178,  -0.92826,  0.00923, -1.27554,
                        -1.08628,   -0.07866,  -0.08676, -0.40695,
                        -0.16735,  0.02409, -0.37154,-0.13505,
                        0.16931,   -1.50627, -0.87416,   -0.52965,
                        -1.03186,  -0.51708, 0.90364, -0.22078,
                        -0.606, 1.62887, -1.29993, 0.89997,
                        0.31241,  -1.8549,  -0.67773))

# calculate the signature scores while excluding genes that are not found in the respective datasets
GS2_os_scores <- gs2_hpv(super_os, unmeasured_genes = gs2_coefs$Symbol[(gs2_coefs$Symbol %in% colnames(super_os))==F])
GS2_dfs_scores <- gs2_hpv(super_dfs, unmeasured_genes = gs2_coefs$Symbol[(gs2_coefs$Symbol %in% colnames(super_dfs))==F])

# scale the scores
GS2_score_os_scaled <- scale(GS2_os_scores)
GS2_score_dfs_scaled <- scale(GS2_dfs_scores)


# Analysis dataset: clinical variables and GS scores ----
super_os_analysis <- super_os %>% select(c(1:20)) # clinical variables
super_os_analysis <- super_os_analysis %>% mutate(GS2_score = GS2_score_os_scaled) #add score

super_dfs_analysis <- super_dfs %>% select(c(1:20))
super_dfs_analysis <- super_dfs_analysis %>% mutate(GS2_score = GS2_score_dfs_scaled)

# new variables needed ----
# os DATA: 0 1 instead of alive dead
super_os_analysis <- super_os_analysis %>% mutate(CensoredStatus = ifelse(
  follow_status_of_patient=="Alive",yes=0,no=1))

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


# Analysis: gs2_os_2y ----

# full model
gs2_os_2y <- coxph(Surv(osdays_2y,osstatus_2y)~ dataset+ GS2_score*HPV_status+tumor_region+clinical_age_at_diagnosis+
                     clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                     surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
                     chemo_chemotherapy_treatment, data=super_os_analysis, x=TRUE, model=T) 
summary(gs2_os_2y)
survMisc::rsq(gs2_os_2y) # R^2 

# Obtain marginal effect using relevel
super_os_analysis2 <- super_os_analysis %>% mutate(HPV_status = fct_relevel(HPV_status, "Negative"))
gs2_os_2y_relevel <- coxph(Surv(osdays_2y,osstatus_2y)~ dataset+ GS2_score*HPV_status+tumor_region+clinical_age_at_diagnosis+
                             clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                             surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
                             chemo_chemotherapy_treatment, data=super_os_analysis2, x=TRUE, model=T) 
summary(gs2_os_2y_relevel)

# equivalent clinical base model, comparable c index
summary( coxph(Surv(osdays_2y,osstatus_2y)~ dataset+ HPV_status+tumor_region+clinical_age_at_diagnosis+
        clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
        surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
        chemo_chemotherapy_treatment, data=super_os_analysis, x=TRUE, model=T) )
## c index 
survMisc::rsq(coxph(Surv(osdays_2y,osstatus_2y)~ dataset+ HPV_status+tumor_region+clinical_age_at_diagnosis+
                      clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                      surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
                      chemo_chemotherapy_treatment, data=super_os_analysis, x=TRUE, model=T) )
# R^2

# gs2_os_2y plots ----

# Plot the effect of GS
## make names for each facet/group, including the number of samples for that group
facet_names <- c(
  "Positive" = paste0("Positive (N = ", nrow(super_os_analysis[super_os_analysis$HPV_status=="Positive",]), ")" )
  , "Negative" = paste0("Negative (N = ", nrow(super_os_analysis[super_os_analysis$HPV_status=="Negative",]), ")" )
)

gs2_os_2y_plot <-  plot_surv_area(time="osdays_2y", status="osstatus_2y", variable="GS2_score", title = "HPV status", 
                                  xlab="Time since diagnosis (days)", ylab= "Overall survival probability",legend.title= "3 clusters HPV \nscore",
                                  group="HPV_status", data=super_os_analysis,model=gs2_os_2y,
                                  facet_args = list(labeller = as_labeller(facet_names))) + ylim(0,1)

gs2_os_2y_plot

# 2 histograms underneath
gs2_os_2y_plot_v2 <- gs2_os_2y_plot + theme(plot.tag = element_text()) + labs(tag = "a")

gg_hist_pos<- super_os_analysis %>% filter(HPV_status=="Positive") %>% ggplot(aes(GS2_score))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.tag = element_text()) +
  geom_histogram( color="#e9ecef", alpha=0.7, position = 'identity')+
  #scale_fill_manual(values=c("#69b3a2", "#404080")) +
  labs(title="HPV-positive",  tag = "b", y = "Count")+ xlim(c(-3, 7.5)) + 
  xlab("3 clusters HPV score") 

gg_hist_neg<- super_os_analysis %>% filter(HPV_status=="Negative") %>% ggplot(aes(GS2_score))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.tag = element_text()) +
  geom_histogram( color="#e9ecef", alpha=0.7, position = 'identity')+
  #scale_fill_manual(values=c("#69b3a2", "#404080")) +
  labs(title="HPV-negative", tag = "c", y = "Count")+ xlim(c(-3, 7.5)) + 
  xlab("3 clusters HPV score") 

gridExtra::grid.arrange(gs2_os_2y_plot_v2, gg_hist_pos, gg_hist_neg, 
                        layout_matrix = rbind(c(1, 1), c(2,3)), heights=c(4,1.5))
# Almost none are in the red area... Try building a model excluding these outliers and then plotting again 

# reduced dataset without patients with GS2 >3
super_os_analysis_reduced <- super_os_analysis %>% filter(GS2_score <= 3) # removes 2 cases
super_os_analysis_reduced_relevel <- super_os_analysis2 %>% filter(GS2_score <= 3) # removes 2 cases

# full model
gs2_os_2y_v2 <- coxph(Surv(osdays_2y,osstatus_2y)~ dataset+ GS2_score*HPV_status+tumor_region+clinical_age_at_diagnosis+
                     clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                     surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
                     chemo_chemotherapy_treatment, data=super_os_analysis_reduced, x=TRUE, model=T) 
summary(gs2_os_2y_v2)

summary(coxph(Surv(osdays_2y,osstatus_2y)~ dataset+ GS2_score*HPV_status+tumor_region+clinical_age_at_diagnosis+
        clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
        surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
        chemo_chemotherapy_treatment, data=super_os_analysis_reduced_relevel, x=TRUE, model=T))

# Very little difference when excluding "outliers" - thus should be ok to make plot with the reduced dataset

gs2_os_2y_v2_plot <-  plot_surv_area(time="osdays_2y", status="osstatus_2y", variable="GS2_score", title = "HPV status", 
                                  xlab="Time since diagnosis (days)", ylab= "Overall survival probability",legend.title= "3 clusters HPV \nscore",
                                  group="HPV_status", data=super_os_analysis_reduced,model=gs2_os_2y_v2,
                                  facet_args = list(labeller = as_labeller(facet_names))) + ylim(0,1)

# 2 histograms underneath
gs2_os_2y_v2_plot_v2 <- gs2_os_2y_v2_plot + theme(plot.tag = element_text()) + labs(tag = "a")

gg_hist_pos_v2<- super_os_analysis_reduced %>% filter(HPV_status=="Positive") %>% ggplot(aes(GS2_score))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.tag = element_text()) +
  geom_histogram( color="#e9ecef", alpha=0.7, position = 'identity')+
  #scale_fill_manual(values=c("#69b3a2", "#404080")) +
  labs(title="HPV-positive",  tag = "b", y = "Count")+ xlim(c(-3, 3)) + 
  xlab("3 clusters HPV score") 

gg_hist_neg_v2<- super_os_analysis_reduced %>% filter(HPV_status=="Negative") %>% ggplot(aes(GS2_score))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.tag = element_text()) +
  geom_histogram( color="#e9ecef", alpha=0.7, position = 'identity')+
  #scale_fill_manual(values=c("#69b3a2", "#404080")) +
  labs(title="HPV-negative", tag = "c", y = "Count")+ xlim(c(-3, 3)) + 
  xlab("3 clusters HPV score") 

gridExtra::grid.arrange(gs2_os_2y_v2_plot_v2, gg_hist_pos_v2, gg_hist_neg_v2, 
                        layout_matrix = rbind(c(1, 1), c(2,3)), heights=c(4,1.5))



# plot for publication ----
# Note: data=super_os_analysis_reduced, model=gs2_os_2y_v2 ; loaded under above plot section

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
gs2_os_2y_final_plot <-  plot_surv_area(time="osdays_2y", status="osstatus_2y", variable="GS2_score", title = "HPV status", 
                                     xlab="Time since diagnosis (days)", ylab= "Overall survival probability",legend.title= "3 clusters HPV \nscore",
                                     group="HPV_status", data=super_os_analysis_reduced,model=gs2_os_2y_v2,
                                     facet_args = list(labeller = new_labeller)) + ylim(0,1)



# 2 histograms underneath

gg_hist_pos_v2<- super_os_analysis_reduced %>% filter(HPV_status=="Positive") %>% ggplot(aes(GS2_score))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.tag = element_text()) +
  geom_histogram( color="#e9ecef", alpha=0.7, position = 'identity')+
  #scale_fill_manual(values=c("#69b3a2", "#404080")) +
  labs(title="HPV-positive",  tag = "(C)", y = "Count")+ xlim(c(-3, 3)) + 
  xlab("3 clusters HPV score") 

gg_hist_neg_v2<- super_os_analysis_reduced %>% filter(HPV_status=="Negative") %>% ggplot(aes(GS2_score))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.tag = element_text()) +
  geom_histogram( color="#e9ecef", alpha=0.7, position = 'identity')+
  #scale_fill_manual(values=c("#69b3a2", "#404080")) +
  labs(title="HPV-negative", tag = "(D)", y = "Count")+ xlim(c(-3, 3)) + 
  xlab("3 clusters HPV score") 

gridExtra::grid.arrange(gs2_os_2y_final_plot, gg_hist_pos_v2, gg_hist_neg_v2, 
                        layout_matrix = rbind(c(1, 1), c(2,3)), heights=c(4,1.5))


# make tiff file
file_path <- "M:/data_directory/GS2_plot.tiff"
tiff(file_path, units="mm", width=120, height=150, res=300)

gridExtra::grid.arrange(gs2_os_2y_final_plot, gg_hist_pos_v2, gg_hist_neg_v2, 
                        layout_matrix = rbind(c(1, 1), c(2,3)), heights=c(4,1.5))

dev.off()




# Analysis: gs1_os_5y ----
# full model
gs2_os_5y <- coxph(Surv(osdays_5y,osstatus_5y)~ dataset+ GS2_score*HPV_status+tumor_region+clinical_age_at_diagnosis+
                     clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                     surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
                     chemo_chemotherapy_treatment, data=super_os_analysis, x=TRUE, model=T) 
summary(gs2_os_5y)

# Obtain marginal effect using relevel
super_os_analysis2 <- super_os_analysis %>% mutate(HPV_status = fct_relevel(HPV_status, "Negative"))
gs2_os_5y_relevel <- coxph(Surv(osdays_5y,osstatus_5y)~ dataset+ GS2_score*HPV_status+tumor_region+clinical_age_at_diagnosis+
                             clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                             surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
                             chemo_chemotherapy_treatment, data=super_os_analysis2, x=TRUE, model=T) 
summary(gs2_os_5y_relevel)

# Analysis: gs2_dfs_2y ----

# full model
gs2_dfs_2y <- coxph(Surv(dfsdays_2y,dfs_event_2y)~ dataset+ GS2_score*HPV_status+tumor_region+clinical_age_at_diagnosis+
                      clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                      surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
                      chemo_chemotherapy_treatment, data=super_dfs_analysis, x=TRUE, model=T) 
summary(gs2_dfs_2y)
survMisc::rsq(gs2_dfs_2y) # r^2 


# Obtain marginal effect using relevel
super_dfs_analysis2 <- super_dfs_analysis %>% mutate(HPV_status = fct_relevel(HPV_status, "Negative"))
gs2_dfs_2y_relevel <- coxph(Surv(dfsdays_2y,dfs_event_2y)~ dataset+ GS2_score*HPV_status+tumor_region+clinical_age_at_diagnosis+
                              clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                              surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
                              chemo_chemotherapy_treatment, data=super_dfs_analysis2, x=TRUE, model=T)
summary(gs2_dfs_2y_relevel)


# equivalent clinical base model, comparable c index
summary( coxph(Surv(dfsdays_2y,dfs_event_2y)~ dataset+ HPV_status+tumor_region+clinical_age_at_diagnosis+
                 clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                 surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
                 chemo_chemotherapy_treatment, data=super_dfs_analysis, x=TRUE, model=T)  )
## c index 
survMisc::rsq(coxph(Surv(dfsdays_2y,dfs_event_2y)~ dataset+ HPV_status+tumor_region+clinical_age_at_diagnosis+
                      clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                      surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
                      chemo_chemotherapy_treatment, data=super_dfs_analysis, x=TRUE, model=T))


# Analysis: gs2_dfs_5y ----

# full model
gs2_dfs_5y <- coxph(Surv(dfsdays_5y,dfs_event_5y)~ dataset+ GS2_score*HPV_status+tumor_region+clinical_age_at_diagnosis+
                      clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                      surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
                      chemo_chemotherapy_treatment, data=super_dfs_analysis, x=TRUE, model=T) 
summary(gs2_dfs_5y)

# Obtain marginal effect using relevel
super_dfs_analysis2 <- super_dfs_analysis %>% mutate(HPV_status = fct_relevel(HPV_status, "Negative"))
gs2_dfs_5y_relevel <- coxph(Surv(dfsdays_5y,dfs_event_5y)~ dataset+ GS2_score*HPV_status+tumor_region+clinical_age_at_diagnosis+
                              clinical_sex+ ctn_disease_extension_diagnosis + smoking_category+
                              surge_undergone_cancer_surgery+radio_radiotherapy_treatment+
                              chemo_chemotherapy_treatment, data=super_dfs_analysis2, x=TRUE, model=T)
summary(gs2_dfs_5y_relevel)

