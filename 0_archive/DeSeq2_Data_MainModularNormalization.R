## DeSeq Normalization Cleaned 02

# The aim of the code is to generate normalized data matrix by using previous analysis.
# It will clean out the prevously detected samples
# it will filter out some genes.


###*****************************
# INITIAL COMMANDS TO RESET THE SYSTEM

rm(list = ls())
if (is.integer(dev.list())){dev.off()}
cat("\014")

# Set Working Directory
setwd('/Users/umut/GitHub/ecoli_multiple_growth_conditions/a_code_dataPreperation_RNA&Protein/') # mac computer
###*****************************





###*****************************
# DOWNLOAD LIBRARIES

library("Biobase") 
library("DESeq2")
library("dplyr")
###*****************************



###*****************************
# LOAD DATA
unnormalized_rna_data=read.csv(file=paste0("../a_results/","rnaMatrix_RNA.csv")) #The main data matrix (includes everything [mRNA , tRNA, RF] )
unnormalized_mrna_data=read.csv(file=paste0("../a_results/","rnaMatrix_mRNA.csv")) #The main unnormalized mRNA matrix
unnormalized_protein_data=read.csv(file=paste0("../a_results/","proteinMatrix.csv")) # The main unnormalized protein matrix
unnormalized_protein_data_wo_NA=read.csv(file=paste0("../a_results/","proteinMatrix_wo_NA.csv")) # The unnormalized protein matrix contains proteins that are measured for all datasets
meta_rna=read.csv(paste0("../a_results/","metaRNA.csv")) #Main Meta data Frame
meta_protein=read.csv(paste0("../a_results/","metaProtein.csv")) #Main Meta data Frame
###*****************************


###*****************************
# FUNCTIONS
source("DeSeq2_Data_Functions.R")
###*****************************


###*****************************
# PREPEARE DATA

# Data matrix preperation for unnormalized rna data frame 
rownames(unnormalized_rna_data)<-unnormalized_rna_data[["gene_ID"]]
unnormalized_rna_data %>% 
  dplyr::select(-c(gene_Type, gene_ID))->unnormalized_rna_data

# Data matrix preperation for unnormalized mRNA data frame
rownames(unnormalized_mrna_data)<-unnormalized_mrna_data[["gene_ID"]]
unnormalized_mrna_data %>% 
  dplyr::select(-c(gene_Type, gene_ID))->unnormalized_mrna_data

# Data matrix preperation for unnormalized protein data frame
# includes decimal data need to be rounded
rownames(unnormalized_protein_data)<-unnormalized_protein_data[["gene_id"]]
unnormalized_protein_data %>% 
  dplyr::select(-c(gene_id))->unnormalized_protein_data

# Data matrix preperation for unnormalized protein data frame without NA
# includes decimal data need to be rounded
rownames(unnormalized_protein_data_wo_NA)<-unnormalized_protein_data_wo_NA[["gene_id"]]
unnormalized_protein_data_wo_NA %>% 
  dplyr::select(-c(gene_id))->unnormalized_protein_data_wo_NA
unnormalized_protein_data_wo_NA=round(unnormalized_protein_data_wo_NA)
unnormalized_protein_data_wo_NAx6=round(unnormalized_protein_data_wo_NA)*6
###*****************************


###*****************************
# Run the FUNCTION
finalData<-generateNormalizedData(unnormalized_rna_Input=unnormalized_rna_data, 
                       unnormalized_mrna_Input=unnormalized_mrna_data,
                       unnormalized_protein_Input=unnormalized_protein_data,
                       unnormalized_protein_Input_wo_NA=unnormalized_protein_data_wo_NA,
                       unnormalized_protein_Input_wo_NAx6=unnormalized_protein_data_wo_NAx6,
                       #meta_Input=meta_rna,  # meta_rna or meta_protein
                       conditionNumberChoice="uniqueCondition",  # "uniqueCondition", "uniqueCondition02"
                       badDataFilterSetChoice="set02", # "set00", set01" , "set02"
                       dataTypeChoice="protein_wo_NA", # can be "rna", "mrna","protein","protein_wo_NA","protein_wo_NAx6"
                       dataTimeChoice="wholeSet", # can be "wholeSet", "exponential", "stationary", "late_stationary"
                       MgLevelChoice="allMg",  # can be "lowMg", "midMg", "highMg", "allMg"
                       NaLevelChoice="allNa",  # can be "lowNa","highNa", "allNa"
                       carbonTypeChoice="SYAN", # a letter combination from the list "SYAN"  S (glucose), Y (glycerol), A (lactate), N (gluconate) 
                       filterTypeChoice="noFilter",  # can be "noFilter", "mean", "sd", "max" "threshold"
                       deSeqNormChoice="p1", # can be "p1", "p6", reg"
                       normalizationMethodChoice= "vst", # can be "vst", "log2", "none"
                       experimentChoice=c("allEx"))  # can be "allEx", 
                                            # or a combination of below
                                            # "Stc" for "glucose_time_course", 
                                            # "Ytc" for "glycerol_time_course", 
                                            # "Nas" for "NaCl_stress", 
                                            # "Agr" for "lactate_growth", 
                                            # "Ngr" for "gluconate_growth", 
                                            # "Mgl" for "MgSO4_stress_low", 
                                            # "Mgh" for "MgSO4_stress_high" 









