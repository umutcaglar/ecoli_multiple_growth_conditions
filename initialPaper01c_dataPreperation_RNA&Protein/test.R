# Comapring the 2 protein experiments

#"protein_wo_NA","protein_wo_NAx6"

###*****************************
# INITIAL COMMANDS TO RESET THE SYSTEM

rm(list = ls())
if (is.integer(dev.list())){dev.off()}
cat("\014")

# Set Working Directory
setwd('/Users/umut/GitHub/AG3C_Analyze/initialPaper01c_dataPreperation_RNA&Protein/') # mac computer
###*****************************


###*****************************
# DOWNLOAD LIBRARIES

library("Biobase") 
library("DESeq2")
library("dplyr")
###*****************************


###*****************************
load(file = "../initialPaper01r/normalized_none_p1_noFilter_protein_wo_NA_wholeSet_allMg_allNa_SYAN_set00_allEx.RData")
load(file = "../initialPaper01r/normalized_none_p6_noFilter_protein_wo_NAx6_wholeSet_allMg_allNa_SYAN_set00_allEx.RData")
###*****************************


###*****************************
q1=normalized_none_p1_noFilter_protein_wo_NA_wholeSet_allMg_allNa_SYAN_set00_allEx
q2=normalized_none_p6_noFilter_protein_wo_NAx6_wholeSet_allMg_allNa_SYAN_set00_allEx
###*****************************


###*****************************
m=hist(as.vector(q2-q1),50)
###*****************************

###*****************************
load(file = "../initialPaper01r/normalized_vst_p1_noFilter_protein_wo_NA_wholeSet_allMg_allNa_SYAN_set00_allEx.RData")
load(file = "../initialPaper01r/normalized_vst_p6_noFilter_protein_wo_NAx6_wholeSet_allMg_allNa_SYAN_set00_allEx.RData")
load(file = "../initialPaper01r/normalized_log2_p6_noFilter_protein_wo_NAx6_wholeSet_allMg_allNa_SYAN_set00_allEx.RData")
###*****************************


###*****************************
q3=normalized_vst_p1_noFilter_protein_wo_NA_wholeSet_allMg_allNa_SYAN_set00_allEx
q4=normalized_vst_p6_noFilter_protein_wo_NAx6_wholeSet_allMg_allNa_SYAN_set00_allEx
q5=normalized_log2_p6_noFilter_protein_wo_NAx6_wholeSet_allMg_allNa_SYAN_set00_allEx
###*****************************


###*****************************
m=hist(as.vector(q3-q4),50)
###*****************************