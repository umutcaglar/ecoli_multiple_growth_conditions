## DeSeq Normalization Cleaned 03

# The aim of the code is to generate normalized data matrix.
# The work flow compses of four parts
# Pick up the samples
# pick up the rows
# calculate size factors
# do the normalization


###*****************************
# INITIAL COMMANDS TO RESET THE SYSTEM
rm(list = ls())
if (is.integer(dev.list())){dev.off()}
cat("\014")
###*****************************


###*****************************
# Set Working Directory
# One needs to arrange the correct pathway if this is not umut's computer ;)
if(as.vector(Sys.info()["effective_user"]=="umut"))
{setwd(paste0("/Users/umut/GitHub/ecoli_multiple_growth_conditions/c_code_change_wrt_variables_RNA&Protein/test_DeSeq2/"))} # mac computer
###*****************************


###*****************************
# DOWNLOAD LIBRARIES
require("Biobase") 
require("DESeq2")
require("dplyr")
require("tidyr")
###*****************************


###*****************************
#Load Functions
source("data_filter_normalization_functions_fakeBatch.R")
###*****************************


###*****************************
saveFiles=TRUE
runDeSeqForDifExp=TRUE
# The data filtering function that controls sub functions.
mainData=filter_data(dataType = "mrna", # can be "rna", "mrna", "protein", "protein_wo_NA"
                     badDataSet = "set00", # can be "set00",set01","set02", "set03"
                     # referenceParameters can be a vector like
                     # c("growthPhase", "Mg_mM_Levels", "Na_mM_Levels", "carbonSource", "experiment")
                     referenceParameters=c("growthPhase",
                                           "Mg_mM_Levels", 
                                           "Na_mM_Levels", 
                                           "carbonSource", 
                                           "experiment"),
                     # referenceLevels can be a vector like
                     # c("exponential", "baseMg", "baseNa", "glucose", "glucose_time_course")
                     referenceLevels=c("exponential",
                                       "baseMg", 
                                       "baseNa", 
                                       "glucose", 
                                       "glucose_time_course"),
                     experimentVector = c("allEx"), # can be "Stc","Ytc","Nas","Agr","Ngr","Mgl","Mgh" // "allEx"
                     carbonSourceVector = "SYAN", # can be any sub combination of "SYAN"
                     MgLevelVector = c("allMg"), # can be "lowMg","baseMg","highMg" // "allMg"
                     NaLevelVector = c("allNa"), # can be "baseNa","highNa" // "allNa"
                     # can be "exponential","stationary","late_stationary" // "allPhase"
                     growthPhaseVector = c("allPhase"), 
                     filterGenes = "noFilter", # can be "noFilter", "meanFilter", "maxFilter", "sdFilter" 
                     threshold=NA, # the threshold value for "meanFilter", "maxFilter", "sdFilter"
                     roundData=TRUE,
                     sumTechnicalReplicates=TRUE,
                     deSeqSfChoice="p1Sf", # can be "regSf", "p1Sf", "noSf"
                     normalizationMethodChoice= "noNorm") # can be "vst", "rlog", "log10", "noNorm"
###*****************************


###*****************************
#Decompose the container
deseq_DataObj=mainData[[1]]
objectName=mainData[[2]]
###*****************************


###*****************************
# Test Multiple conditions to understand how deseq2 behave 
# (Numbers refer to excell sheet)
browser()

# Test DeSeq on batch

###*****************************
# Do the DeSeq2 test
# sample 1
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
test_for="Mg_mM_Levels"
DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ",test_for))
differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
(res_Mg <- DESeq2::results(object = differentialGeneAnalResults, pAdjustMethod ="fdr"))
mcols(res_Mg, use.names=TRUE)  
res_df_Mg<-as.data.frame(res_Mg)

DESeq2::summary.DESeqResults(object = res_Mg,alpha = 0.05)
###*****************************




###*****************************
# Do the DeSeq2 test
# sample 4
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
test_for="Mg_mM_Levels"
DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ",test_for," + batchNumber"))
differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
(res_MgB <- DESeq2::results(object = differentialGeneAnalResults, pAdjustMethod ="fdr"))
mcols(res_MgB, use.names=TRUE)  
res_df_MgB<-as.data.frame(res_MgB)

DESeq2::summary.DESeqResults(object = res_MgB,alpha = 0.05)
###*****************************



###*****************************
# Do the DeSeq2 test
# sample 5
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
test_for="Mg_mM_Levels"
DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ",test_for," + batchNumber"))
differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
(res_MgBconhb <- DESeq2::results(object = differentialGeneAnalResults, pAdjustMethod ="fdr", contrast = c(test_for,"highMg","baseMg")))
mcols(res_MgBconhb, use.names=TRUE)  
res_df_MgBconhb<-as.data.frame(res_MgBconhb)

DESeq2::summary.DESeqResults(object = res_MgBconhb,alpha = 0.05)
###*****************************


###*****************************
# Do the DeSeq2 test
# sample 6
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
test_for="Mg_mM_Levels"
DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ",test_for," + batchNumber"))
differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
(res_MgBconlb <- DESeq2::results(object = differentialGeneAnalResults, pAdjustMethod ="fdr", contrast = c(test_for,"lowMg","baseMg")))
mcols(res_MgBconlb, use.names=TRUE)  
res_df_MgBconlb<-as.data.frame(res_MgBconlb)

DESeq2::summary.DESeqResults(object = res_MgBconlb,alpha = 0.05)
###*****************************







###*****************************
# Do the DeSeq2 test
# sample 8
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
test_for="Na_mM_Levels"
DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ",test_for))
differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
(res_Na <- DESeq2::results(object = differentialGeneAnalResults, pAdjustMethod ="fdr"))
mcols(res_Na, use.names=TRUE)  
res_df_Na<-as.data.frame(res_Na)

DESeq2::summary.DESeqResults(object = res_Na,alpha = 0.05)
###*****************************


###*****************************
# Do the DeSeq2 test
# sample 9
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
test_for="Na_mM_Levels"
DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ",test_for))
differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
(res_NaCon <- DESeq2::results(object = differentialGeneAnalResults, pAdjustMethod ="fdr", contrast = c(test_for,"highNa","baseNa")))
mcols(res_NaCon, use.names=TRUE)  
res_df_NaCon<-as.data.frame(res_NaCon)

DESeq2::summary.DESeqResults(object = res_NaCon,alpha = 0.05)
###*****************************




###*****************************
# Do the DeSeq2 test
# sample 10
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
test_for="Na_mM_Levels"
DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ",test_for," + batchNumber"))
differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
(res_NaB <- DESeq2::results(object = differentialGeneAnalResults, pAdjustMethod ="fdr"))
mcols(res_NaB, use.names=TRUE)  
res_df_NaB<-as.data.frame(res_NaB)

DESeq2::summary.DESeqResults(object = res_NaB,alpha = 0.05)
###*****************************


###*****************************
# Do the DeSeq2 test
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
# sample 11
test_for="Na_mM_Levels"
DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ",test_for," + batchNumber"))
differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
(res_NaBCon <- DESeq2::results(object = differentialGeneAnalResults, pAdjustMethod ="fdr", contrast = c(test_for,"highNa","baseNa")))
mcols(res_NaBCon, use.names=TRUE)  
res_df_NaBCon<-as.data.frame(res_NaBCon)

DESeq2::summary.DESeqResults(object = res_NaBCon,alpha = 0.05)
###*****************************


###*****************************
# Do the DeSeq2 test
# sample 12
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
test_for="Na_mM_Levels"
DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ","batchNumber + ",test_for))
differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
(res_NaBrevCon <- DESeq2::results(object = differentialGeneAnalResults, pAdjustMethod ="fdr", contrast = c(test_for,"highNa","baseNa")))
mcols(res_NaBrevCon, use.names=TRUE)  
res_df_NaBrevCon<-as.data.frame(res_NaBrevCon)

DESeq2::summary.DESeqResults(object = res_NaBrevCon,alpha = 0.05)
###*****************************



###*****************************
# Do the DeSeq2 test
# sample 13
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
test_for="Na_mM_Levels"
DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ",test_for," + batchNumber"))
differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
(res_NaBConrev <- DESeq2::results(object = differentialGeneAnalResults, pAdjustMethod ="fdr", contrast = c(test_for,"baseNa","highNa")))
mcols(res_NaBConrev, use.names=TRUE)  
res_df_NaBConrev<-as.data.frame(res_NaBConrev)

DESeq2::summary.DESeqResults(object = res_NaBConrev,alpha = 0.05)
###*****************************








###*****************************
# Do the DeSeq2 test
# sample 16
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
test_for="carbonSource"
DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ",test_for))
differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
(res_CarbonCon <- DESeq2::results(object = differentialGeneAnalResults, pAdjustMethod ="fdr", contrast = c(test_for,"glycerol","glucose")))
mcols(res_CarbonCon, use.names=TRUE)  
res_df_CarbonCon<-as.data.frame(res_CarbonCon)

DESeq2::summary.DESeqResults(object = res_CarbonCon,alpha = 0.05)
###*****************************


###*****************************
# Do the DeSeq2 test
# sample 20
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
test_for="carbonSource"
DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ",test_for," + batchNumber"))
differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
(res_CarbonBCon <- DESeq2::results(object = differentialGeneAnalResults, pAdjustMethod ="fdr", contrast = c(test_for,"glycerol","glucose")))
mcols(res_CarbonBCon, use.names=TRUE)  
res_df_CarbonBCon<-as.data.frame(res_CarbonBCon)

DESeq2::summary.DESeqResults(object = res_CarbonBCon,alpha = 0.05)
###*****************************


###*****************************
# Do the DeSeq2 test
# sample 22
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
test_for="carbonSource"
DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ",test_for," + fake_batch"))
differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
(res_CarbonfBCon <- DESeq2::results(object = differentialGeneAnalResults, pAdjustMethod ="fdr", contrast = c(test_for,"glycerol","glucose")))
mcols(res_CarbonfBCon, use.names=TRUE)  
res_df_CarbonfBCon<-as.data.frame(res_CarbonfBCon)

DESeq2::summary.DESeqResults(object = res_CarbonfBCon,alpha = 0.05)
###*****************************
