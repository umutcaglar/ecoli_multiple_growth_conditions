## DeSeq Name generate function

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
{setwd('/Users/umut/GitHub/ecoli_multiple_growth_conditions/c_code_change_wrt_variables_RNA&Protein/')} # mac computer
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
source("../a_code_dataPreperation_RNA&Protein/data_naming_functions.R")
###*****************************


###*****************************
# The main data naming function that controls sub functions.
dataName=name_data(initialValue="resDf", # can be c("genes0.05","genes_P0.05Fold2","resDf")
                   dataType = "mrna", # can be "rna", "mrna", "protein", "protein_wo_NA"
                   badDataSet = "set02", # can be "set00",set01","set02", "set03"
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
                   carbonSourceVector = "S", # can be any sub combination of "SYAN"
                   MgLevelVector = c("baseMg"), # can be "lowMg","baseMg","highMg" // "allMg"
                   NaLevelVector = c("allNa"), # can be "baseNa","highNa" // "allNa"
                   growthPhaseVector = c("exponential"), # can be "exponential","stationary","late_stationary" // "allPhase"
                   filterGenes = "noFilter", # can be "noFilter", "meanFilter", "maxFilter", "sdFilter" 
                   threshold=NA, # the threshold value for "meanFilter", "maxFilter", "sdFilter"
                   roundData=TRUE,
                   sumTechnicalReplicates=TRUE,
                   deSeqSfChoice="p1Sf", # can be "regSf", "p1Sf", "noSf"
                   normalizationMethodChoice= "noNorm", # can be "vst", "rlog", "log10", "noNorm"
                   test_for = "noTest")  # works only if normalizationMethodChoice == noNorm
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource", "noTest")
###*****************************

dataNameDF=as.data.frame(dataName[1])

metaDataName=dataNameDF
metaDataName$objectName.initial="metaData"
treeDataName=dataNameDF
treeDataName$objectName.initial="treeData"
heatMapName=dataNameDF
heatMapName$objectName.initial="heatMap"

dataName=paste(dataNameDF,collapse = "_")
metaDataName=paste(metaDataName,collapse = "_")
treeDataName=paste(treeDataName,collapse = "_")
heatMapName=paste(heatMapName,collapse = "_")

mainDataFrame=read.csv(file = paste0("../a_results/",dataName,".csv"),header = TRUE,row.names = 1)
condition=read.csv(file = paste0("../a_results/",metaDataName,".csv"),header = TRUE)


