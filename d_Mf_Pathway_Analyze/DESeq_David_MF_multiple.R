# Analayze Go Annotations Molecular function DESeq + DAVID

# Aim of the code is to find the genes in kegg pathways and send files to generate figures


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
{setwd(paste0("/Users/umut/GitHub/ecoli_multiple_growth_conditions/",
              "d_Mf_Pathway_Analyze/"))} # mac computer
###*****************************


###*****************************
# DOWNLOAD LIBRARIES
require("Biobase") 
require("DESeq2")
require("dplyr")
require("tidyr")
require("ggplot2")
require("RColorBrewer")
require("scales")
require("cowplot")
require("ggrepel")
###*****************************


###*****************************
# Load Files
source("../a_code_dataPreperation_RNA&Protein/data_naming_functions.R")
source("../b_code_histogram_RNA&Protein/replace_fun.R")
###*****************************


###*****************************
# Add Input List
inputList1=list(initialValueInput="genes_P0.05Fold2",dataTypeInput = "mrna",badDataSetInput = "set02",
                referenceParametersInput=c("growthPhase","Mg_mM_Levels","Na_mM_Levels","carbonSource","experiment"),
                referenceLevelsInput=c("exponential","baseMg","baseNa","glucose","glucose_time_course"),
                experimentVectorInput = c("allEx"),
                carbonSourceVectorInput = "S",growthPhaseVectorInput = c("exponential"),
                MgLevelVectorInput = c("baseMg","highMg"),NaLevelVectorInput = c("baseNa"), 
                filterGenesInput = "noFilter",
                thresholdInput=NA, 
                roundDataInput=TRUE,sumTechnicalReplicatesInput=TRUE,deSeqSfChoiceInput="p1Sf",
                normalizationMethodChoiceInput= "noNorm",
                test_forInput = "Mg_mM_Levels")

inputList2=list(initialValueInput="genes_P0.05Fold2",dataTypeInput = "mrna",badDataSetInput = "set02",
                referenceParametersInput=c("growthPhase","Mg_mM_Levels","Na_mM_Levels","carbonSource","experiment"),
                referenceLevelsInput=c("exponential","baseMg","baseNa","glucose","glucose_time_course"),
                experimentVectorInput = c("allEx"),
                carbonSourceVectorInput = "S",growthPhaseVectorInput = c("exponential"),
                MgLevelVectorInput = c("baseMg","lowMg"),NaLevelVectorInput = c("baseNa"), 
                filterGenesInput = "noFilter",
                thresholdInput=NA, 
                roundDataInput=TRUE,sumTechnicalReplicatesInput=TRUE,deSeqSfChoiceInput="p1Sf",
                normalizationMethodChoiceInput= "noNorm",
                test_forInput = "Mg_mM_Levels")

inputList3=list(initialValueInput="genes_P0.05Fold2",dataTypeInput = "mrna",badDataSetInput = "set02",
                referenceParametersInput=c("growthPhase","Mg_mM_Levels","Na_mM_Levels","carbonSource","experiment"),
                referenceLevelsInput=c("exponential","baseMg","baseNa","glucose","glucose_time_course"),
                experimentVectorInput = c("allEx"),
                carbonSourceVectorInput = "S",growthPhaseVectorInput = c("exponential"),
                MgLevelVectorInput = c("baseMg"),NaLevelVectorInput = c("allNa"), 
                filterGenesInput = "noFilter",
                thresholdInput=NA, 
                roundDataInput=TRUE,sumTechnicalReplicatesInput=TRUE,deSeqSfChoiceInput="p1Sf",
                normalizationMethodChoiceInput= "noNorm",
                test_forInput = "Na_mM_Levels")

inputList4=list(initialValueInput="genes_P0.05Fold2",dataTypeInput = "mrna",badDataSetInput = "set02",
                referenceParametersInput=c("growthPhase","Mg_mM_Levels","Na_mM_Levels","carbonSource","experiment"),
                referenceLevelsInput=c("exponential","baseMg","baseNa","glucose","glucose_time_course"),
                experimentVectorInput = c("allEx"),
                carbonSourceVectorInput = "SY",growthPhaseVectorInput = c("exponential"),
                MgLevelVectorInput = c("baseMg"),NaLevelVectorInput = c("baseNa"), 
                filterGenesInput = "noFilter",
                thresholdInput=NA, 
                roundDataInput=TRUE,sumTechnicalReplicatesInput=TRUE,deSeqSfChoiceInput="p1Sf",
                normalizationMethodChoiceInput= "noNorm",
                test_forInput = "carbonSource")

inputList5=list(initialValueInput="genes_P0.05Fold2",dataTypeInput = "mrna",badDataSetInput = "set02",
                referenceParametersInput=c("growthPhase","Mg_mM_Levels","Na_mM_Levels","carbonSource","experiment"),
                referenceLevelsInput=c("exponential","baseMg","baseNa","glucose","glucose_time_course"),
                experimentVectorInput = c("allEx"),
                carbonSourceVectorInput = "SN",growthPhaseVectorInput = c("exponential"),
                MgLevelVectorInput = c("baseMg"),NaLevelVectorInput = c("baseNa"), 
                filterGenesInput = "noFilter",
                thresholdInput=NA, 
                roundDataInput=TRUE,sumTechnicalReplicatesInput=TRUE,deSeqSfChoiceInput="p1Sf",
                normalizationMethodChoiceInput= "noNorm",
                test_forInput = "carbonSource")

inputList6=list(initialValueInput="genes_P0.05Fold2",dataTypeInput = "mrna",badDataSetInput = "set02",
                referenceParametersInput=c("growthPhase","Mg_mM_Levels","Na_mM_Levels","carbonSource","experiment"),
                referenceLevelsInput=c("exponential","baseMg","baseNa","glucose","glucose_time_course"),
                experimentVectorInput = c("allEx"),
                carbonSourceVectorInput = "SA",growthPhaseVectorInput = c("exponential"),
                MgLevelVectorInput = c("baseMg"),NaLevelVectorInput = c("baseNa"), 
                filterGenesInput = "noFilter",
                thresholdInput=NA, 
                roundDataInput=TRUE,sumTechnicalReplicatesInput=TRUE,deSeqSfChoiceInput="p1Sf",
                normalizationMethodChoiceInput= "noNorm",
                test_forInput = "carbonSource")


inputList7=list(initialValueInput="genes_P0.05Fold2",dataTypeInput = "mrna",badDataSetInput = "set02",
                referenceParametersInput=c("growthPhase","Mg_mM_Levels","Na_mM_Levels","carbonSource","experiment"),
                referenceLevelsInput=c("stationary","baseMg","baseNa","glucose","glucose_time_course"),
                experimentVectorInput = c("allEx"),
                carbonSourceVectorInput = "S",growthPhaseVectorInput = c("stationary"),
                MgLevelVectorInput = c("baseMg","highMg"),NaLevelVectorInput = c("baseNa"), 
                filterGenesInput = "noFilter",
                thresholdInput=NA, 
                roundDataInput=TRUE,sumTechnicalReplicatesInput=TRUE,deSeqSfChoiceInput="p1Sf",
                normalizationMethodChoiceInput= "noNorm",
                test_forInput = "Mg_mM_Levels")

inputList8=list(initialValueInput="genes_P0.05Fold2",dataTypeInput = "mrna",badDataSetInput = "set02",
                referenceParametersInput=c("growthPhase","Mg_mM_Levels","Na_mM_Levels","carbonSource","experiment"),
                referenceLevelsInput=c("stationary","baseMg","baseNa","glucose","glucose_time_course"),
                experimentVectorInput = c("allEx"),
                carbonSourceVectorInput = "S",growthPhaseVectorInput = c("stationary"),
                MgLevelVectorInput = c("baseMg","lowMg"),NaLevelVectorInput = c("baseNa"), 
                filterGenesInput = "noFilter",
                thresholdInput=NA, 
                roundDataInput=TRUE,sumTechnicalReplicatesInput=TRUE,deSeqSfChoiceInput="p1Sf",
                normalizationMethodChoiceInput= "noNorm",
                test_forInput = "Mg_mM_Levels")

inputList9=list(initialValueInput="genes_P0.05Fold2",dataTypeInput = "mrna",badDataSetInput = "set02",
                referenceParametersInput=c("growthPhase","Mg_mM_Levels","Na_mM_Levels","carbonSource","experiment"),
                referenceLevelsInput=c("stationary","baseMg","baseNa","glucose","glucose_time_course"),
                experimentVectorInput = c("allEx"),
                carbonSourceVectorInput = "S",growthPhaseVectorInput = c("stationary"),
                MgLevelVectorInput = c("baseMg"),NaLevelVectorInput = c("allNa"), 
                filterGenesInput = "noFilter",
                thresholdInput=NA, 
                roundDataInput=TRUE,sumTechnicalReplicatesInput=TRUE,deSeqSfChoiceInput="p1Sf",
                normalizationMethodChoiceInput= "noNorm",
                test_forInput = "Na_mM_Levels")

inputList10=list(initialValueInput="genes_P0.05Fold2",dataTypeInput = "mrna",badDataSetInput = "set02",
                 referenceParametersInput=c("growthPhase","Mg_mM_Levels","Na_mM_Levels","carbonSource","experiment"),
                 referenceLevelsInput=c("stationary","baseMg","baseNa","glucose","glucose_time_course"),
                 experimentVectorInput = c("allEx"),
                 carbonSourceVectorInput = "SY",growthPhaseVectorInput = c("stationary"),
                 MgLevelVectorInput = c("baseMg"),NaLevelVectorInput = c("baseNa"), 
                 filterGenesInput = "noFilter",
                 thresholdInput=NA, 
                 roundDataInput=TRUE,sumTechnicalReplicatesInput=TRUE,deSeqSfChoiceInput="p1Sf",
                 normalizationMethodChoiceInput= "noNorm",
                 test_forInput = "carbonSource")

inputList11=list(initialValueInput="genes_P0.05Fold2",dataTypeInput = "mrna",badDataSetInput = "set02",
                 referenceParametersInput=c("growthPhase","Mg_mM_Levels","Na_mM_Levels","carbonSource","experiment"),
                 referenceLevelsInput=c("stationary","baseMg","baseNa","glucose","glucose_time_course"),
                 experimentVectorInput = c("allEx"),
                 carbonSourceVectorInput = "SN",growthPhaseVectorInput = c("stationary"),
                 MgLevelVectorInput = c("baseMg"),NaLevelVectorInput = c("baseNa"), 
                 filterGenesInput = "noFilter",
                 thresholdInput=NA, 
                 roundDataInput=TRUE,sumTechnicalReplicatesInput=TRUE,deSeqSfChoiceInput="p1Sf",
                 normalizationMethodChoiceInput= "noNorm",
                 test_forInput = "carbonSource")

inputList12=list(initialValueInput="genes_P0.05Fold2",dataTypeInput = "mrna",badDataSetInput = "set02",
                 referenceParametersInput=c("growthPhase","Mg_mM_Levels","Na_mM_Levels","carbonSource","experiment"),
                 referenceLevelsInput=c("stationary","baseMg","baseNa","glucose","glucose_time_course"),
                 experimentVectorInput = c("allEx"),
                 carbonSourceVectorInput = "SA",growthPhaseVectorInput = c("stationary"),
                 MgLevelVectorInput = c("baseMg"),NaLevelVectorInput = c("baseNa"), 
                 filterGenesInput = "noFilter",
                 thresholdInput=NA, 
                 roundDataInput=TRUE,sumTechnicalReplicatesInput=TRUE,deSeqSfChoiceInput="p1Sf",
                 normalizationMethodChoiceInput= "noNorm",
                 test_forInput = "carbonSource")

inputList13=list(initialValueInput="genes_P0.05Fold2",dataTypeInput = "protein",badDataSetInput = "set00",
                 referenceParametersInput=c("growthPhase","Mg_mM_Levels","Na_mM_Levels","carbonSource","experiment"),
                 referenceLevelsInput=c("exponential","baseMg","baseNa","glucose","glucose_time_course"),
                 experimentVectorInput = c("allEx"),
                 carbonSourceVectorInput = "S",growthPhaseVectorInput = c("exponential"),
                 MgLevelVectorInput = c("baseMg","highMg"),NaLevelVectorInput = c("baseNa"), 
                 filterGenesInput = "noFilter",
                 thresholdInput=NA, 
                 roundDataInput=TRUE,sumTechnicalReplicatesInput=TRUE,deSeqSfChoiceInput="p1Sf",
                 normalizationMethodChoiceInput= "noNorm",
                 test_forInput = "Mg_mM_Levels")

inputList14=list(initialValueInput="genes_P0.05Fold2",dataTypeInput = "protein",badDataSetInput = "set00",
                 referenceParametersInput=c("growthPhase","Mg_mM_Levels","Na_mM_Levels","carbonSource","experiment"),
                 referenceLevelsInput=c("exponential","baseMg","baseNa","glucose","glucose_time_course"),
                 experimentVectorInput = c("allEx"),
                 carbonSourceVectorInput = "S",growthPhaseVectorInput = c("exponential"),
                 MgLevelVectorInput = c("baseMg","lowMg"),NaLevelVectorInput = c("baseNa"), 
                 filterGenesInput = "noFilter",
                 thresholdInput=NA, 
                 roundDataInput=TRUE,sumTechnicalReplicatesInput=TRUE,deSeqSfChoiceInput="p1Sf",
                 normalizationMethodChoiceInput= "noNorm",
                 test_forInput = "Mg_mM_Levels")

inputList15=list(initialValueInput="genes_P0.05Fold2",dataTypeInput = "protein",badDataSetInput = "set00",
                 referenceParametersInput=c("growthPhase","Mg_mM_Levels","Na_mM_Levels","carbonSource","experiment"),
                 referenceLevelsInput=c("exponential","baseMg","baseNa","glucose","glucose_time_course"),
                 experimentVectorInput = c("allEx"),
                 carbonSourceVectorInput = "S",growthPhaseVectorInput = c("exponential"),
                 MgLevelVectorInput = c("baseMg"),NaLevelVectorInput = c("allNa"), 
                 filterGenesInput = "noFilter",
                 thresholdInput=NA, 
                 roundDataInput=TRUE,sumTechnicalReplicatesInput=TRUE,deSeqSfChoiceInput="p1Sf",
                 normalizationMethodChoiceInput= "noNorm",
                 test_forInput = "Na_mM_Levels")

inputList16=list(initialValueInput="genes_P0.05Fold2",dataTypeInput = "protein",badDataSetInput = "set00",
                 referenceParametersInput=c("growthPhase","Mg_mM_Levels","Na_mM_Levels","carbonSource","experiment"),
                 referenceLevelsInput=c("exponential","baseMg","baseNa","glucose","glucose_time_course"),
                 experimentVectorInput = c("allEx"),
                 carbonSourceVectorInput = "SY",growthPhaseVectorInput = c("exponential"),
                 MgLevelVectorInput = c("baseMg"),NaLevelVectorInput = c("baseNa"), 
                 filterGenesInput = "noFilter",
                 thresholdInput=NA, 
                 roundDataInput=TRUE,sumTechnicalReplicatesInput=TRUE,deSeqSfChoiceInput="p1Sf",
                 normalizationMethodChoiceInput= "noNorm",
                 test_forInput = "carbonSource")

inputList17=list(initialValueInput="genes_P0.05Fold2",dataTypeInput = "protein",badDataSetInput = "set00",
                 referenceParametersInput=c("growthPhase","Mg_mM_Levels","Na_mM_Levels","carbonSource","experiment"),
                 referenceLevelsInput=c("exponential","baseMg","baseNa","glucose","glucose_time_course"),
                 experimentVectorInput = c("allEx"),
                 carbonSourceVectorInput = "SN",growthPhaseVectorInput = c("exponential"),
                 MgLevelVectorInput = c("baseMg"),NaLevelVectorInput = c("baseNa"), 
                 filterGenesInput = "noFilter",
                 thresholdInput=NA, 
                 roundDataInput=TRUE,sumTechnicalReplicatesInput=TRUE,deSeqSfChoiceInput="p1Sf",
                 normalizationMethodChoiceInput= "noNorm",
                 test_forInput = "carbonSource")

inputList18=list(initialValueInput="genes_P0.05Fold2",dataTypeInput = "protein",badDataSetInput = "set00",
                 referenceParametersInput=c("growthPhase","Mg_mM_Levels","Na_mM_Levels","carbonSource","experiment"),
                 referenceLevelsInput=c("exponential","baseMg","baseNa","glucose","glucose_time_course"),
                 experimentVectorInput = c("allEx"),
                 carbonSourceVectorInput = "SA",growthPhaseVectorInput = c("exponential"),
                 MgLevelVectorInput = c("baseMg"),NaLevelVectorInput = c("baseNa"), 
                 filterGenesInput = "noFilter",
                 thresholdInput=NA, 
                 roundDataInput=TRUE,sumTechnicalReplicatesInput=TRUE,deSeqSfChoiceInput="p1Sf",
                 normalizationMethodChoiceInput= "noNorm",
                 test_forInput = "carbonSource")


inputList19=list(initialValueInput="genes_P0.05Fold2",dataTypeInput = "protein",badDataSetInput = "set00",
                 referenceParametersInput=c("growthPhase","Mg_mM_Levels","Na_mM_Levels","carbonSource","experiment"),
                 referenceLevelsInput=c("stationary","baseMg","baseNa","glucose","glucose_time_course"),
                 experimentVectorInput = c("allEx"),
                 carbonSourceVectorInput = "S",growthPhaseVectorInput = c("stationary"),
                 MgLevelVectorInput = c("baseMg","highMg"),NaLevelVectorInput = c("baseNa"), 
                 filterGenesInput = "noFilter",
                 thresholdInput=NA, 
                 roundDataInput=TRUE,sumTechnicalReplicatesInput=TRUE,deSeqSfChoiceInput="p1Sf",
                 normalizationMethodChoiceInput= "noNorm",
                 test_forInput = "Mg_mM_Levels")

inputList20=list(initialValueInput="genes_P0.05Fold2",dataTypeInput = "protein",badDataSetInput = "set00",
                 referenceParametersInput=c("growthPhase","Mg_mM_Levels","Na_mM_Levels","carbonSource","experiment"),
                 referenceLevelsInput=c("stationary","baseMg","baseNa","glucose","glucose_time_course"),
                 experimentVectorInput = c("allEx"),
                 carbonSourceVectorInput = "S",growthPhaseVectorInput = c("stationary"),
                 MgLevelVectorInput = c("baseMg","lowMg"),NaLevelVectorInput = c("baseNa"), 
                 filterGenesInput = "noFilter",
                 thresholdInput=NA, 
                 roundDataInput=TRUE,sumTechnicalReplicatesInput=TRUE,deSeqSfChoiceInput="p1Sf",
                 normalizationMethodChoiceInput= "noNorm",
                 test_forInput = "Mg_mM_Levels")

inputList21=list(initialValueInput="genes_P0.05Fold2",dataTypeInput = "protein",badDataSetInput = "set00",
                 referenceParametersInput=c("growthPhase","Mg_mM_Levels","Na_mM_Levels","carbonSource","experiment"),
                 referenceLevelsInput=c("stationary","baseMg","baseNa","glucose","glucose_time_course"),
                 experimentVectorInput = c("allEx"),
                 carbonSourceVectorInput = "S",growthPhaseVectorInput = c("stationary"),
                 MgLevelVectorInput = c("baseMg"),NaLevelVectorInput = c("allNa"), 
                 filterGenesInput = "noFilter",
                 thresholdInput=NA, 
                 roundDataInput=TRUE,sumTechnicalReplicatesInput=TRUE,deSeqSfChoiceInput="p1Sf",
                 normalizationMethodChoiceInput= "noNorm",
                 test_forInput = "Na_mM_Levels")

inputList22=list(initialValueInput="genes_P0.05Fold2",dataTypeInput = "protein",badDataSetInput = "set00",
                 referenceParametersInput=c("growthPhase","Mg_mM_Levels","Na_mM_Levels","carbonSource","experiment"),
                 referenceLevelsInput=c("stationary","baseMg","baseNa","glucose","glucose_time_course"),
                 experimentVectorInput = c("allEx"),
                 carbonSourceVectorInput = "SY",growthPhaseVectorInput = c("stationary"),
                 MgLevelVectorInput = c("baseMg"),NaLevelVectorInput = c("baseNa"), 
                 filterGenesInput = "noFilter",
                 thresholdInput=NA, 
                 roundDataInput=TRUE,sumTechnicalReplicatesInput=TRUE,deSeqSfChoiceInput="p1Sf",
                 normalizationMethodChoiceInput= "noNorm",
                 test_forInput = "carbonSource")

inputList23=list(initialValueInput="genes_P0.05Fold2",dataTypeInput = "protein",badDataSetInput = "set00",
                 referenceParametersInput=c("growthPhase","Mg_mM_Levels","Na_mM_Levels","carbonSource","experiment"),
                 referenceLevelsInput=c("stationary","baseMg","baseNa","glucose","glucose_time_course"),
                 experimentVectorInput = c("allEx"),
                 carbonSourceVectorInput = "SN",growthPhaseVectorInput = c("stationary"),
                 MgLevelVectorInput = c("baseMg"),NaLevelVectorInput = c("baseNa"), 
                 filterGenesInput = "noFilter",
                 thresholdInput=NA, 
                 roundDataInput=TRUE,sumTechnicalReplicatesInput=TRUE,deSeqSfChoiceInput="p1Sf",
                 normalizationMethodChoiceInput= "noNorm",
                 test_forInput = "carbonSource")

inputList24=list(initialValueInput="genes_P0.05Fold2",dataTypeInput = "protein",badDataSetInput = "set00",
                 referenceParametersInput=c("growthPhase","Mg_mM_Levels","Na_mM_Levels","carbonSource","experiment"),
                 referenceLevelsInput=c("stationary","baseMg","baseNa","glucose","glucose_time_course"),
                 experimentVectorInput = c("allEx"),
                 carbonSourceVectorInput = "SA",growthPhaseVectorInput = c("stationary"),
                 MgLevelVectorInput = c("baseMg"),NaLevelVectorInput = c("baseNa"), 
                 filterGenesInput = "noFilter",
                 thresholdInput=NA, 
                 roundDataInput=TRUE,sumTechnicalReplicatesInput=TRUE,deSeqSfChoiceInput="p1Sf",
                 normalizationMethodChoiceInput= "noNorm",
                 test_forInput = "carbonSource")


inputList<-rbind(inputList1,inputList2,inputList3,inputList4,inputList5,inputList6,
                 inputList7,inputList8,inputList9,inputList10,inputList11,inputList12,
                 inputList13,inputList14,inputList15,inputList16,inputList17,inputList18,
                 inputList19,inputList20,inputList21,inputList22,inputList23,inputList24)
###*****************************




###*****************************
# Beginning of the loop
for(counter02 in c(1:5,6,7,8,9:24))
{
  ###*****************************
  # Download the DAVID input and output
  dataName=name_data(initialValue=inputList[[counter02,"initialValueInput"]], # can be "genes0.05", "genes_P0.05Fold2"
                     dataType = inputList[[counter02,"dataTypeInput"]], # can be "rna", "mrna", "protein", "protein_wo_NA"
                     badDataSet = inputList[[counter02,"badDataSetInput"]], # can be "set00",set01","set02", "set03"
                     # referenceParameters can be a vector like
                     # c("growthPhase", "Mg_mM_Levels", "Na_mM_Levels", "carbonSource", "experiment")
                     referenceParameters=inputList[[counter02,"referenceParametersInput"]],
                     # referenceLevels can be a vector like
                     # c("exponential", "baseMg", "baseNa", "glucose", "glucose_time_course")
                     referenceLevels=inputList[[counter02,"referenceLevelsInput"]],
                     experimentVector = inputList[[counter02,"experimentVectorInput"]], 
                     # can be "Stc","Ytc","Nas","Agr","Ngr","Mgl","Mgh" // "allEx"
                     carbonSourceVector = inputList[[counter02,"carbonSourceVectorInput"]], 
                     # can be any sub combination of "SYAN"
                     MgLevelVector = inputList[[counter02,"MgLevelVectorInput"]], 
                     # can be "lowMg","baseMg","highMg" // "allMg"
                     NaLevelVector = inputList[[counter02,"NaLevelVectorInput"]], 
                     # can be "baseNa","highNa" // "allNa"
                     growthPhaseVector = inputList[[counter02,"growthPhaseVectorInput"]], 
                     # can be "exponential","stationary","late_stationary" // "allPhase"
                     filterGenes = inputList[[counter02,"filterGenesInput"]], 
                     # can be "noFilter", "meanFilter", "maxFilter", "sdFilter" 
                     threshold=inputList[[counter02,"thresholdInput"]], 
                     # the threshold value for "meanFilter", "maxFilter", "sdFilter"
                     roundData=inputList[[counter02,"roundDataInput"]],
                     sumTechnicalReplicates=inputList[[counter02,"sumTechnicalReplicatesInput"]],
                     deSeqSfChoice=inputList[[counter02,"deSeqSfChoiceInput"]], # can be "regSf", "p1Sf"
                     normalizationMethodChoice= inputList[[counter02,"normalizationMethodChoiceInput"]], 
                     # can be "vst", "rlog", "log10", "noNorm"
                     test_for = inputList[[counter02,"test_forInput"]])  
  # works only if normalizationMethodChoice == noNorm
  # c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
  
  objectName=paste(dataName$objectName,collapse = "_")
  ###*****************************
  
  
  ###*****************************
  if(paste0(objectName,"_mf.txt")%in%dir("../c_results/david_results/"))
  {
    GoMF_input<-read.csv(file = paste0("../c_results/DeSeq2_diffGene_Results/",objectName,".csv"),
                         header = TRUE)
    objectName_df=dataName$objectName
    objectName_df$initial="resDf"
    objectName_df=paste(objectName_df,collapse = "_")
    GoMF_input_df<-read.csv(file = paste0("../c_results/DeSeq2_diffGene_Results/",objectName_df,".csv"),
                            header = TRUE)
    
    GoMF_result<-read.table(file=paste0("../c_results/david_results/",objectName,"_mf.txt"),
                            sep = "\t",header = TRUE)
    
    names(GoMF_result)[names(GoMF_result) == 'Term'] <- 'MF_Name'
    ###*****************************
    
    
    ###*****************************
    # load reference libraries
    dataType=dataName$objectName$pick_data
    
    if(dataType %in% c("rna","mrna"))
    {
      MF_david_ref_2009_tidy<-read.table(file = paste0("../d_Mf_Pathway_Analyze/ReferenceFiles/",
                                                       "MF_david_mrna_ref_2009_tidy.txt"),
                                         sep = "\t",header = TRUE,fill = TRUE, quote = "")
    }
    
    
    if(dataType %in% c("protein_wo_NA", "protein"))
    {
      MF_david_ref_2009_tidy<-read.table(paste0("../d_Mf_Pathway_Analyze/ReferenceFiles/",
                                                "MF_david_protein_ref_2009_tidy.txt"),
                                         sep = "\t",header = TRUE,fill = TRUE, quote = "")
    }
    ###*****************************
    
    
    ###*****************************
    # Find MF output in ref 
    # do the consistency check if everything is in
    
    for(counter01 in 1:nrow(GoMF_result))
    {
      theGoMF=as.vector(GoMF_result[["MF_Name"]][counter01])
      MF_david_ref_2009_tidy %>%
        dplyr::filter(MF_Name==theGoMF)->temp
      
      m2<-as.vector(temp$ID)
      m1<-strsplit(x=as.vector(GoMF_result$Genes[counter01]),split = ",")
      m1<-sub(" ","",noquote(m1[[1]]))
      m1<-paste0(tolower(substr(start = 1,stop = 3,x = m1)),substr(start = 4,stop = 4,x = m1))
      print(all(m1 %in% m2))
    }
    ###*****************************
    
    
    ###*****************************
    # Find go annotations (MF) // 
    # select go annotation (MF) related genes 
    # // find the genes in our data related with Go annotation//
    GoMF_result ->GoMF_result_short
    GoMF_list=as.vector(GoMF_result_short$MF_Name)
    inputGenes=unique(as.vector(GoMF_input_df$gene_name))
    GoMF_input_df %>%
      dplyr::select(ID=gene_name,padj_gene=padj, log2=log2FoldChange, signChange)%>%
      dplyr::mutate(score_gene=-signChange*log10(padj_gene))->GoMF_input_df_narrow
    
    GoMF_result_short %>%
      dplyr::select(MF_Name, FDR_GoMF=FDR)%>%
      dplyr::mutate(FDR_GoMF=FDR_GoMF/100)->GoMF_result_narrow
    
    MF_david_ref_2009_tidy %>%
      dplyr::select(ID, MF_Name)%>%
      dplyr::filter(MF_Name %in% GoMF_list)%>%
      dplyr::group_by(MF_Name)%>%
      dplyr::summarise(numGenesInCat=length(MF_Name))->numGenesInCat_df
    
    
    MF_david_ref_2009_tidy %>%
      dplyr::select(ID, MF_Name)%>%
      dplyr::filter(MF_Name %in% GoMF_list,
                    ID %in% inputGenes) %>%
      dplyr::left_join(.,GoMF_input_df_narrow) %>%
      dplyr::left_join(.,GoMF_result_narrow) %>%
      dplyr::left_join(.,numGenesInCat_df) %>%
      dplyr::distinct()%>%
      dplyr::filter(!is.na(padj_gene))%>%
      dplyr::group_by(MF_Name)%>%
      dplyr::arrange(desc(score_gene))->selectedDf
    
    
    
    # add limitations for figures
    selectedDf %>%
      dplyr::filter(padj_gene<0.05,FDR_GoMF<0.05, abs(log2)>1) ->selectedDf
    
    if(nrow(selectedDf)!=0)
    {
      selectedDf %>%
        dplyr::mutate(abs_score=abs(score_gene))%>%
        dplyr::group_by(MF_Name,signChange)%>%
        dplyr::arrange(abs_score)%>%
        dplyr::mutate(rank=signChange*seq(1,n()))%>%
        dplyr::group_by(MF_Name)%>%
        dplyr::mutate(numSigP=(max(rank)+abs(max(rank)))/2,
                      numSigN=abs(min(rank)))%>%
        dplyr::group_by(ID,MF_Name)%>%
        dplyr::mutate(MF_Name_long=paste0(sub(".*~","",MF_Name),
                                          "\n padj:",
                                          sprintf("%.5f", FDR_GoMF),
                                          " N( -",numSigN,"/ +",numSigP,"/ ",numGenesInCat,")"))%>%
        dplyr::mutate(MF_Name_short=paste0(sub(".*~","",MF_Name)))%>%
        dplyr::group_by(MF_Name)%>%
        dplyr::arrange(desc(score_gene)) -> selectedDf
      
      selectedDf %>%
        dplyr::group_by(MF_Name_long)%>%
        dplyr::summarise(FDR_GoMF=unique(FDR_GoMF))%>%
        dplyr::arrange(FDR_GoMF)->summary_df
      
      as.vector(summary_df$MF_Name_long)
      
      
      selectedDf$MF_Name_long <- factor(selectedDf$MF_Name_long, 
                                        levels = rev(as.vector(summary_df$MF_Name_long)))
      ###*****************************
      
      
      ###*****************************
      # Generate simple Data Frame
      # Additional Parameters
      maxPathway=20
      maxGene=15
      
      if(length(unique(as.vector(selectedDf$FDR_GoMF)))<maxPathway)
      {maxPathway=length(unique(as.vector(selectedDf$FDR_GoMF)))}
      
      FDR_GoMFTopn=sort(unique(as.vector(selectedDf$FDR_GoMF)))[maxPathway]
      selectedDf %>%
        dplyr::group_by()%>%
        dplyr::filter(FDR_GoMF<=FDR_GoMFTopn) %>%
        dplyr::group_by(MF_Name) %>%
        dplyr::arrange(desc(abs_score))%>%
        dplyr::top_n(n=maxGene, wt = abs_score)%>%
        dplyr::group_by(MF_Name,signChange)%>%
        dplyr::arrange(abs_score)%>%
        dplyr::mutate(rank=signChange*seq(1,n()))->selectedDf_simp
      

      selectedDf_simp %>%
        dplyr::group_by()%>%
        dplyr::filter(selectedDf_simp$MF_Name %in% as.vector(GoMF_result$MF_Name))%>%
        dplyr::group_by(MF_Name_short)%>%
        dplyr::summarise(FDR_GoMF=unique(FDR_GoMF))%>%
        dplyr::arrange(FDR_GoMF)->summary_df_simp
      
      as.vector(summary_df_simp$MF_Name_short)
      
      
      selectedDf_simp$MF_Name_short <- factor(selectedDf_simp$MF_Name_short, 
                                              levels = rev(as.vector(summary_df_simp$MF_Name_short)))
      ###*****************************
      # if(all(dataName$objectName$pick_data=="mrna",
      #        dataName$objectName$growthPhase_names=="Sta")){browser()}
      
      ###*****************************
      # Generate Figures 
      # a) Complex figure
      scaleHigh=max(abs(selectedDf$score_gene))
      scaleMid=0
      scaleLow=-max(abs(selectedDf$score_gene))
      
      fig01=ggplot( selectedDf, aes( x=rank,y=MF_Name_long)) +
        geom_tile(aes(fill=score_gene))+
        scale_fill_gradientn(colours=c("Blue","Grey50","Red"),
                             values=rescale(c(scaleLow,scaleMid,scaleHigh)),
                             limits=c(scaleLow,scaleHigh),
                             guide = guide_colorbar(title = "-sign(cor)*P_log10"))+
        geom_text(aes(label=ID),size=3, colour="White", fontface="bold")+
        theme_bw()+
        scale_x_continuous(breaks=min(selectedDf$rank):max(selectedDf$rank))+
        ggtitle(paste0(objectName,"_mf"))+
        theme(axis.line.y = element_blank(),
              legend.position="bottom",
              axis.title.y = element_blank(),
              panel.grid.minor=element_blank(),
              panel.grid.major.x=element_blank())
      
      #print(fig01)
      
      
      
      
      # b) simple figure with geom point
      GoMF_result_narrow %>%
        dplyr::arrange(FDR_GoMF)%>%
        dplyr::mutate(MF_Name_short=gsub("*.*~","",MF_Name))%>%
        dplyr::filter(MF_Name_short%in%as.vector(unique(selectedDf_simp$MF_Name_short)))->GoMF_result_narrow
      
      newLabels<-replace_fun(input_vector = gsub("*.*:","",rev(GoMF_result_narrow$MF_Name_short)),
                             initialVal = c("hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds, in linear amides",
                                            "Biosynthesis of siderophore group nonribosomal peptides",
                                            "di−, tri−valent inorganic cation transmembrane transporter activity"),
                             finalVal = c("hydrolase activity, acting on carbon-nitrogen\n (but not peptide) bonds, in linear amides",
                                          "Biosynthesis of siderophore group\n nonribosomal peptides",
                                          "di−, tri−valent inorganic cation\n transmembrane transporter activity"))
      
      
      minimumFold=min(selectedDf_simp$log2)
      if(minimumFold>-1){minimumFold=-1}
      maximumFold=max(selectedDf_simp$log2)
      if(maximumFold<1){maximumFold=1}
      
      
      fig03=ggplot(selectedDf_simp, aes( x=log2,y=MF_Name_short)) +
        geom_point(colour="blue", size=2.5)+
        geom_vline(xintercept = c(log2(1/2),log2(2)), colour="orange", linetype = "longdash")+
        geom_vline(xintercept = c(log2(1)), colour="black", linetype = "longdash")+
        geom_text_repel(aes(label=ID),size=3, colour="Black", fontface="plain")+
        theme_bw()+
        scale_x_continuous(breaks=seq(floor(minimumFold),ceiling(maximumFold)))+
        scale_y_discrete(labels=newLabels)+
        xlab("Log2 Fold Change")+
        theme(axis.line.y = element_blank(),
              legend.position="bottom",
              axis.title.y = element_blank(),
              panel.grid.minor=element_blank(),
              panel.grid.major.x=element_blank(),
              strip.text.x = element_text(size = 16),
              strip.text.y = element_text(size = 16),
              axis.text.x=element_text(size=10),
              axis.text.y=element_text(size=12),
              axis.title.x=element_text(size=16),
              axis.title.y=element_text(size=16),
              legend.title=element_text(size=14),
              legend.text=element_text(size=14))
      
      
      ###*****************************
      
      
      ###*****************************
      # Save Files
      selectedDf<-cbind(selectedDf, 
                        unique(GoMF_input_df[,c("pick_data","growthPhase","test_for","vs")]),
                        df_category="MF")
      write.csv(x = selectedDf, file = paste0("../d_results/",objectName,"_mf.csv"))
      ###*****************************
      
      
      ###*****************************
      # Save Figures
      
      # Detailed Figure
      colWidth=ifelse(max(selectedDf$rank)-min(selectedDf$rank)+1<16, 
                      16, 
                      max(selectedDf$rank)-min(selectedDf$rank))
      rowWidth=ifelse(nrow(summary_df)*1<3,3,nrow(summary_df)*1)
      cowplot::save_plot(filename = paste0("../d_figures/",objectName,"_mf.pdf"),
                         plot = fig01,
                         nrow=1.5,
                         ncol=2,
                         base_height = rowWidth,
                         base_width = colWidth,
                         limitsize = FALSE)
      
      # Save simple figure
      if(nrow(summary_df_simp)<5)
      {
        rowWidth=7
      }
      else if(nrow(summary_df_simp)<10)
      {
        rowWidth=9
      }
      else
      {
        rowWidth=nrow(summary_df_simp)*.9 
      }
      
      if(objectName=="genes_P0.05Fold2_mrna_trT_set02_StcNasMglMgh_S_baseMglowMg_baseNa_Sta_noFilter_p1Sf_noNorm__Mg_mM_Levels_mf"){browser()}
      cowplot::save_plot(filename = paste0("../d_figures/simple",objectName,"_mf.pdf"),
                         plot = fig03,
                         base_height = rowWidth,
                         ncol=1.4,
                         nrow=1.2,
                         limitsize = FALSE)
      
    }
  }
}
###*****************************	
# End of the loop






