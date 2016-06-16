# Draw heat map by hand for proteins
# Just an internal check for not metioned proteins in different csv files. 
# They will be assumed as 0 in further analyzeses

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
              "b_code_histogram_RNA&Protein/"))} # mac computer
###*****************************


###*****************************
# INSTALL LIBRARIES
require("dplyr")
require("tidyr")

# For Plotting
require("ggplot2")
require("cowplot")


###*****************************


###*****************************
#Load Functions
source("../a_code_dataPreperation_RNA&Protein/data_naming_functions.R")
source("replace_fun.R")
###*****************************

###*****************************
# Find the csv files that need to be imported
dataName=name_data(initialValue=c("resDf"), # can be c("genes0.05","genes_P0.05Fold2","resDf")
                   dataType = "protein", # can be "rna", "mrna", "protein", "protein_wo_NA"
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
                   growthPhaseVector = c("allPhase"), # can be "exponential","stationary","late_stationary" // "allPhase"
                   filterGenes = "noFilter", # can be "noFilter", "meanFilter", "maxFilter", "sdFilter" 
                   threshold=NA, # the threshold value for "meanFilter", "maxFilter", "sdFilter"
                   roundData=TRUE,
                   sumTechnicalReplicates=TRUE,
                   deSeqSfChoice="p1Sf", # can be "regSf", "p1Sf"
                   normalizationMethodChoice= "vst", # can be "vst", "rlog", "log10", "noNorm"
                   test_for = "noTest")  # works only if normalizationMethodChoice == noNorm
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource", "noTest")

dataName=as.data.frame(dataName[1])

metaDataName=dataName
metaDataName$objectName.initial="metaData"
treeDataName=dataName
treeDataName$objectName.initial="treeData"
heatMapName=dataName
heatMapName$objectName.initial="heatMap"

dataName=paste(dataName,collapse = "_")
metaDataName=paste(metaDataName,collapse = "_")
treeDataName=paste(treeDataName,collapse = "_")
heatMapName=paste(heatMapName,collapse = "_")

mainDataFrame=read.csv(file = paste0("../a_results/proteinMatrix_w_NA.csv"),header = TRUE,row.names = 1)
condition=read.csv(file = paste0("../a_results/",metaDataName,".csv"),header = TRUE)
###*****************************

###*****************************
# Make mainDataFrame binary
mainDataFrame[!is.na(mainDataFrame)]<-1
mainDataFrame[is.na(mainDataFrame)]<-0
###*****************************


###*****************************
# find not existed genes
mainDataFrame%>%
  dplyr::add_rownames(var = "gene_id")%>%
  dplyr::filter(MURI_089_0==0)%>%
  dplyr::select(gene_id)%>%
  t(.)%>%
  as.vector(.)->lactate

mainDataFrame%>%
  dplyr::add_rownames(var = "gene_id")%>%
  dplyr::filter(MURI_095_0==0)%>%
  dplyr::select(gene_id)%>%
  t(.)%>%
  as.vector(.)->gluconate

mainDataFrame%>%
  dplyr::add_rownames(var = "gene_id")%>%
  dplyr::filter(MURI_134_0==0)%>%
  dplyr::select(gene_id)%>%
  t(.)%>%
  as.vector(.)->Mg

mainDataFrame%>%
  dplyr::add_rownames(var = "gene_id")%>%
  dplyr::filter(MURI_098_0==0)%>%
  dplyr::select(gene_id)%>%
  t(.)%>%
  as.vector(.)->GlucoseTimeCourse


mainDataFrame%>%
  dplyr::add_rownames(var = "gene_id")%>%
  dplyr::filter(MURI_079_0==0)%>%
  dplyr::select(gene_id)%>%
  t(.)%>%
  as.vector(.)->Na

mainDataFrame%>%
  dplyr::add_rownames(var = "gene_id")%>%
  dplyr::filter(MURI_053_0==0)%>%
  dplyr::select(gene_id)%>%
  t(.)%>%
  as.vector(.)->GlycerolTimeCourse


maxLenght=max(length(lactate),
              length(gluconate),
              length(Mg),
              length(Na),
              length(GlucoseTimeCourse),
              length(GlycerolTimeCourse))

length(lactate)<-maxLenght
length(gluconate)<-maxLenght
length(Mg)<-maxLenght
length(Na)<-maxLenght
length(GlucoseTimeCourse)<-maxLenght
length(GlycerolTimeCourse)<-maxLenght

cbind(lactate,gluconate,Mg,Na,GlucoseTimeCourse,GlycerolTimeCourse)->list
write.csv(x = list,file = "../a_results/missingProteins.csv")
###*****************************


###*****************************
# Draw Heatmap
heatmap(x=as.matrix(mainDataFrame), scale = "none", col = c("red","blue"), main = "Gene existance heatmap",xlab = "samples", ylab="genes",cexRow=0.5,cexCol = 0.5)
###*****************************


