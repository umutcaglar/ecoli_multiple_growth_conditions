# generation of cophenetic distance matrix

###*****************************
# INITIAL COMMANDS TO RESET THE SYSTEM
rm(list = ls())
if (is.integer(dev.list())){dev.off()}
cat("\014")
###*****************************


###*****************************
# Set seed
set.seed(1023)
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
library("dplyr")
library("tidyr")
###*****************************


###*****************************
#Load Functions
source("../a_code_dataPreperation_RNA&Protein/data_naming_functions.R")
source("distance_table_functions.R")
###*****************************


###*****************************
# Parameters

# Number of Runs
numRuns=2000
# the distance type (cophenetic / euclidean)
distanceType="cophenetic"
###*****************************


###*****************************
# Find the csv files that need to be imported
dataName=name_data(initialValue=c("treeData"), # can be c("genes0.05","genes_P0.05Fold2","resDf")
                   dataType = "mrna", # can be "rna", "mrna", "protein", "protein_wo_NA"
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
                   deSeqSfChoice="p1Sf", # can be "regSf", "p1Sf"
                   normalizationMethodChoice= "vst", # can be "vst", "rlog", "log10", "noNorm"
                   test_for = "noTest")  # works only if normalizationMethodChoice == noNorm
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource", "noTest")

dataNameDF=as.data.frame(dataName[1])
dataName=paste(dataNameDF,collapse = "_")
load(file = paste0("../b_results/",dataName,".RData"))

clusteringDataNameDF=dataNameDF
clusteringDataNameDF$objectName.initial="clustering"
clusteringDataName=paste(clusteringDataNameDF,collapse = "_")
###*****************************


###*****************************
# Genrate Distance Matrices
cophenetic_distance_df_real=as.data.frame(as.matrix(cophenetic(METree.reorder_x_d)))
real_distance_df=as.data.frame(as.matrix(dist_x))


if (distanceType=="euclidean"){distance_df=real_distance_df}
if (distanceType=="cophenetic"){distance_df=cophenetic_distance_df_real}

###*****************************


###*****************************
# Generate additional columns for condition summary
colNamesCDM=colnames(distance_df)
rowNamesCDM=row.names(distance_df)
if(any(rowNamesCDM!=colNamesCDM)){stop("names are not consistent")}

conditionSummary<-conditionSummary[match(colNamesCDM, conditionSummary$dataSet),]
#********************************


#********************************
# real data analyze
growth_real<-calculateClusterMeans(variableInput="growthPhase",
                                   metaInput=conditionSummary,
                                   mainDataFrame=distance_df)

carbon_real<-calculateClusterMeans(variableInput="carbonSource",
                                   metaInput=conditionSummary,
                                   mainDataFrame=distance_df)

Mg_real<-calculateClusterMeans(variableInput="Mg_mM_Levels",
                               metaInput=conditionSummary,
                               mainDataFrame=distance_df)

Na_real<-calculateClusterMeans(variableInput="Na_mM_Levels",
                               metaInput=conditionSummary,
                               mainDataFrame=distance_df)

batch_real<-calculateClusterMeans(variableInput="batchNumber",
                                  metaInput=conditionSummary,
                                  mainDataFrame=distance_df)
#********************************


#********************************
# test for growth phase
for (counter01 in 1:numRuns){
  print(paste0("growth  ",counter01," / ",numRuns))
  conditionSummary %>%
    dplyr::mutate(growthPhase_fake=sample(growthPhase))->conditionSummary_fake
  
  growth_fake<-calculateClusterMeans(variableInput="growthPhase_fake",
                                     metaInput=conditionSummary_fake,
                                     mainDataFrame=distance_df)
  if(counter01==1){growth_fakeList=growth_fake}
  if(counter01!=1){growth_fakeList=rbind(growth_fakeList,growth_fake)}
}

z_growth=(growth_real-colMeans(growth_fakeList))/apply(growth_fakeList, 2, sd)
z_growth<-reShapeZScores(variableInput="growthPhase",
                            metaInput=conditionSummary,
                            z_data=z_growth)


# test for carbon source
for (counter01 in 1:numRuns){
  print(paste0("carbon  ",counter01," / ",numRuns))
  conditionSummary %>%
    dplyr::mutate(carbonSource_fake=sample(carbonSource))->conditionSummary_fake
  
  carbon_fake<-calculateClusterMeans(variableInput="carbonSource_fake",
                                     metaInput=conditionSummary_fake,
                                     mainDataFrame=distance_df)
  if(counter01==1){carbon_fakeList=carbon_fake}
  if(counter01!=1){carbon_fakeList=rbind(carbon_fakeList,carbon_fake)}
}

z_carbon=(carbon_real-colMeans(carbon_fakeList))/apply(carbon_fakeList, 2, sd)
z_carbon<-reShapeZScores(variableInput="carbonSource",
                            metaInput=conditionSummary,
                            z_data=z_carbon)

# test for Mg_mM Levels
for (counter01 in 1:numRuns){
  print(paste0("Mg  ",counter01," / ",numRuns))
  conditionSummary %>%
    dplyr::mutate(Mg_mM_Levels_fake=sample(Mg_mM_Levels))->conditionSummary_fake
  
  Mg_fake<-calculateClusterMeans(variableInput="Mg_mM_Levels_fake",
                                 metaInput=conditionSummary_fake,
                                 mainDataFrame=distance_df)
  if(counter01==1){Mg_fakeList=Mg_fake}
  if(counter01!=1){Mg_fakeList=rbind(Mg_fakeList,Mg_fake)}
}

z_Mg=(Mg_real-colMeans(Mg_fakeList))/apply(Mg_fakeList, 2, sd)
z_Mg<-reShapeZScores(variableInput="Mg_mM_Levels",
                        metaInput=conditionSummary,
                        z_data=z_Mg)


# test for Na_mM Levels
for (counter01 in 1:numRuns){
  print(paste0("Na  ",counter01," / ",numRuns))
  conditionSummary %>%
    dplyr::mutate(Na_mM_Levels_fake=sample(Na_mM_Levels))->conditionSummary_fake
  
  Na_fake<-calculateClusterMeans(variableInput="Na_mM_Levels_fake",
                                 metaInput=conditionSummary_fake,
                                 mainDataFrame=distance_df)
  if(counter01==1){Na_fakeList=Na_fake}
  if(counter01!=1){Na_fakeList=rbind(Na_fakeList,Na_fake)}
}

z_Na=(Na_real-colMeans(Na_fakeList))/apply(Na_fakeList, 2, sd)
z_Na<-reShapeZScores(variableInput="Na_mM_Levels",
                        metaInput=conditionSummary,
                        z_data=z_Na)


# test for batch number
for (counter01 in 1:numRuns){
  print(paste0("Batch  ",counter01," / ",numRuns))
  conditionSummary %>%
    dplyr::mutate(batchNumber_fake=sample(batchNumber))->conditionSummary_fake
  
  batch_fake<-calculateClusterMeans(variableInput="batchNumber_fake",
                                    metaInput=conditionSummary_fake,
                                    mainDataFrame=distance_df)
  if(counter01==1){batch_fakeList=batch_fake}
  if(counter01!=1){batch_fakeList=rbind(batch_fakeList,batch_fake)}
}

z_batch=(batch_real-colMeans(batch_fakeList))/apply(batch_fakeList, 2, sd)
z_batch<-reShapeZScores(variableInput="batchNumber",
                           metaInput=conditionSummary,
                           z_data=z_batch)
#********************************


#********************************
# Combine results into a single DF
z_scores<-rbind(z_growth,z_carbon,z_Mg,z_Na,z_batch)
z_scores$Overall_Z_score<-sprintf("%.2f", z_scores$Overall_Z_score)
z_scores$Z_score<-sprintf("%.2f", z_scores$Z_score)
#********************************


#********************************
write.csv(x = z_scores,file = paste0("../b_results/",clusteringDataName,"_",distanceType,".csv"),
          row.names = F,quote = TRUE)
#********************************



