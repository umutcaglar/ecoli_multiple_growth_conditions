# generation of cophenetic distance matrix

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
library("dplyr")
library("tidyr")
###*****************************


###*****************************
#Load Functions
source("../a_code_dataPreperation_RNA&Protein/data_naming_functions.R")
source("cophenetic_distance_functions.R")
###*****************************

###*****************************
# Find the csv files that need to be imported
dataName=name_data(initialValue=c("treeData"), # can be c("genes0.05","genes_P0.05Fold2","resDf")
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

clusteringDataName=dataName
clusteringDataName$objectName.initial="clustering"

dataName=paste(dataName,collapse = "_")
clusteringDataName=paste(clusteringDataName,collapse = "_")

load(file = paste0("../b_results/",dataName,".RData"))
###*****************************


###*****************************
# Make condition summary tidy
tidyr::gather(conditionSummary, category, condition, growthPhase:batchNumber)->conditionSummary_tidy
###*****************************
browser()

###*****************************
# Number of Runs
numRuns=1000
###*****************************


###*****************************
# Calculate Cophenetic Distance
cophenetic_distance_matrix=as.matrix(cophenetic(METree.reorder_x_d))
distance_matrix=as.matrix(dist_x)
###*****************************


###*****************************
# Correlation in between distance matrix and cophenetic distance matrix
cor(t(rbind(as.vector(cophenetic_distance_matrix),as.vector(distance_matrix))))
###*****************************


###*****************************
# calculation of median 

# real values
main_growthPhase=mainDistanceCategory(cophenetic_distance_matrix, 
                                      conditionSummary_tidy, 
                                      "growthPhase")

main_carbonSource=mainDistanceCategory(cophenetic_distance_matrix, 
                                       conditionSummary_tidy, 
                                       "carbonSource")

main_Mg=mainDistanceCategory(cophenetic_distance_matrix, 
                             conditionSummary_tidy, 
                             "Mg_mM_Levels")

main_Na=mainDistanceCategory(cophenetic_distance_matrix, 
                             conditionSummary_tidy, 
                             "Na_mM_Levels")

main_batch=mainDistanceCategory(cophenetic_distance_matrix, 
                             conditionSummary_tidy, 
                             "batchNumber")
#******************************


#******************************
# calculation of fakes
for (counter01 in 1:numRuns){
  print(counter01)
  main_Fake=mainDistanceCategoryF(cophenetic_distance_matrix, 
                                          conditionSummary_tidy, 
                                          "growthPhase")
  main_Fake=dplyr::mutate(main_Fake,iteration=counter01)
  if(counter01==1){main_FakeList=main_Fake}
  if(counter01!=1){main_FakeList=rbind(main_FakeList,main_Fake)}
}
main_growthPhaseF=main_FakeList


for (counter01 in 1:numRuns){
  print(counter01)
  main_Fake=mainDistanceCategoryF(cophenetic_distance_matrix, 
                                  conditionSummary_tidy, 
                                  "carbonSource")
  main_Fake=dplyr::mutate(main_Fake,iteration=counter01)
  if(counter01==1){main_FakeList=main_Fake}
  if(counter01!=1){main_FakeList=rbind(main_FakeList,main_Fake)}
}
main_carbonSourceF=main_FakeList

for (counter01 in 1:numRuns){
  print(counter01)
  main_Fake=mainDistanceCategoryF(cophenetic_distance_matrix, 
                                  conditionSummary_tidy, 
                                  "Mg_mM_Levels")
  main_Fake=dplyr::mutate(main_Fake,iteration=counter01)
  if(counter01==1){main_FakeList=main_Fake}
  if(counter01!=1){main_FakeList=rbind(main_FakeList,main_Fake)}
}
main_MgF=main_FakeList

for (counter01 in 1:numRuns){
  print(counter01)
  main_Fake=mainDistanceCategoryF(cophenetic_distance_matrix, 
                                  conditionSummary_tidy, 
                                  "Na_mM_Levels")
  main_Fake=dplyr::mutate(main_Fake,iteration=counter01)
  if(counter01==1){main_FakeList=main_Fake}
  if(counter01!=1){main_FakeList=rbind(main_FakeList,main_Fake)}
}
main_NaF=main_FakeList

for (counter01 in 1:numRuns){
  print(counter01)
  main_Fake=mainDistanceCategoryF(cophenetic_distance_matrix, 
                                  conditionSummary_tidy, 
                                  "batchNumber")
  main_Fake=dplyr::mutate(main_Fake,iteration=counter01)
  if(counter01==1){main_FakeList=main_Fake}
  if(counter01!=1){main_FakeList=rbind(main_FakeList,main_Fake)}
}
main_batchF=main_FakeList
#******************************


#******************************
# Calculation of z scores
main_growthPhaseF %>%
  dplyr::mutate(overall_fakeMean=mean(meanVal),
                overall_fakeStd=sd(meanVal))%>%
  dplyr::group_by(condition) %>%
  dplyr::summarize(fakeMean=mean(meanVal),
                   fakeStd=sd(meanVal),
                   overall_fakeMean=unique(overall_fakeMean),
                   overall_fakeStd=unique(overall_fakeStd),
                   numRepeats=length(condition),
                   category=unique(category))->FakeSummary
dplyr::left_join(main_growthPhase,FakeSummary)->main_growthPhase
main_growthPhase%>%
  dplyr::mutate(z_score=(meanVal-fakeMean)/fakeStd,
                overall_z_score=(overal_Mean-overall_fakeMean)/overall_fakeStd)->main_growthPhase

main_carbonSourceF %>%
  dplyr::mutate(overall_fakeMean=mean(meanVal),
                overall_fakeStd=sd(meanVal))%>%
  dplyr::group_by(condition) %>%
  dplyr::summarize(fakeMean=mean(meanVal),
                   fakeStd=sd(meanVal),
                   overall_fakeMean=unique(overall_fakeMean),
                   overall_fakeStd=unique(overall_fakeStd),
                   numRepeats=length(condition),
                   category=unique(category))->FakeSummary
dplyr::left_join(main_carbonSource,FakeSummary)->main_carbonSource
main_carbonSource%>%
  dplyr::mutate(z_score=(meanVal-fakeMean)/fakeStd,
                overall_z_score=(overal_Mean-overall_fakeMean)/overall_fakeStd)->main_carbonSource

main_MgF %>%
  dplyr::mutate(overall_fakeMean=mean(meanVal),
                overall_fakeStd=sd(meanVal))%>%
  dplyr::group_by(condition) %>%
  dplyr::summarize(fakeMean=mean(meanVal),
                   fakeStd=sd(meanVal),
                   overall_fakeMean=unique(overall_fakeMean),
                   overall_fakeStd=unique(overall_fakeStd),
                   numRepeats=length(condition),
                   category=unique(category))->FakeSummary
dplyr::left_join(main_Mg,FakeSummary)->main_Mg
main_Mg%>%
  dplyr::mutate(z_score=(meanVal-fakeMean)/fakeStd,
                overall_z_score=(overal_Mean-overall_fakeMean)/overall_fakeStd)->main_Mg

main_NaF %>%
  dplyr::mutate(overall_fakeMean=mean(meanVal),
                overall_fakeStd=sd(meanVal))%>%
  dplyr::group_by(condition) %>%
  dplyr::summarize(fakeMean=mean(meanVal),
                   fakeStd=sd(meanVal),
                   overall_fakeMean=unique(overall_fakeMean),
                   overall_fakeStd=unique(overall_fakeStd),
                   numRepeats=length(condition),
                   category=unique(category))->FakeSummary
dplyr::left_join(main_Na,FakeSummary)->main_Na
main_Na%>%
  dplyr::mutate(z_score=(meanVal-fakeMean)/fakeStd,
                overall_z_score=(overal_Mean-overall_fakeMean)/overall_fakeStd)->main_Na

main_batchF %>%
  dplyr::mutate(overall_fakeMean=mean(meanVal),
                overall_fakeStd=sd(meanVal))%>%
  dplyr::group_by(condition) %>%
  dplyr::summarize(fakeMean=mean(meanVal),
                   fakeStd=sd(meanVal),
                   overall_fakeMean=unique(overall_fakeMean),
                   overall_fakeStd=unique(overall_fakeStd),
                   numRepeats=length(condition),
                   category=unique(category))->FakeSummary
dplyr::left_join(main_batch,FakeSummary)->main_batch
main_batch%>%
  dplyr::mutate(z_score=(meanVal-fakeMean)/fakeStd,
                overall_z_score=(overal_Mean-overall_fakeMean)/overall_fakeStd)->main_batch

all_results=rbind_all(list(main_growthPhase,main_carbonSource,main_Mg,main_Na,main_batch))

all_results %>%
  dplyr::select(category,condition,z_score, overall_z_score)->all_results_Summary
#******************************


#******************************
new_name_all_results=paste0("all_results_",dataName)
assign(x=new_name_all_results,all_results)

new_name_all_results_Summary=paste0("all_results_Summary_",dataName)
assign(x=new_name_all_results_Summary,all_results_Summary)
#******************************


#******************************
save(list = c(new_name_all_results,new_name_all_results_Summary),
     file = paste0("../b_results/",clusteringDataName,".RData"))

write.csv(x = all_results_Summary, 
          file = paste0("../b_results/",clusteringDataName,".csv"),
          row.names = F)
#******************************
