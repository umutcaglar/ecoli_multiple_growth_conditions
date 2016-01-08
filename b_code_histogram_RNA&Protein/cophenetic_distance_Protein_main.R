# generation of cophenetic distance matrix

###*****************************
# INITIAL COMMANDS TO RESET THE SYSTEM

rm(list = ls())
if (is.integer(dev.list())){dev.off()}
cat("\014")

# Set Working Directory
setwd('/Users/umut/GitHub/ecoli_multiple_growth_conditions/b_code_histogram_RNA&Protein/') # mac computer
###*****************************


###*****************************
# INSTALL LIBRARIES
library("dplyr")
library("tidyr")
###*****************************


###*****************************
# PARAMETERS FOR DATA
conditionNumberChoice="uniqueCondition"
dataTypeChoice="protein_wo_NA" # can be "rna", "mrna","protein","protein_wo_NA","protein_wo_NAx6"
badDataFilterSetChoice="set00" # "set00", "set01", "set02"
dataTimeChoice="wholeSet" # exponential / stationary / late_stationary / wholeSet
MgLevelChoice="allMg" # allMg highMg midMg lowMg
NaLevelChoice="allNa" # allNa lowNa highNa
carbonTypeChoice="SY" # a letter combination from the list "SYAN"  
#S (glucose), Y (glycerol), A (lactate), N (gluconate)
filterTypeChoice="noFilter" # mean/sd/ max/ noFilter
deSeqNormChoice="p1"
normalizationMethodChoice="vst" # "vst" , "log2" 
experimentChoice=c("allEx") # can be "allEx", 
# or a combination of below
# "Stc" for "glucose_time_course", 
# "Ytc" for "glycerol_time_course", 
# "Nas" for "NaCl_stress", 
# "Agr" for "lactate_growth", 
# "Ngr" for "gluconate_growth", 
# "Mgh" for "MgSO4_stress_low", 
# "Mgl" for "MgSO4_stress_high"
###*****************************


###*****************************
# Order Experiment List
experimentListOrder=c("allEx","Stc","Ytc","Nas","Agr","Ngr","Mgl","Mgh")
currentOrder=experimentChoice
experimentChoice=currentOrder[order(match(currentOrder,experimentListOrder))]
#*******************************


###*****************************
# FILE NAME GENERATE

experimentList=paste(experimentChoice, collapse = '')
longNameList=paste0(dataTypeChoice,"_",dataTimeChoice,
                    "_",MgLevelChoice,"_",NaLevelChoice,
                    "_",carbonTypeChoice,"_",badDataFilterSetChoice,
                    "_",experimentList)
step02=paste0("unnormalized_",filterTypeChoice,"_",longNameList)
step03=paste0("normalized_",normalizationMethodChoice,"_",deSeqNormChoice,"_",
              filterTypeChoice,"_",longNameList) 
conditionName=paste0("condition_",normalizationMethodChoice,"_",deSeqNormChoice,"_",
                     filterTypeChoice,"_",longNameList) 
###*****************************



###*****************************
# Load data
load(file = paste0("../b_results/treeFile_",step03,".RData"))
tidyr::gather(conditionSummary, category, condition, growthPhase:batchNumber)->conditionSummary_tidy
###*****************************


###*****************************
# Install Functions
source("cophenetic_distance_Protein_functions.R")
###*****************************


###*****************************
# Number of Runs
numRuns=100
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
new_name_all_results=paste0("all_results_",step03)
assign(x=new_name_all_results,all_results)

new_name_all_results_Summary=paste0("all_results_Summary_",step03)
assign(x=new_name_all_results_Summary,all_results_Summary)
#******************************


#******************************
save(list = c(new_name_all_results,new_name_all_results_Summary),
     file = paste0("../b_results/clustering_",step03,".RData"))

write.csv(x = all_results_Summary, 
          file = paste0("../b_results/clustering_",step03,".csv"),
          row.names = F)
#******************************