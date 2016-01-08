# Generate Heatmap by hand

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
require("reshape2")
library("lazyeval")
require("flashClust")

# For Plotting
library("ggplot2")
library("RColorBrewer")
library("grid")
library("gridExtra")
library("cowplot")
require("ggdendro")
require("scales")
require("gtable")
###*****************************


###*****************************
# PARAMETERS FOR RNA DATA
conditionNumberChoice="uniqueCondition"
dataTypeChoice="mrna" # can be "rna", "mrna","protein","protein_wo_NA","protein_wo_NAx6"
badDataFilterSetChoice="set02" # "set00", "set01", "set02"
dataTimeChoice="wholeSet" # exponential / stationary / late_stationary / wholeSet
MgLevelChoice="allMg" # allMg highMg midMg lowMg
NaLevelChoice="allNa" # allNa lowNa highNa
carbonTypeChoice="SYAN" # a letter combination from the list "SYAN"  
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
load(file = paste0("../a_results/",step03,".RData"))
assign(x = "mainDataFrameRNA",value = get(step03))
assign(x = "conditionRNA",value = get(conditionName))
###*****************************

###*****************************
# Modify Condition RNA
conditionRNA %>%
  dplyr::mutate(dataSet2d=ifelse(carbonSource=="glucose","S",NA)) %>%
  dplyr::mutate(dataSet2d=ifelse(carbonSource=="glycerol","Y",dataSet2d)) %>%
  dplyr::mutate(dataSet2d=ifelse(carbonSource=="lactate","A",dataSet2d)) %>%
  dplyr::mutate(dataSet2d=ifelse(carbonSource=="gluconate","N",dataSet2d))->conditionRNA

conditionRNA %>%
  dplyr::mutate(dataSet2a=ifelse(growthPhase=="exponential","exp",NA)) %>%
  dplyr::mutate(dataSet2a=ifelse(growthPhase=="stationary","sta",dataSet2a)) %>%
  dplyr::mutate(dataSet2a=ifelse(growthPhase=="late_stationary","lat",dataSet2a)) ->conditionRNA

conditionRNA %>%
  dplyr::mutate(dataSet2e=ifelse(experiment=="glucose_time_course","Stc",NA)) %>%
  dplyr::mutate(dataSet2e=ifelse(experiment=="glycerol_time_course","Ytc",dataSet2e)) %>%
  dplyr::mutate(dataSet2e=ifelse(experiment=="NaCl_stress","Nas",dataSet2e)) %>%
  dplyr::mutate(dataSet2e=ifelse(experiment=="lactate_growth","Agr",dataSet2e)) %>%
  dplyr::mutate(dataSet2e=ifelse(experiment=="gluconate_growth","Ngr",dataSet2e)) %>%
  dplyr::mutate(dataSet2e=ifelse(experiment=="MgSO4_stress_low","Mgl",dataSet2e)) %>%
  dplyr::mutate(dataSet2e=ifelse(experiment=="MgSO4_stress_high","Mgh",dataSet2e)) -> conditionRNA

conditionRNA %>%
  mutate(dataSet2=paste0(dataSet,"_",
                         dataSet2a,"_",
                         Mg_mM_Levels,"_",
                         Na_mM_Levels,"_",
                         dataSet2d,"_",
                         dataSet2e))->conditionRNA
###*****************************


###*****************************
# PARAMETERS FOR PROTEIN DATA
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
load(file = paste0("../a_results/",step03,".RData")) #Used for parameters data
assign(x = "mainDataFrameProtein",value = get(step03))
assign(x = "conditionProtein",value = get(conditionName))
###*****************************


###*****************************
# Calculate the average of repeated conditions
mainDataFrameProtein<-as.data.frame(mainDataFrameProtein)
mainDataFrameProtein %>% mutate(gene_id=row.names(mainDataFrameProtein))->mainDataFrameProteinTidy
endVal=ncol(mainDataFrameProtein)

mainDataFrameProteinTidy %>%
  tidyr::gather(conditionProtein,read,1:endVal)%>%
  dplyr::mutate(conditionShort=gsub("_[[:digit:]]$","",conditionProtein))->mainDataFrameProteinTidy

mainDataFrameProteinTidy %>% 
  group_by(gene_id,conditionShort)%>%
  dplyr::summarise(read=mean(read))->mainDataFrameProteinSummary

mainDataFrameProteinSummary %>%
  dplyr::group_by()%>%  
  tidyr::spread(key = conditionShort, value = read)->mainDataFrameProteinSummary

row.names(mainDataFrameProteinSummary)<-mainDataFrameProteinSummary$gene_id
mainDataFrameProteinSummary %>%
  dplyr::select(-gene_id)->mainDataFrameProteinSummary
mainDataFrameProtein=as.matrix(mainDataFrameProteinSummary)
###*****************************


###*****************************
#Reduce number of rows in condition and get rid of repeats
conditionProtein=conditionProtein[which(substr(as.vector(conditionProtein$dataSet),10,10) %in% c(0)),]
conditionProtein$dataSet<-gsub("_[[:digit:]]$","",as.vector(conditionProtein$dataSet))
conditionProtein$sampleNum<-gsub("_[[:digit:]]$","",as.vector(conditionProtein$sampleNum))
###*****************************


###*****************************
# Generate a new column to condition df that uniquely defines the condition named dataSet2
conditionProtein %>%
  dplyr::mutate(dataSet2d=ifelse(carbonSource=="glucose","S",NA)) %>%
  dplyr::mutate(dataSet2d=ifelse(carbonSource=="glycerol","Y",dataSet2d)) %>%
  dplyr::mutate(dataSet2d=ifelse(carbonSource=="lactate","A",dataSet2d)) %>%
  dplyr::mutate(dataSet2d=ifelse(carbonSource=="gluconate","N",dataSet2d))->conditionProtein

conditionProtein %>%
  dplyr::mutate(dataSet2a=ifelse(growthPhase=="exponential","exp",NA)) %>%
  dplyr::mutate(dataSet2a=ifelse(growthPhase=="stationary","sta",dataSet2a)) %>%
  dplyr::mutate(dataSet2a=ifelse(growthPhase=="late_stationary","lat",dataSet2a)) ->conditionProtein

conditionProtein %>%
  dplyr::mutate(dataSet2e=ifelse(experiment=="glucose_time_course","Stc",NA)) %>%
  dplyr::mutate(dataSet2e=ifelse(experiment=="glycerol_time_course","Ytc",dataSet2e)) %>%
  dplyr::mutate(dataSet2e=ifelse(experiment=="NaCl_stress","Nas",dataSet2e)) %>%
  dplyr::mutate(dataSet2e=ifelse(experiment=="lactate_growth","Agr",dataSet2e)) %>%
  dplyr::mutate(dataSet2e=ifelse(experiment=="gluconate_growth","Ngr",dataSet2e)) %>%
  dplyr::mutate(dataSet2e=ifelse(experiment=="MgSO4_stress_low","Mgl",dataSet2e)) %>%
  dplyr::mutate(dataSet2e=ifelse(experiment=="MgSO4_stress_high","Mgh",dataSet2e)) -> conditionProtein


# "Stc" for "glucose_time_course", 
# "Ytc" for "glycerol_time_course", 
# "Nas" for "NaCl_stress", 
# "Agr" for "lactate_growth", 
# "Ngr" for "gluconate_growth", 
# "Mgl" for "MgSO4_stress_low", 
# "Mgh" for "MgSO4_stress_high" 


conditionProtein %>%
  dplyr::mutate(dataSet2=paste0(dataSet,"_",
                         dataSet2a,"_",
                         Mg_mM_Levels,"_",
                         Na_mM_Levels,"_",
                         dataSet2d,"_",
                         dataSet2e))->conditionProtein
###*****************************


###*****************************
# Load Dictionaries
dictionaryProtein=read.csv(paste0("../generateDictionary/","nameDictionary_Protein.csv"))
dictionaryRNA=read.csv(paste0("../generateDictionary/","nameDictionary_RNA_barrick.csv")) 
###*****************************


###*****************************
# find intersected conditions
conditionProtein$sampleNum=as.integer(conditionProtein$sampleNum)
conditionRNA$sampleNum=as.integer(conditionRNA$sampleNum)
conditionIntersection<-dplyr::semi_join(conditionProtein,conditionRNA)
###*****************************


###*****************************
# Pick intersected conditions and generate a data tidy df for RNA
mainDataFrameRNA<-as.data.frame(mainDataFrameRNA)

colNums <- match(as.vector(conditionIntersection$dataSet),names(mainDataFrameRNA))
mainDataFrameRNA %>% 
  dplyr::select(colNums)%>%
  dplyr::mutate(gene_id=row.names(mainDataFrameRNA))->mainDataFrameRNATidy
endVal=ncol(mainDataFrameRNATidy)

mainDataFrameRNATidy %>%
  tidyr::gather(conditionRNA,read,1:(endVal-1))->mainDataFrameRNATidy
###*****************************

###*****************************
# add dictionary to RNA
mainDataFrameRNATidy<-dplyr::left_join(mainDataFrameRNATidy,dictionaryRNA)
colnames(mainDataFrameRNATidy)[3]<-"readRNA"
colnames(mainDataFrameRNATidy)[2]<-"dataSet"
###*****************************

###*****************************
# Pick intersected conditions and generate a data tidy df for Proteins
mainDataFrameProtein<-as.data.frame(mainDataFrameProtein)

colNums <- match(as.vector(conditionIntersection$dataSet),names(mainDataFrameProtein))
mainDataFrameProtein %>% 
  dplyr::select(colNums)%>%
  dplyr::mutate(gene_id=row.names(mainDataFrameProtein))->mainDataFrameProteinTidy
endVal=ncol(mainDataFrameProteinTidy)

mainDataFrameProteinTidy %>%
  tidyr::gather(conditionProtein,read,1:(endVal-1))->mainDataFrameProteinTidy
###*****************************


###*****************************
# add dictionary to Protein
mainDataFrameProteinTidy %>%
  dplyr::left_join(x=., y = dictionaryProtein, by = c("gene_id" = "Protein_id"))->mainDataFrameProteinTidy
colnames(mainDataFrameProteinTidy)[3]<-"readProtein"
colnames(mainDataFrameProteinTidy)[2]<-"dataSet"
colnames(mainDataFrameProteinTidy)[1]<-"protein_id"
###*****************************


###*****************************
# Combine Data
dplyr::inner_join(x=mainDataFrameRNATidy,
                 y=mainDataFrameProteinTidy)->mainDataFrameTidy
###*****************************


###*****************************
# Calculate Correlation
mainDataFrameTidy %>%
  dplyr::group_by(dataSet) %>%
  dplyr::summarize(correlation=cor(x=readRNA, y=readProtein))->summaryResult

dplyr::left_join(summaryResult, conditionIntersection)->summaryResult
 
ggplot(summaryResult, aes(x = correlation)) + 
  geom_density(fill="blue") +
  theme_classic()+
  xlim(0,1)+
  theme( axis.text.x=element_text(size=16),
         axis.text.y=element_text(size=16),
         axis.title.x=element_text(size=16),
         axis.title.y=element_text(size=16),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14))

mean(summaryResult$correlation)

###*****************************
# Calculate Correlation Stationary
summaryResult %>%
  dplyr::filter(growthPhase=="stationary")->summaryResultStationary

ggplot(summaryResultStationary, aes(x = correlation)) + 
  geom_density(fill="blue") +
  theme_classic()+
  xlim(0,1)+
  theme( axis.text.x=element_text(size=16),
         axis.text.y=element_text(size=16),
         axis.title.x=element_text(size=16),
         axis.title.y=element_text(size=16),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14))

###*****************************
# Calculate Correlation Exponential
summaryResult %>%
  dplyr::filter(growthPhase=="exponential")->summaryResultExponential

ggplot(summaryResultExponential, aes(x = correlation)) + 
  geom_density(fill="blue") +
  theme_classic()+
  xlim(0,1)+
  theme( axis.text.x=element_text(size=16),
         axis.text.y=element_text(size=16),
         axis.title.x=element_text(size=16),
         axis.title.y=element_text(size=16),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14))

mean(summaryResultExponential$correlation)

