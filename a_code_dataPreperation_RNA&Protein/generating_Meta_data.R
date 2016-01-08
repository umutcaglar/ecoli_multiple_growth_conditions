# GENERATION META_DATA FILE 
# this is the main start file for the pipe line initial paper.
# Outputs of the file will be meta data file in form of CSV.
# The data frames for meta file for representing information about each column of data
#***************************************


#***************************************
# Initial Command to Reset the System & set working directory
rm(list = ls())
if (is.integer(dev.list())){dev.off()}
cat("\014")

# Set Working Directory
setwd('/Users/umut/GitHub/ecoli_multiple_growth_conditions/a_code_dataPreperation_RNA&Protein/') # mac computer
#***************************************

#***************************************
#Load Libraries and functions
require(dplyr)
require("grid")
require("gridExtra")
library("ggplot2")

#order function
source("fun_order.R")
#***************************************


# *******************************
# GENERATING META DATA FILE 


# Downloading Meta File
locationOfRNADataSets="../RNA_reads_03_04_2015/"
temp<-read.csv(file=paste0(locationOfRNADataSets,"sampleInformation.csv"), 
         header = TRUE, 
         sep = ",",
         dec = ".", fill = TRUE)

temp -> meta_data
remove(temp)

# add 2 columns indicating if we have data for RNA AND PROTEIN

# learn the list of rna files we have
locationOfRNADataSets="../RNA_reads_03_04_2015/"
RNA_FileList=dir(path=locationOfRNADataSets)

RawRNA_FileList=grep("raw_rna_count", RNA_FileList, value = TRUE)
RawRNA_FileList<-RawRNA_FileList[-grep("ND", RawRNA_FileList)] #get rid of not depleted 
RawRNA_colNames=gsub("_raw_rna_count.txt","",RawRNA_FileList)
RawRNA_colNames=sprintf("MURI_%03d",
                        as.numeric(gsub("MURI_","",RawRNA_colNames)))

RawRNA_colNames=sort(as.numeric(gsub("MURI_","",RawRNA_colNames)))
RawRNA_Repeats=as.data.frame(table(RawRNA_colNames))
colnames(RawRNA_Repeats)<-c("Sample.","RNA_Data_Freq")

# learn the list of Protein files we have
locationOfProteinDataSets="../Protein_reads_10_28_2015/"
Protein_FileList=dir(path=locationOfProteinDataSets)

RawProtein_FileList=grep("raw_protein_count", Protein_FileList, value = TRUE)
RawProtein_colNames=gsub("_raw_protein_count.txt","",RawProtein_FileList)
RawProtein_colNames=sort(as.numeric(gsub("_1|_0","",gsub("MURI_","",RawProtein_colNames))))

RawProtein_Repeats=as.data.frame(table(RawProtein_colNames))
colnames(RawProtein_Repeats)<-c("Sample.","Protein_Data_Freq")

merge(meta_data,RawRNA_Repeats,all.x = T)->meta_data
merge(meta_data,RawProtein_Repeats,all.x = T)->meta_data
meta_data$RNA_Data_Freq[is.na(meta_data$RNA_Data_Freq)] <- 0
meta_data$Protein_Data_Freq[is.na(meta_data$Protein_Data_Freq)] <- 0

# artificially make RNA data for MURI_083 "0"
meta_data$RNA_Data_Freq[which(meta_data$Sample.==83)]<-0
# artificially make Protein data for MURI_082 "0"
meta_data$Protein_Data_Freq[which(meta_data$Sample.==82)]<-0

# Downloading Growth Conditions
temp<-read.csv(file=paste0(locationOfRNADataSets,"AG3C_doubling_times.csv"), 
               header = TRUE, 
               sep = ",",
               dec = ".", fill = TRUE)

temp -> doubling_times
remove(temp)

# Extending doubling times
doubling_times %>% 
  mutate(experiment=ifelse(name=="Gluconate.tab","gluconate_growth",NA)) %>%
  mutate(experiment=ifelse(name=="Glucose.tab","glucose_time_course",experiment)) %>%
  mutate(experiment=ifelse(name=="Glycerol.tab","glycerol_time_course",experiment)) %>%
  mutate(experiment=ifelse(name=="Lactate.tab","lactate_growth",experiment)) %>%
  mutate(experiment=ifelse(name=="MgSO4_000.080_mM.tab","MgSO4_stress_high",experiment)) %>%
  mutate(experiment=ifelse(name=="MgSO4_000.800_mM.tab","MgSO4_stress_high",experiment)) %>%
  mutate(experiment=ifelse(name=="MgSO4_008.000_mM.tab","MgSO4_stress_high",experiment)) %>%
  mutate(experiment=ifelse(name=="MgSO4_050.000_mM.tab","MgSO4_stress_high",experiment)) %>%
  mutate(experiment=ifelse(name=="MgSO4_200.000_mM.tab","MgSO4_stress_high",experiment)) %>%
  mutate(experiment=ifelse(name=="MgSO4_400.000_mM.tab","MgSO4_stress_high",experiment)) %>%
  mutate(experiment=ifelse(name=="MgSO4-2_000.005_mM.tab","MgSO4_stress_low",experiment)) %>%
  mutate(experiment=ifelse(name=="MgSO4-2_000.010_mM.tab","MgSO4_stress_low",experiment)) %>%
  mutate(experiment=ifelse(name=="MgSO4-2_000.020_mM.tab","MgSO4_stress_low",experiment)) %>%
  mutate(experiment=ifelse(name=="MgSO4-2_000.040_mM.tab","MgSO4_stress_low",experiment)) %>%
  mutate(experiment=ifelse(name=="MgSO4-2_000.080_mM.tab","MgSO4_stress_low",experiment)) %>%
  mutate(experiment=ifelse(name=="NaCl_005_mM.tab","NaCl_stress",experiment)) %>%
  mutate(experiment=ifelse(name=="NaCl_100_mM.tab","NaCl_stress",experiment)) %>%
  mutate(experiment=ifelse(name=="NaCl_200_mM.tab","NaCl_stress",experiment)) %>%
  mutate(experiment=ifelse(name=="NaCl_300_mM.tab","NaCl_stress",experiment)) -> doubling_times

'%!in%' <- function(x,y)!('%in%'(x,y))
doubling_times %>% 
  mutate(Mg_mM=ifelse(experiment %!in% c("MgSO4_stress_low","MgSO4_stress_high"),0.8,NA)) %>%
  mutate(Mg_mM=ifelse(experiment=="MgSO4_stress_high",gsub("_mM.tab","",gsub("MgSO4_","",.$name)),Mg_mM)) %>%
  mutate(Mg_mM=ifelse(experiment=="MgSO4_stress_low",gsub("_mM.tab","",gsub("MgSO4-2_","",.$name)),Mg_mM)) %>%
  mutate(Mg_mM=as.numeric(Mg_mM))%>% 
  mutate(Na_mM=ifelse(experiment != "NaCl_stress",5,NA)) %>%
  mutate(Na_mM=ifelse(experiment=="NaCl_stress",gsub("_mM.tab","",gsub("NaCl_","",.$name)),Na_mM)) %>%
  mutate(Na_mM=as.numeric(Na_mM)) %>%
  select(-name)->doubling_times

names(doubling_times)[which(names(doubling_times)=="doubling.time.minutes")]<-"doublingTimeMinutes"
names(doubling_times)[which(names(doubling_times)=="doubling.time.minutes.95m")]<-"doublingTimeMinutes.95m"
names(doubling_times)[which(names(doubling_times)=="doubling.time.minutes.95p")]<-"doublingTimeMinutes_95p"
names(doubling_times)[which(names(doubling_times)=="r.squared")]<-"rSquared"

# EXTENDING meta_data

# addition of "MURI_" column with 3 decimal digits
meta_data %>% dplyr::mutate(dataSet=sprintf("MURI_%03d",Sample.))->meta_data
temp<-moveme(names(meta_data), "dataSet first")
meta_data <- meta_data[temp]

# generating SampleNum column 
meta_data %>% mutate(sampleNum=as.numeric(Sample.))%>% select(-Sample.)->meta_data
temp<-moveme(names(meta_data), "sampleNum first")
meta_data <- meta_data[temp]

# generate batchNumber
meta_data %>%
  mutate(batchNumber=as.character(substr(.[["Batch"]],1,3)))->meta_data

# generate experiment column
names(meta_data)[which(names(meta_data)=="Experiment")]<-"experiment"
meta_data %>% mutate(experiment=ifelse(experiment=="Lactate_growth",
                                      "lactate_growth", as.character(experiment)))->meta_data
meta_data %>% mutate(experiment=ifelse(experiment=="Gluconate_growth",
                                      "gluconate_growth", as.character(experiment)))->meta_data
meta_data %>% mutate(experiment=ifelse(experiment=="MgSO4_stress" & batchNumber %in% c("022","023","024"),
                                      "MgSO4_stress_high", as.character(experiment)))->meta_data
meta_data %>% mutate(experiment=ifelse(experiment=="MgSO4_stress" & batchNumber %in% c("025","026","027"),
                                      "MgSO4_stress_low", as.character(experiment)))->meta_data

# generate carbonSource Column
meta_data %>%
  mutate(carbonSource=GrowthConditions)->meta_data
meta_data$carbonSource=gsub("\\+NaCl","",meta_data$carbonSource)
meta_data$carbonSource=gsub("\\+MgSO4","",meta_data$carbonSource)

# generate MgSO4_mM Column
meta_data %>% 
  mutate(Mg_mM=ifelse(grepl("[E,S][0-9]",.[["Batch"]]) & experiment %in% c("MgSO4_stress_high","MgSO4_stress_low"),
                              gsub("[0-9][0-9][0-9][E,S]","",as.character(.[["Batch"]])),
                              NA))%>%
  mutate(Mg_mM=ifelse(is.na(Mg_mM),0.8,as.numeric(Mg_mM))) %>%
  mutate(Mg_mM=ifelse(batchNumber %in% c("025","026","027"),Mg_mM/1000,Mg_mM))-> meta_data

# generate NaCl_mM
meta_data %>% 
  mutate(Na_mM=ifelse(grepl("[E,S][0-9]",.[["Batch"]]) & experiment=="NaCl_stress",
                         gsub("[0-9][0-9][0-9][E,S]","",as.character(.[["Batch"]])),
                         NA)) %>%
  mutate(Na_mM=ifelse(is.na(Na_mM),5,as.numeric(Na_mM))) -> meta_data


# generate concentration column
meta_data %>% 
  group_by() %>% 
  mutate(concentration=ifelse(grepl("[E,S][0-9]",.[["Batch"]]),
                              gsub("[0-9][0-9][0-9][E,S]","",as.character(.[["Batch"]])),
                              NA)) %>%
  mutate(concentration=ifelse(is.na(concentration),NA,as.numeric(concentration)))-> meta_data

# generate growthTime_hr column
names(meta_data)[which(names(meta_data)=="Growthtime.hr.")]<-"growthTime_hr"

# generate growthPhase column lag/ exponential / stationary / late stationary
meta_data %>% 
  mutate(growthPhase=ifelse(grepl("E",.[["Batch"]]),"exponential",NA)) %>%
  mutate(growthPhase=ifelse(grepl("S",.[["Batch"]]),"stationary",growthPhase)) %>%
  mutate(growthPhase=ifelse(grepl("[t,R]",.[["Batch"]]) & growthTime_hr<=14,
                            "exponential",growthPhase))%>%
  mutate(growthPhase=ifelse(grepl("[t,R]",.[["Batch"]]) & growthTime_hr>14 & growthTime_hr<=48,
                            "stationary",growthPhase))%>%
  mutate(growthPhase=ifelse(grepl("[t,R]",.[["Batch"]]) & growthTime_hr>48,
                            "late_stationary",growthPhase))->meta_data

# GENERATE Mg_mM_Levels & Na_mM_Levels
meta_data %>%
dplyr::mutate(Mg_mM_Levels = ifelse(Mg_mM<0.8, "lowMg",NA),
                Mg_mM_Levels = ifelse(Mg_mM==0.8, "baseMg",Mg_mM_Levels),
                Mg_mM_Levels = ifelse(Mg_mM>0.8, "highMg",Mg_mM_Levels))->meta_data

meta_data %>%
  dplyr::mutate(Na_mM_Levels = ifelse(Na_mM==5, "baseNa",NA),
                Na_mM_Levels = ifelse(Na_mM>5, "highNa",Na_mM_Levels))->meta_data

# Generate Unique Condition Column
meta_data %>%
  group_by(growthPhase,carbonSource,Mg_mM,Na_mM) %>%
  summarize(replicaNum=length(growthPhase)) ->sampleSizeDf

sampleSizeDf%>%
  group_by()%>%
  mutate(uniqueCondition=sprintf("unique_condition_%02d",as.numeric(rownames(.)))) -> sampleSizeDf


sampleSizeDf%>%
  select(-replicaNum)%>%
  left_join(meta_data,.)->meta_data

# Generate Unique Condition 02 Column
meta_data %>%
  group_by(growthPhase,carbonSource,Mg_mM_Levels,Na_mM_Levels) %>%
  summarize(replicaNum=length(growthPhase)) ->sampleSizeDf02

sampleSizeDf02%>%
  group_by()%>%
  mutate(uniqueCondition02=sprintf("unique_condition_%02d",as.numeric(rownames(.)))) -> sampleSizeDf02



sampleSizeDf02%>%
  select(-replicaNum)%>%
  left_join(meta_data,.)->meta_data


# generate harvestDate, notesNGS, cellTotal, cellsPerTube, tubesInStorage columns
names(meta_data)[which(names(meta_data)=="HarvestDate")]<-"harvestDate"
names(meta_data)[which(names(meta_data)=="Notes.NGS.")]<-"notesNGS"
names(meta_data)[which(names(meta_data)=="Cell_total")]<-"cellTotal"
names(meta_data)[which(names(meta_data)=="cells_per_tube")]<-"cellsPerTube"
names(meta_data)[which(names(meta_data)=="X.tubesinstorage")]<-"tubesInStorage"

meta_data %>% select(-tubesInStorage, -notesNGS)->meta_data




meta_data %>% 
  select(-Batch,-GrowthConditions, -concentration) %>%
  arrange(dataSet)->meta_data

meta_data %>% 
  left_join(.,doubling_times)->meta_data
#********************************


#********************************
# Generate meta_rna and meta_protein
meta_data %>%
  dplyr::filter(RNA_Data_Freq!=0) %>%
  dplyr::select(-RNA_Data_Freq,-Protein_Data_Freq)->meta_rna

meta_data %>%
  dplyr::filter(Protein_Data_Freq==2) %>%
  dplyr::mutate(dataSet=paste0(dataSet,"_1")) %>%
  dplyr::mutate(sampleNum=paste0(sampleNum,"_1"))-> meta_protein_addition

meta_data$sampleNum=as.character(meta_data$sampleNum)
meta_data$sampleNum=paste0(meta_data$sampleNum,"_0")
meta_data$dataSet=paste0(meta_data$dataSet,"_0")

dplyr::rbind_list(meta_data,meta_protein_addition)->meta_protein


meta_protein %>%
  dplyr::arrange(dataSet)%>%
  dplyr::filter(Protein_Data_Freq!=0) %>%
  dplyr::select(-RNA_Data_Freq,-Protein_Data_Freq)->meta_protein
#********************************


# *******************************
savedFilename=paste0("../a_results/","metaRawData.Rda")
savingList=c("meta_data",
             "meta_rna",
             "meta_protein",
             "sampleSizeDf")
save(list=savingList,file=savedFilename)



# Save as csv files
savedFilename=paste0("../a_results/","metaData.csv")
write.csv(meta_data, file = savedFilename, row.names = FALSE)

savedFilename=paste0("../a_results/","metaRNA.csv")
write.csv(meta_rna, file = savedFilename, row.names = FALSE)

savedFilename=paste0("../a_results/","metaProtein.csv")
write.csv(meta_protein, file = savedFilename, row.names = FALSE)

savedFilename=paste0("../a_results/","sampleSizeDf.csv")
write.csv(sampleSizeDf, file = savedFilename, row.names = FALSE)

