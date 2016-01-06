# trend-GO analysis of MgSO4 

###*****************************
# INITIAL COMMANDS TO RESET THE SYSTEM

rm(list = ls())
if (is.integer(dev.list())){dev.off()}
cat("\014")

# Set Working Directory
setwd('/Users/umut/GitHub/AG3C_Analyze/initialPaper03c_change_wrt_parameters_RNA&Protein/') # mac computer
###*****************************




###*****************************
# INSTALL LIBRARIES

library("Biobase") 
library("DESeq2")
library("dplyr")
library("tidyr")
library("lazyeval")

# For Machine Learning & correlations
library("e1071")
library("caret")
library("ICC")
library("gtools")

# For Plotting
library("ggplot2")
library("RColorBrewer")
library("grid")
library("gridExtra")
library("cowplot")
###*****************************


###*****************************
# PARAMETERS FOR DATA
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

test_for="lactate" # "growthTime_hr", 
                     # "exponential","stationary","late_stationary",
                    # "Mg_mM", "Na_mM"
                    # "glucose", "glycerol", "gluconate", "lactate"
                    # "carbonSource"
test_type="spearman" # "spearman", "ICC"

dictionaryPick="BL"
# there are 4 alternatves for dictionary
# barricks own dictionary -Barrick Lab (BL)
# what I found on web (WEB)
# a combination of 2 -Combined (COM)
# Protein dictionary (PD)

filterThreshold_p=1 # give some data direct 0 If the value is 1 we include all data

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
extenedFileName=paste0("extened_",normalizationMethodChoice,"_",deSeqNormChoice,"_",
              filterTypeChoice,"_",longNameList,"_",test_type,
              "_test_",test_for,"_dic",dictionaryPick) 
goInputFileName=paste0("goInput_",normalizationMethodChoice,"_",deSeqNormChoice,"_",
                       filterTypeChoice,"_",longNameList,"_",test_type,"_test_",
                       test_for,"_dic",dictionaryPick) 
goDetailedInputFileName=paste0("goDetailedInput_",normalizationMethodChoice,"_",deSeqNormChoice,"_",
                       filterTypeChoice,"_",longNameList,"_",test_type,"_test_",test_for,
                       "_dic",dictionaryPick) 
conditionName=paste0("condition_",normalizationMethodChoice,"_",deSeqNormChoice,"_",
                     filterTypeChoice,"_",longNameList) 
###*****************************


###*****************************
#LOAD FILES
loadedFilename=paste0("../initialPaper01r/",step03,".RData")
load(file=loadedFilename)

if (dataTypeChoice=="mrna"|dataTypeChoice=="rna"){
  if(dictionaryPick=="BL"){
    dictionary=read.csv(paste0("../generateDictionary/","nameDictionary_RNA_barrick.csv")) 
  }
  if(dictionaryPick=="WEB"){
    dictionary=read.csv(paste0("../generateDictionary/","nameDictionary_RNA.internet.csv")) 
  }
  if(dictionaryPick=="COM"){
    # gene name dictionary comb
    dictionary=read.csv(paste0("../generateDictionary/","nameDictionaryCombined_RNA.csv")) 
    
  }
  names(dictionary)[1]<-"ID"
}

if (dataTypeChoice=="protein"|dataTypeChoice=="protein_wo_NA"|dataTypeChoice=="protein_wo_NAx6"){
  dictionary=read.csv(paste0("../generateDictionary/","nameDictionary_Protein.csv")) 
  names(dictionary)[1]<-"ID"
}


#*********************************


#*********************************
# Install Go Functions
source("GE01_GO_Functions.R")
#*********************************


#*********************************
#Rename Loaded Files
assign(x = "mainData", value = get(step03))
assign(x = "condition", value = get(conditionName))
#*********************************


###*****************************
#Get rid of exact replicasin proteins and check if it hurts rna data 
# Calculate the average of repeated conditions
mainDataFrame<-as.data.frame(mainData)
mainDataFrame %>% mutate(ID=row.names(mainDataFrame))->mainDataFrameTidy
endVal=ncol(mainDataFrame)
mainDataFrameTidy %>%
  tidyr::gather(condition,read,1:endVal)%>%
  dplyr::mutate(conditionShort=gsub("_[[:digit:]]$","",condition))->mainDataFrameTidy
mainDataFrameTidy %>% 
  group_by(ID,conditionShort)%>%
  dplyr::summarise(read=mean(read))->mainDataFrameSummary

mainDataFrameSummary %>%
  dplyr::group_by()%>%  
  tidyr::spread(key = conditionShort, value = read)->mainDataFrameSummary

row.names(mainDataFrameSummary)<-mainDataFrameSummary$ID
mainDataFrameSummary %>%
  dplyr::select(-ID)->mainDataFrameSummary
mainData=as.matrix(mainDataFrameSummary)
###*****************************


###*****************************
#Get rid of exact replicasin proteins and check if it hurts rna data 
#Reduce number of rows in condition and get rid of repeats

if(unique(nchar(as.vector(condition$dataSet))==10)){
  condition=condition[which(substr(as.vector(condition$dataSet),10,10) %in% c(0)),]
  condition$dataSet<-gsub("_[[:digit:]]$","",as.vector(condition$dataSet))
  condition$sampleNum<-gsub("_[[:digit:]]$","",as.vector(condition$sampleNum))}
###*****************************


#*********************************
# Combined DF. includes condition, mainData, dictionary 
mainDataDF=as.data.frame(t(mainData))
geneNames=colnames(mainDataDF)
dataSetList=row.names(mainDataDF)
mainDataDF%>%
  dplyr::mutate(dataSet=dataSetList)->mainDataDF
extendedDataFrame=dplyr::left_join(condition,mainDataDF)
extendedDataFrame %>%
  dplyr::mutate(glucose=(carbonSource=="glucose")+0) %>%
  dplyr::mutate(glycerol=(carbonSource=="glycerol")+0) %>%
  dplyr::mutate(lactate=(carbonSource=="lactate")+0) %>%
  dplyr::mutate(gluconate=(carbonSource=="gluconate")+0) %>%
  dplyr::mutate(exponential=(growthPhase=="exponential")+0) %>%
  dplyr::mutate(stationary=(growthPhase=="stationary")+0) %>%
  dplyr::mutate(late_stationary=(growthPhase=="late_stationary")+0)->extendedDataFrame


geneLocationList=match(geneNames,colnames(extendedDataFrame))
extendedDataFrame %>%
  tidyr::gather(ID, amount, geneLocationList) -> extendedDataFrame_tidy
dplyr::left_join(extendedDataFrame_tidy,dictionary) -> extendedDataFrame_tidy

remove(mainDataDF,extendedDataFrame, dataSetList)
#*********************************


#*********************************
# Calculate Spearman correlation p values
if(test_type=="spearman"){
  summaryData_detailed=spearman_function_v2(main_data=extendedDataFrame_tidy,
                                            x_axis_input=test_for,
                                            y_axis_input="amount")
}

if(test_type=="ICC"){
  summaryData_detailed=ICC_function(main_data=extendedDataFrame_tidy,
                                    x_axis_input=test_for,
                                    y_axis_input="amount")
}
#*********************************


#*********************************
# Join the results to extened, rename and save big file
extendedDataFrame_tidy %>%
  dplyr::left_join(.,summaryData_detailed)->extendedDataFrame_tidy

assign(x = extenedFileName, value = extendedDataFrame_tidy)


saveFileName=paste0("../initialPaper03r/",extenedFileName,".RData")
save(list = c(step02,extenedFileName),file = saveFileName)
#*********************************


#*********************************
# Clean Out Summary and add new rescaled columns

summaryData_detailed=summaryData_detailed[which(!grepl("ECB_",summaryData_detailed$gene_name)),]
summaryData_detailed=summaryData_detailed[which(!grepl("^$",summaryData_detailed$gene_name)),]

#get rid of repeated columns
summaryData_detailed %>% 
  dplyr::group_by(gene_name) %>%
  dplyr::arrange(ID)%>%
  dplyr::mutate(index=seq(1,length(gene_name))) %>%
  dplyr::group_by() %>%
  dplyr::mutate(gene_name=paste0(gene_name,"_",index))%>%
  dplyr::select(-index)->summaryData_detailed


if(test_type=="spearman"){
  summaryData_detailed %>% 
    dplyr::mutate(p_related=p_related*p_related_sign) %>%
    dplyr::select(gene_name,p_related) %>%
    dplyr::mutate(p_related=ifelse(p_related> -log10(filterThreshold_p) |
                                              p_related< log10(filterThreshold_p), 
                                            p_related, 0)) %>% # new try
    dplyr::arrange(desc(p_related)) %>%
    dplyr::group_by(p_related) %>% # new try
    dplyr::mutate(gene_name=sample(gene_name))->summaryData # new try
  
  colnames(summaryData)<-c("gene","p_related_wrep")
}

if(test_type=="ICC"){
}


#*********************************


#*********************************
# # addFake to summary data
# summaryData %>%
#   dplyr::group_by() %>%
#   dplyr::mutate(gene=sample(gene))->summaryData
#*********************************


#*********************************


saveFileName=paste0("../initialPaper03r/",goInputFileName,".csv")
write.csv(summaryData,file=saveFileName,row.names=F,quote=F)

saveFileName=paste0("../initialPaper03r/",goDetailedInputFileName,".csv")
write.csv(summaryData_detailed,file=saveFileName,row.names=F,quote=F)
#*********************************


#*********************************
# OPTIONAL FIGURE
# most decrease #asr

chosenGeneName="idnO"
extendedDataFrame_tidy %>%
  dplyr::filter(gene_name==chosenGeneName) ->temp

fig01<-ggplot(temp, aes_string(x=test_for, y="amount")) +
  geom_point(size=3)+
 # geom_vline(xintercept = 0.8, colour="red", linetype = "dashed")+
  ggtitle(paste0(chosenGeneName,"_",unique(temp$p.adj_related)))+
  xlab(test_for)+
  ylab("normalized mRNA reads")+
  #scale_x_continuous(breaks=c(0,1))+
  scale_x_discrete(labels=unique(temp[test_for]))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
print(fig01)

