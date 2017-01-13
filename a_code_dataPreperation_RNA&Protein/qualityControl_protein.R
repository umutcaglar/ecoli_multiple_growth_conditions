# Data Quality control

# This file will generate quality control related figures / results for given data sets
# The idea is to look at triplets for  

#***************************
# Initial Command to Reset the System
rm(list = ls())
if (is.integer(dev.list())){dev.off()}
cat("\014")
#***************************


###*****************************
# Set Working Directory
# One needs to arrange the correct pathway if this is not umut's computer ;)
if(as.vector(Sys.info()["effective_user"]=="umut"))
{setwd(paste0("/Users/umut/GitHub/ecoli_multiple_growth_conditions/",
              "a_code_dataPreperation_RNA&Protein//"))} # mac computer
###*****************************


#***************************************
#Load Libraries and functions

# data organization
require("dplyr")
require("tidyr")

# Figures
require("ggplot2")
require("cowplot")
#***************************************


###*****************************
# Load Functions
source(file = "data_naming_functions.R")
###*****************************


#***************************************
# We will use
# * protein data
# * we use all samples
# * we sum technical replicates we round the raw data
# * we do a +1 normalization for size factors
# * No filtering

# The main data naming function that controls sub functions.
dataName=name_data(initialValue="resDf", # can be c("genes0.05","genes_P0.05Fold2","resDf")
                   dataType = "protein", # can be "rna", "mrna", "protein", "protein_wo_NA" # Use protein instead of "protein_wo_NA"
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
                   filterGenes = "noMatchFilter", # can be either "noFilter", or any combination of c("meanFilter", "maxFilter", "sdFilter", "noMatchFilter")
                   threshold=NA, # the threshold value for "meanFilter", "maxFilter", "sdFilter" can be  c(meanFilter=5,maxFilter=3,sdFilter=7)
                   roundData=TRUE,
                   sumTechnicalReplicates=TRUE,
                   deSeqSfChoice="p1Sf", # can be "regSf", "p1Sf", "noSf"
                   normalizationMethodChoice= "noNorm", # can be "vst", "rlog", "log10", "noNorm"
                   test_for = "noTest")  # works only if normalizationMethodChoice == noNorm
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource", "noTest")

dataNameDF=as.data.frame(dataName[1])

metaDataName=dataNameDF
metaDataName$objectName.initial="metaData"

dataName=paste(dataNameDF,collapse = "_")
metaDataName=paste(metaDataName,collapse = "_")

mainDataFrame=read.csv(file = paste0("../a_results/",dataName,".csv"),header = TRUE,row.names = 1)
condition=read.csv(file = paste0("../a_results/",metaDataName,".csv"),header = TRUE)
###*****************************

###*****************************
# Prepeare Data for grouping

# Divide the condition into two parts as time series and others
# in time series data group with respect to sample hour
# in other data group with respect to growth phase
condition %>% 
  dplyr::filter(experiment %in% c("glucose_time_course","glycerol_time_course")) -> condition_timeC
condition %>%
  dplyr::filter(!experiment %in% c("glucose_time_course","glycerol_time_course")) -> condition_other

condition_timeC %>%
  dplyr::mutate(uniqueCondition= paste0(carbonSource, "_",
                                        Mg_mM, "mMMg_",
                                        Na_mM, "mMNa_",
                                        sprintf("%03d", growthTime_hr),"hr"))->condition_timeC 

condition_other %>%
  dplyr::mutate(uniqueCondition= paste0(carbonSource,"_", 
                                        Mg_mM,"mMMg_", 
                                        Na_mM, "mMNa_",
                                        growthPhase))->condition_other

rbind_list(condition_timeC, condition_other)->condition02

condition02 %>%
  dplyr::group_by(uniqueCondition) %>%
  dplyr::summarise(numSamples=n())->condition_summary
###*****************************


###*****************************
# For each group do log2((sample_A+1)/(sample_B+1)) for each pair
counter04=0
for (counter01 in 1: nrow(condition_summary))
{
  condition02 %>%
    dplyr::filter(condition02$uniqueCondition==condition_summary$uniqueCondition[counter01])->subCondition
  subConditionVector=as.vector(subCondition$dataSet)
  
  if (condition_summary$numSamples[counter01]!=1)
  {
    for (counter02 in 1: (length(subConditionVector)-1))
    {
      for (counter03 in (counter02+1):length(subConditionVector))
      {
        counter04=counter04+1
        
        sample_A=mainDataFrame[[subConditionVector[counter02]]]
        sample_B=mainDataFrame[[subConditionVector[counter03]]]
        
        sample_A[is.na(sample_A)]<-0
        sample_B[is.na(sample_B)]<-0
        
        vs_name=gsub("MURI_","",paste(subConditionVector[counter02],subConditionVector[counter03],sep="_vs_"))
        uniqueCondition=unique(subCondition$uniqueCondition)
        if(uniqueCondition=="glucose_0.8mMMg_100mMNa_stationary"){browser()}
        
        ratio=log2((sample_A+1)/(sample_B+1))
        if(sum(is.na(ratio))>0){browser()}
        ratioDf=as.data.frame(x=ratio,row.names = as.vector(mainDataFrame$gene_ID))
        ratioDf%>%
          dplyr::mutate(vs_name=vs_name,
                        uniqueCondition=uniqueCondition)->ratioDf
        
        if(counter04==1){ratioDfL=ratioDf}
        if(counter04!=1){ratioDfL=rbind_list(ratioDfL,ratioDf)}
        
      }
    }
  }
}
###*****************************


###*****************************
# Generate summary tables
ratioDfL %>%
  dplyr::mutate(sample_A=gsub("_vs_*.*","",vs_name),
                sample_B=gsub("*.*_vs_","",vs_name))->ratioDfL

ratioDfL %>%
  dplyr::group_by(vs_name) %>%
  dplyr::summarise(meanVal=mean(ratio), std=sd(ratio), medianVal=median(ratio),
                   uniqueCondition=unique(uniqueCondition)) %>%
  dplyr::mutate(sample_A=gsub("_vs_*.*","",vs_name),
                sample_B=gsub("*.*_vs_","",vs_name))->ratioDfL_summary

ratioDfL_summary %>%
  dplyr::group_by(uniqueCondition) %>%
  dplyr::summarise(medianMax=max(abs(medianVal)),# associated wih size factors
                   stdMax=max(std), # associated with variety of condition
                   sd_sd=sd(std))->ratioDfL_sumsummary # associated with sample problems

# the function finds the farest sample from all the other samples un a given unique condition
# then calculate the average distance of this point to other points in the space (say A)
# after that it calculates the average distance between remainin points (say B)
# Finally it calculates ratioValue=A/B
# It returns both ratioValue and the name of the farest point.

orderDifferencesFunction<-function(x)
{
  x%>%
    dplyr::mutate(sample_C=sample_A,
                  sample_A=sample_B,
                  sample_B=sample_C)%>%
    dplyr::select(-sample_C)->x2
  
  dplyr:::rbind_list(x, x2)->x_comb
  remove(x2)
  
  x_comb %>%
    dplyr::group_by(sample_A) %>%
    dplyr::summarise(stdMean=mean(std)) %>%
    dplyr::arrange(desc(stdMean))->x_summary
  
  x_comb %>% 
    dplyr::filter(sample_A==x_summary$sample_A[1])->x_upperDf
  
  x %>% 
    dplyr::filter(sample_A!=x_summary$sample_A[1] & sample_B!=x_summary$sample_A[1])->x_lowerDf
  
  x_upper=x_upperDf$std
  x_lower=x_lowerDf$std
  
  x_upper_mean=mean(x_upper)
  x_lower_mean=mean(x_lower)
  
  ratioValue=x_upper_mean/x_lower_mean
  oddOne=paste0("MURI_",x_summary$sample_A[1])
  return(list(ratioValue=ratioValue, oddOne=oddOne))
}

ratioDfL_summary %>% 
  dplyr::group_by(uniqueCondition) %>%
  dplyr::do(data.frame(returnObj=orderDifferencesFunction(.)))->ratioDfL_ratioValue

colnames(ratioDfL_ratioValue)<-gsub("returnObj.","",colnames(ratioDfL_ratioValue))


ratioDfL_sumsummary <- left_join(ratioDfL_sumsummary,condition_summary)
ratioDfL_sumsummary <- left_join(ratioDfL_sumsummary,ratioDfL_ratioValue)
remove(ratioDfL_ratioValue)

ratioDfL_sumsummary %>%
  dplyr::arrange(desc(ratioValue))->ratioDfL_sumsummary

ratioDfL$uniqueCondition <- factor(ratioDfL$uniqueCondition,
                                   levels = as.vector(ratioDfL_sumsummary$uniqueCondition))

ratioDfL_summary$uniqueCondition <- factor(ratioDfL_summary$uniqueCondition,
                                           levels = as.vector(ratioDfL_sumsummary$uniqueCondition))
###*****************************


###*****************************
# Draw histograms associated with quality control data

# arrange colors
ratioDfL_summary %>%
  dplyr::group_by(uniqueCondition)%>%
  dplyr::mutate(colorVariation=paste0("color",sprintf("%02.f",seq(1:n()))))->ratioDfL_summaryC

ratioDfL<-left_join(ratioDfL,ratioDfL_summaryC)

manuelColors=c("#9d007f",
               "#64c007",
               "#3695ff",
               "#ff5839",
               "#28bdbc",
               "#ff2562",
               "#01734b",
               "#e09e3f",
               "#642b0e",
               "#5f7400")
#remove(ratioDfL_summaryC)

# draw figure
fig01<-ggplot2::ggplot(ratioDfL,aes(x=ratio, group=vs_name, color=colorVariation))+
  geom_line(aes(y=..density..), stat="density",adjust=5,size=2) +
  scale_colour_manual(values = manuelColors)+
  facet_wrap(~ uniqueCondition, ncol = 5,scales = "free")+
  theme_bw()+ #can be theme_bw or theme_classic
  scale_x_continuous(limits = c(-10,10))+
  geom_vline(xintercept = 0, colour="red", linetype = "longdash")+
  theme( axis.text.x=element_text(size=16),
         axis.text.y=element_text(size=16),
         axis.title.x=element_text(size=16),
         axis.title.y=element_text(size=16),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14) )

cowplot::save_plot(filename = "../a_figures/qualityControl_protein.pdf",plot = fig01,ncol = 5,nrow = 7)

ratioDfL_summary%>%
  dplyr::mutate(sample_C=sample_A,
                sample_A=sample_B,
                sample_B=sample_C)%>%
  dplyr::select(-sample_C)->ratioDfL_summary_B

dplyr:::rbind_list(ratioDfL_summary, ratioDfL_summary_B)->ratioDfL_summary_J

# Draw the figure that shows the distribution of standart devaiation between divisions of samples
fig02<-ggplot2::ggplot(ratioDfL_summary_J,aes(x=sample_A, y=sample_B))+
  geom_tile(aes(fill=std), color="grey")+
  geom_text(aes(label=sprintf("%1.2f", std)))+
  facet_wrap(~ uniqueCondition, ncol = 5, scales = "free")+
  scale_fill_gradient(low = "White", high = "Black",limits=c(0,5),name = "Std")+
  theme_classic()+ #can be theme_bw or theme_classic
  theme( axis.text.x=element_text(size=16),
         axis.text.y=element_text(size=16),
         axis.title.x=element_text(size=16),
         axis.title.y=element_text(size=16),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14) )

print(fig02)

cowplot::save_plot(filename = "../a_figures/qualityControl_protein_heatMap.pdf",plot = fig02,ncol = 5,nrow = 7)
###*****************************


###*****************************
# save the data frame that shows the odd samples and how odd they are
write.csv(x = ratioDfL_sumsummary,file = "../a_results/odd_protein_samples.csv")
###*****************************