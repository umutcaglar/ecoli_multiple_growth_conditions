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


### DO FOR mRNA DATA ###
###*****************************
# Load Data
condition<-read.csv(file = "../a_results/metaRNA.csv")
mainDataFrame<-read.csv(file = "../a_results/rnaMatrix_mRNA.csv")
###*****************************


###*****************************
# Prepeare Data for grouping

# A) mainDataFrame
mainDataFrame %>% dplyr::select(-gene_Type)->mainDataFrame# remove the firts column named as "gene_Type"

# B) condition
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
                   uniqueCondition=unique(uniqueCondition))%>%
  dplyr::mutate(sample_A=gsub("_vs_*.*","",vs_name),
                sample_B=gsub("*.*_vs_","",vs_name))->ratioDfL_summary

ratioDfL_summary %>%
  dplyr::group_by(uniqueCondition) %>%
  dplyr::summarise(medianMax=max(abs(medianVal)),# associated wih size factors
                   stdMax=max(std), # associated with variety of condition
                   sd_sd=sd(std))->ratioDfL_sumsummary # associated with sample problems

ratioDfL_sumsummary <- left_join(ratioDfL_sumsummary,condition_summary)
###*****************************


###*****************************
# Draw histograms associated with quality control data
fig01<-ggplot2::ggplot(ratioDfL,aes(x=ratio, group=vs_name))+
  geom_line(aes(y=..density..), stat="density") +
  facet_wrap(~ uniqueCondition, ncol = 4)+
  theme_bw()+ #can be theme_bw or theme_classic
  xlim(-10,10)+
  geom_vline(xintercept = 0, colour="red", linetype = "longdash")+
  theme( axis.text.x=element_text(size=16),
         axis.text.y=element_text(size=16),
         axis.title.x=element_text(size=16),
         axis.title.y=element_text(size=16),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14) )

print(fig01)
cowplot::save_plot(filename = "../a_figures/qualityControl_mRNA.pdf",plot = fig01,ncol = 4,nrow = 12)


ratioDfL_summary%>%
  dplyr::mutate(sample_C=sample_A,
                sample_A=sample_B,
                sample_B=sample_C)%>%
  dplyr::select(-sample_C)->ratioDfL_summary_B

dplyr:::rbind_list(ratioDfL_summary, ratioDfL_summary_B)->ratioDfL_summary_J

fig02<-ggplot2::ggplot(ratioDfL_summary_J,aes(x=sample_A, y=sample_B))+
  geom_tile(aes(fill=std), color="grey")+
  geom_text(aes(label=sprintf("%1.2f", std)))+
  facet_wrap(~ uniqueCondition, ncol = 4, scales = "free")+
  scale_fill_gradient(low = "White", high = "Black",limits=c(0,5),name = "Std")+
  theme_classic()+ #can be theme_bw or theme_classic
  theme( axis.text.x=element_text(size=16),
         axis.text.y=element_text(size=16),
         axis.title.x=element_text(size=16),
         axis.title.y=element_text(size=16),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14) )

print(fig02)

cowplot::save_plot(filename = "../a_figures/qualityControl_mRNA_heatMap.pdf",plot = fig02,ncol = 4,nrow = 12)
###*****************************