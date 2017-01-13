# results DF summary

###*****************************
# INITIAL COMMANDS TO RESET THE SYSTEM
rm(list = ls())
if (is.integer(dev.list())){dev.off()}
cat("\014")
set.seed(14159)
###*****************************


###*****************************
# Set Working Directory
# One needs to arrange the correct pathway if this is not umut's computer ;)
if(as.vector(Sys.info()["effective_user"]=="umut"))
{setwd(paste0("/Users/umut/GitHub/ecoli_multiple_growth_conditions/c_results/"))} # mac computer
###*****************************


###*****************************
# REQUIRED LIBRARIES
# Data tracking
require("DESeq2")
require("dplyr")
require("tidyr")

# Graphing
require("ggplot2")
require("cowplot")

# Text manipulation
require("stringr")
###*****************************


###*****************************
# Combine for Genes
# get file names convert into df divide into parts
dir() %>%
  data.frame(fullFileName=.) %>%
  dplyr::mutate(ffN=fullFileName) %>%
  dplyr::group_by(fullFileName) %>%
  dplyr::filter(grepl(pattern = "genes_",x = fullFileName)) %>%
  dplyr::mutate(ffN=gsub(pattern = "_mM_Levels","",ffN)) %>%
  dplyr::mutate(ffN=gsub(pattern = ".csv","",ffN)) %>%
  tidyr::separate(data = ., col = ffN, into = c("v1","thresholds_for_input_genes","dataType", 
                                                "sumTechnicalReplicates", "preEleminationOfData", "experiments",
                                                "carbonSource","Mg","Na",
                                                "growthPhase","filter","plus1SizeFactor",
                                                "normalization","empty","investigatedEffect",
                                                "empty2","testVSbase","database","v19"),sep = "_") %>%
  dplyr::group_by() %>%
  dplyr::select(-v1, -sumTechnicalReplicates, -preEleminationOfData, -experiments,
                -filter, -plus1SizeFactor, -normalization, -empty, -empty2, -v19, -database)->allResultFiles


for(counter01 in 1:nrow(allResultFiles))
{
  print(counter01)
  m<-read.csv(file = as.vector(allResultFiles$fullFileName[counter01]))
  #colnames(m)<-gsub(pattern = "KEGG_", replacement = "KEGG_MF_", colnames(m))
  #colnames(m)<-gsub(pattern = "MF_", replacement = "KEGG_MF_", colnames(m))
  #colnames(m)<-gsub(pattern = "_MF", replacement = "", colnames(m))
  m %>% dplyr::mutate(fullFileName = as.vector(allResultFiles$fullFileName[counter01])) ->m
  
  if(nrow(m)!=0)
  {
    if(counter01==1){fullList<-m}
    if(counter01!=1){fullList=bind_rows(fullList,m)}
  }
}
###*****************************


###*****************************
# Combine with experiment properties
dplyr::left_join(fullList, allResultFiles)->fullList
###*****************************


###*****************************
write.csv(x = fullList, file = "combinedDifferentiallyExpressedGenes_DeSeq.csv")
###*****************************


###*****************************
# Combine for Data Frames
# get file names convert into df divide into parts
dir() %>%
  data.frame(fullFileName=.) %>%
  dplyr::mutate(ffN=fullFileName) %>%
  dplyr::group_by(fullFileName) %>%
  dplyr::filter(grepl(pattern = "resDf_",x = fullFileName)) %>%
  dplyr::mutate(ffN=gsub(pattern = "_mM_Levels","",ffN)) %>%
  dplyr::mutate(ffN=gsub(pattern = ".csv","",ffN)) %>%
  tidyr::separate(data = ., col = ffN, into = c("v1","thresholds_for_input_genes","dataType", 
                                                "sumTechnicalReplicates", "preEleminationOfData", "experiments",
                                                "carbonSource","Mg","Na",
                                                "growthPhase","filter","plus1SizeFactor",
                                                "normalization","empty","investigatedEffect",
                                                "empty2","testVSbase","database","v19"),sep = "_") %>%
  dplyr::group_by() %>%
  dplyr::select(-v1, -sumTechnicalReplicates, -preEleminationOfData, -experiments,
                -filter, -plus1SizeFactor, -normalization, -empty, -empty2, -v19, -database)->allResultFiles


for(counter01 in 1:nrow(allResultFiles))
{
  print(counter01)
  m<-read.csv(file = as.vector(allResultFiles$fullFileName[counter01]))
  #colnames(m)<-gsub(pattern = "KEGG_", replacement = "KEGG_MF_", colnames(m))
  #colnames(m)<-gsub(pattern = "MF_", replacement = "KEGG_MF_", colnames(m))
  #colnames(m)<-gsub(pattern = "_MF", replacement = "", colnames(m))
  m %>% dplyr::mutate(fullFileName = as.vector(allResultFiles$fullFileName[counter01])) ->m
  
  if(nrow(m)!=0)
  {
    if(counter01==1){fullList<-m}
    if(counter01!=1){fullList=bind_rows(fullList,m)}
  }
}
###*****************************


###*****************************
# Combine with experiment properties
dplyr::left_join(fullList, allResultFiles)->fullList
###*****************************

###*****************************
write.csv(x = fullList, file = "combinedOutputDF_DeSeq.csv")
###*****************************