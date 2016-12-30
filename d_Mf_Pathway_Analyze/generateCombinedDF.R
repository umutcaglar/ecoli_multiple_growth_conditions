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
{setwd(paste0("/Users/umut/GitHub/ecoli_multiple_growth_conditions/d_results/"))} # mac computer
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
# get file names convert into df divide into parts
dir() %>%
  data.frame(fullFileName=.) %>%
  dplyr::mutate(ffN=fullFileName) %>%
  dplyr::group_by(fullFileName) %>%
  dplyr::mutate(ffN=gsub(pattern = "_mM_Levels","",ffN)) %>%
  dplyr::mutate(ffN=gsub(pattern = ".csv","",ffN)) %>%
  dplyr::filter(!grepl(pattern = "_o",x = fullFileName)) %>%
  tidyr::separate(data = ., col = ffN, into = c("v1","thresholds_for_input_genes","dataType", 
                                                "sumTechnicalReplicates", "preEleminationOfData", "experiments",
                                                "carbonSource","Mg","Na",
                                                "growthPhase","filter","plus1SizeFactor",
                                                "normalization","empty","investigatedEffect",
                                                "empty2","testVSbase","database","v19"),sep = "_") %>%
  dplyr::group_by() %>%
  dplyr::filter(v1=="ez") %>%
  dplyr::select(-v1, -sumTechnicalReplicates, -preEleminationOfData, -experiments,
                -filter, -plus1SizeFactor, -normalization, -empty, -empty2, -v19, -database)->allResultFiles

for(counter01 in 1:nrow(allResultFiles))
{
  print(counter01)
  m<-read.csv(file = as.vector(allResultFiles$fullFileName[counter01]),row.names = 1)
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
fullList %>%dplyr::select(-X., gene_number, -abs_score, -numSigP, -numSigN)->fullList
###*****************************


###*****************************
write.csv(x = fullList, file = "combinedResultList_DAVID.csv")
###*****************************
###*****************************