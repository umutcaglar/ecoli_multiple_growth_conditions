# Analayze Kegg Pathways DESeq + DAVID

# Aim of the code is to find the genes in kegg pathways and send files to generate figures


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
              "d_Mf_Pathway_Analyze/"))} # mac computer
###*****************************


###*****************************
# DOWNLOAD LIBRARIES
require("Biobase") 
require("DESeq2")
require("dplyr")
require("tidyr")
###*****************************


###*****************************
# Load Files
source("../c_code_change_wrt_variables_RNA&Protein/data_naming_functions.R")
###*****************************


###*****************************
# load reference libraries
kegg_pathway_david_mrna_ref_2009_tidy<-read.table(file = paste0("../d_Mf_Pathway_Analyze/ReferenceFiles/",
                                                                "kegg_pathway_david_mrna_ref_2009_tidy.txt"),
                                                  sep = "\t",header = TRUE,fill = TRUE, quote = "")


kegg_pathway_david_protein_ref_2009_tidy<-read.table(paste0("../d_Mf_Pathway_Analyze/ReferenceFiles/",
                                                            "kegg_pathway_david_protein_ref_2009_tidy.txt"),
                                                     sep = "\t",header = TRUE,fill = TRUE, quote = "")
###*****************************


###*****************************
# Download the DAVID input and output
dataName=name_data(initialValue="genes0.05", # can be "genes0.05"
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
                   carbonSourceVector = "S", # can be any sub combination of "SYAN"
                   MgLevelVector = c("baseMg"), # can be "lowMg","baseMg","highMg" // "allMg"
                   NaLevelVector = c("allNa"), # can be "baseNa","highNa" // "allNa"
                   growthPhaseVector = c("exponential"), # can be "exponential","stationary","late_stationary" // "allPhase"
                   filterGenes = "noFilter", # can be "noFilter", "meanFilter", "maxFilter", "sdFilter" 
                   threshold=NA, # the threshold value for "meanFilter", "maxFilter", "sdFilter"
                   roundData=TRUE,
                   sumTechnicalReplicates=TRUE,
                   deSeqSfChoice="p1Sf", # can be "regSf", "p1Sf"
                   normalizationMethodChoice= "noNorm", # can be "vst", "rlog", "log10", "noNorm"
                   test_for = "Na_mM_Levels")  # works only if normalizationMethodChoice == noNorm
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")

objectName=paste(dataName$objectName,collapse = "_")

kegg_input<-read.csv(file = paste0("../c_results/DeSeq2_diffGene_Results/",objectName,".csv"),
                     header = TRUE)
kegg_result<-read.table(file=paste0("../c_results/david_results/",objectName,"_kegg.txt"),
                        sep = "\t",header = TRUE)
kegg_result<-kegg_result[grep("eco",as.vector(kegg_result$Term)),]
###*****************************


###*****************************
# Find KEGG output pathways in ref 
# do the consistency check if everything is in

for(counter01 in 1:nrow(kegg_result))
{
  thePathway=as.vector(kegg_result[["Term"]][counter01])
  kegg_pathway_david_mrna_ref_2009_tidy %>%
    dplyr::filter(KEGG_PATH_Name==thePathway)->temp
  
  m2<-as.vector(temp$ID)
  m1<-strsplit(x=as.vector(kegg_result$Genes[counter01]),split = ",")
  m1<-sub(" ","",noquote(m1[[1]]))
  m1<-paste0(tolower(substr(start = 1,stop = 3,x = m1)),substr(start = 4,stop = 4,x = m1))
  print(all(m1 %in% m2))
}





