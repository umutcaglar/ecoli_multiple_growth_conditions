## DeSeq Normalization Cleaned 03

# The aim of the code is to generate normalized data matrix.
# The work flow compses of four parts
# Pick up the samples
# pick up the rows
# calculate size factors
# do the normalization


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
              "c_code_change_wrt_variables_RNA&Protein/"))} # mac computer
###*****************************


###*****************************
# DOWNLOAD LIBRARIES
require("Biobase") 
require("DESeq2")
require("dplyr")
require("tidyr")
###*****************************


###*****************************
#Load Functions
source("data_normalization_functions.R")
###*****************************


###*****************************
saveFiles=FALSE
# The data filtering function that controls sub functions.
mainData=filter_data(dataType = "protein_wo_NA", # can be "rna", "mrna", "protein", "protein_wo_NA"
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
                     referenceLevels=c("stationary",
                                       "baseMg", 
                                       "baseNa", 
                                       "glucose", 
                                       "glucose_time_course"),
                     experimentVector = c("allEx"), # can be "Stc","Ytc","Nas","Agr","Ngr","Mgl","Mgh" // "allEx"
                     carbonSourceVector = "SY", # can be any sub combination of "SYAN"
                     MgLevelVector = c("baseMg"), # can be "lowMg","baseMg","highMg" // "allMg"
                     NaLevelVector = c("baseNa"), # can be "baseNa","highNa" // "allNa"
                     # can be "exponential","stationary","late_stationary" // "allPhase"
                     growthPhaseVector = c("stationary"), 
                     filterGenes = "noFilter", # can be "noFilter", "meanFilter", "maxFilter", "sdFilter" 
                     threshold=NA, # the threshold value for "meanFilter", "maxFilter", "sdFilter"
                     roundData=TRUE,
                     sumTechnicalReplicates=TRUE,
                     deSeqSfChoice="p1Sf", # can be "regSf", "p1Sf"
                     normalizationMethodChoice= "noNorm") # can be "vst", "rlog", "log10", "noNorm"
###*****************************


###*****************************
#Decompose the container
deseq_DataObj=mainData[[1]]
objectName=mainData[[2]]

###*****************************

###*****************************
# To do
# Change level boundries for MgLevel, NaLevel, growthPhase
###*****************************


###*****************************
# Run DESeq2 test For
if(objectName$normalizationMethodChoice=="noNorm")
{
  ###*****************************
  # Do the DeSeq2 test
  test_for="carbonSource"
  DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ",test_for))
  differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
  (res <- DESeq2::results(object = differentialGeneAnalResults, pAdjustMethod ="fdr"))
  mcols(res, use.names=TRUE)  
  DESeq2::summary.DESeqResults(object = res,alpha = 0.05)
  ###*****************************
  
  
  ###*****************************
  # Prepeare dictionaries
  if(objectName$pick_data %in% c("rna","mrna"))
  {
    dictionary=read.csv(file = "../generateDictionary/nameDictionary_RNA_barrick.csv")
    colnames(dictionary)[1]<-"id"
  }
  if(objectName$pick_data %in% c("protein","protein_wo_NA"))
  {
    dictionary=read.csv(file = "../generateDictionary/nameDictionary_Protein.csv")
    colnames(dictionary)[1]<-"id"
  }
  ###*****************************
  
  
  ###*****************************
  #Update objectName
  objectName$pick_data=as.character(objectName$pick_data)
  objectName$test_for=paste0("_",test_for)
  ###*****************************
  
  ###*****************************
  # Generate the metaData 
  metaData<-as.data.frame(colData(deseq_DataObj))
  ###*****************************
  
  
  ###*****************************
  # Prepeare p value data frame
  vs=paste(as.vector(unique(metaData[[test_for]])),collapse = "");
  
  res_df<-as.data.frame(res)
  res_df<-cbind(id=rownames(res_df),res_df)
  res_df%>%
    dplyr::left_join(.,dictionary)%>%
    dplyr::mutate(signChange=sign(log2FoldChange))-> res_df
  
  res_df %>%dplyr::mutate(pick_data=as.vector(objectName$pick_data),
                          growthPhase=as.vector(objectName$growthPhase_names),
                          test_for=test_for,
                          vs=vs)->res_df
  
  res_df %>%
    dplyr::filter(padj<0.05, abs(log2FoldChange)>1)%>%
    dplyr::arrange(padj)->res_df_filtered
  
  genes_0.05=as.vector(res_df_filtered$gene_name)
  genes_0.05=genes_0.05[-grep(pattern = "^[[:blank:]]*$",x = noquote(genes_0.05))]
  ###*****************************
  
  
  ###*****************************
  # SAVE FILES
  if(saveFiles){
    # save genes0.05
    objectName$initial="genes0.05"
    fileName=paste(objectName,collapse = "_")
    write.table(x = genes_0.05, 
                file = paste0("../c_results/",fileName,".csv"),
                row.names = FALSE,
                col.names = "genes",
                quote = FALSE)
    
    # save resDF
    objectName$initial="resDf"
    fileName=paste(objectName,collapse = "_")
    write.csv(x = res_df, 
              file = paste0("../c_results/",fileName,".csv"),
              row.names = FALSE,
              quote = FALSE)
    
    # save metaData
    objectName$initial="metaData"
    fileName=paste(objectName,collapse = "_")
    write.csv(x = metaData, 
              file = paste0("../c_results/",fileName,".csv"),
              row.names = FALSE,
              quote = FALSE)
  }
  ###*****************************
}

if(objectName$normalizationMethodChoice!="noNorm")
{
  res_df<-as.data.frame(assay(deseq_DataObj))
  metaData<-as.data.frame(colData(deseq_DataObj))
  
  ###*****************************
  # SAVE FILES
  if(saveFiles){
    
    # save resDF
    objectName$initial="resDf"
    fileName=paste(objectName,collapse = "_")
    write.csv(x = res_df, 
              file = paste0("../c_results/",fileName,".csv"),
              row.names = FALSE,
              quote = FALSE)
    
    # save metaData
    objectName$initial="metaData"
    fileName=paste(objectName,collapse = "_")
    write.csv(x = metaData, 
              file = paste0("../c_results/",fileName,".csv"),
              row.names = FALSE,
              quote = FALSE)
  }
  ###*****************************
}