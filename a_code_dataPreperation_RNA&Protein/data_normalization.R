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
{setwd(paste0("/Users/umut/GitHub/ecoli_multiple_growth_conditions/a_code_dataPreperation_RNA&Protein/"))} # mac computer
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
source("../a_code_dataPreperation_RNA&Protein/data_filter_normalization_functions.R")
###*****************************


###*****************************
saveFiles=TRUE
runDeSeqForDifExp=FALSE
# The data filtering function that controls sub functions.
mainData=filter_data(dataType = "mrna", # can be "rna", "mrna", "protein", "protein_wo_NA"
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
                     # can be "exponential","stationary","late_stationary" // "allPhase"
                     growthPhaseVector = c("allPhase"), 
                     filterGenes = "noFilter", # can be "noFilter", "meanFilter", "maxFilter", "sdFilter" 
                     threshold=NA, # the threshold value for "meanFilter", "maxFilter", "sdFilter"
                     roundData=TRUE,
                     sumTechnicalReplicates=TRUE,
                     deSeqSfChoice="p1Sf", # can be "regSf", "p1Sf", "noSf"
                     normalizationMethodChoice= "vst") # can be "vst", "rlog", "log10", "noNorm"

# DeSeq2 parameters
test_for="carbonSource"
test_base="glucose"
test_contrast="glycerol"
###*****************************


###*****************************
#Decompose the container
deseq_DataObj=mainData[[1]]
objectName=mainData[[2]]
###*****************************


###*****************************
# Run DESeq2 test For
if(objectName$normalizationMethodChoice=="noNorm" & runDeSeqForDifExp)
{
  ###*****************************
  # Do the DeSeq2 test
  # c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
  DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ","batchNumber + ",test_for))
  differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
  res <- DESeq2::results(object = differentialGeneAnalResults, 
                         pAdjustMethod ="fdr",
                         contrast = c(test_for,test_contrast,test_base))
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
  objectName$test_for=paste0("_batchNumber+",test_for)
  objectName$contrast=paste0("_",test_contrast,"VS",test_base)
  ###*****************************
  
  ###*****************************
  # Generate the metaData 
  metaData<-as.data.frame(colData(deseq_DataObj))
  ###*****************************
  
  
  ###*****************************
  # Prepeare p value data frame
  res_df<-as.data.frame(res)
  res_df<-cbind(id=rownames(res_df),res_df)
  res_df%>%
    dplyr::left_join(.,dictionary)%>%
    dplyr::mutate(signChange=sign(log2FoldChange))-> res_df
  
  res_df %>%dplyr::mutate(pick_data=as.vector(objectName$pick_data),
                          growthPhase=as.vector(objectName$growthPhase_names),
                          test_for=test_for,
                          contrast=test_contrast,
                          base=test_base)->res_df
  
  res_df %>%
    dplyr::filter(padj<0.05, abs(log2FoldChange)>1)%>%
    dplyr::arrange(padj)->res_df_filtered
  
  genes_0.05=as.vector(res_df_filtered$gene_name)
  listOfEmptyCells<-grep(pattern = "^[[:blank:]]*$",x = noquote(genes_0.05))
  listOfFilledCells=setdiff(seq(1,length(genes_0.05)),listOfEmptyCells)
  genes_0.05<-genes_0.05[listOfFilledCells]
  ###*****************************
  

  ###*****************************
  # SAVE FILES
  if(saveFiles){
    # save genes0.05
    objectName$initial="genes_P0.05Fold2"
    fileName=paste(objectName,collapse = "_")
    write.table(x = genes_0.05, 
                file = paste0("../a_results/",fileName,".csv"),
                row.names = FALSE,
                col.names = "genes",
                quote = FALSE)
    
    # save resDF
    objectName$initial="resDf"
    fileName=paste(objectName,collapse = "_")
    write.csv(x = res_df, 
              file = paste0("../a_results/",fileName,".csv"),
              row.names = TRUE,
              quote = FALSE)
    
    # save metaData
    objectName$initial="metaData"
    fileName=paste(objectName,collapse = "_")
    write.csv(x = metaData, 
              file = paste0("../a_results/",fileName,".csv"),
              row.names = FALSE,
              quote = FALSE)
  }
  ###*****************************
}

if(objectName$normalizationMethodChoice!="noNorm" & runDeSeqForDifExp)
{
  stop("TO RUN DeSeq2 one needs to use RAW data; please chose \n \"normalizationMethodChoice==\"noNorm\" if runDeSeqForDifExp==TRUE")
}

if(runDeSeqForDifExp==FALSE)
{
  if(objectName$normalizationMethodChoice=="noNorm")
  {
    res_df<-counts(deseq_DataObj, normalized=TRUE)
  }
  if(objectName$normalizationMethodChoice!="noNorm")
  {
    res_df<-as.data.frame(assay(deseq_DataObj))
  }
  
  metaData<-as.data.frame(colData(deseq_DataObj))
  
  ###*****************************
  # SAVE FILES
  if(saveFiles){
    
    # save resDF
    objectName$initial="resDf"
    fileName=paste(objectName,collapse = "_")
    write.csv(x = res_df, 
              file = paste0("../a_results/",fileName,".csv"),
              row.names = TRUE,
              quote = FALSE)
    
    # save metaData
    objectName$initial="metaData"
    fileName=paste(objectName,collapse = "_")
    write.csv(x = metaData, 
              file = paste0("../a_results/",fileName,".csv"),
              row.names = FALSE,
              quote = FALSE)
  }
  ###*****************************
}