###*****************************
# Set Up java
dyn.load('/Library/Java/JavaVirtualMachines/jdk1.8.0_25.jdk/Contents/Home/jre/lib/server/libjvm.dylib')
###*****************************


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
{setwd(paste0("/Users/umut/GitHub/ecoli_multiple_growth_conditions/c_code_change_wrt_variables_RNA&Protein/"))} # mac computer
###*****************************


###*****************************
# DOWNLOAD LIBRARIES
require("Biobase") 
require("dplyr")
require("tidyr")
require("DESeq2")
require("RDAVIDWebService")
require("R.utils")
require(org.Hs.eg.db)
###*****************************


###*****************************
#Load Functions
source("../a_code_dataPreperation_RNA&Protein/data_naming_functions.R")
###*****************************


###*****************************
RUN01=c(dataTypeChoice="mrna", phaseChoice="exponential", phaseBaseChoice="exponential", badDataSetChoice="set00", 
        test_forChoice="carbonSource", test_contrastChoice="glycerol",test_baseChoice="glucose")
RUN02=c(dataTypeChoice="mrna", phaseChoice="exponential", phaseBaseChoice="exponential", badDataSetChoice="set00", 
        test_forChoice="carbonSource", test_contrastChoice="gluconate",test_baseChoice="glucose")
RUN03=c(dataTypeChoice="mrna", phaseChoice="exponential", phaseBaseChoice="exponential", badDataSetChoice="set00", 
        test_forChoice="carbonSource", test_contrastChoice="lactate",test_baseChoice="glucose")

RUN04=c(dataTypeChoice="protein", phaseChoice="exponential", phaseBaseChoice="exponential", badDataSetChoice="set00", 
        test_forChoice="carbonSource", test_contrastChoice="glycerol",test_baseChoice="glucose")
RUN05=c(dataTypeChoice="protein", phaseChoice="exponential", phaseBaseChoice="exponential", badDataSetChoice="set00", 
        test_forChoice="carbonSource", test_contrastChoice="gluconate",test_baseChoice="glucose")
RUN06=c(dataTypeChoice="protein", phaseChoice="exponential", phaseBaseChoice="exponential", badDataSetChoice="set00", 
        test_forChoice="carbonSource", test_contrastChoice="lactate",test_baseChoice="glucose")

RUN07=c(dataTypeChoice="mrna", phaseChoice="stationary", phaseBaseChoice="stationary", badDataSetChoice="set00", 
        test_forChoice="carbonSource", test_contrastChoice="glycerol",test_baseChoice="glucose")
RUN08=c(dataTypeChoice="mrna", phaseChoice="stationary", phaseBaseChoice="stationary", badDataSetChoice="set00", 
        test_forChoice="carbonSource", test_contrastChoice="gluconate",test_baseChoice="glucose")
RUN09=c(dataTypeChoice="mrna", phaseChoice="stationary", phaseBaseChoice="stationary", badDataSetChoice="set00", 
        test_forChoice="carbonSource", test_contrastChoice="lactate",test_baseChoice="glucose")

RUN10=c(dataTypeChoice="protein", phaseChoice="stationary", phaseBaseChoice="stationary", badDataSetChoice="set00", 
        test_forChoice="carbonSource", test_contrastChoice="glycerol",test_baseChoice="glucose")
RUN11=c(dataTypeChoice="protein", phaseChoice="stationary", phaseBaseChoice="stationary", badDataSetChoice="set00", 
        test_forChoice="carbonSource", test_contrastChoice="gluconate",test_baseChoice="glucose")
RUN12=c(dataTypeChoice="protein", phaseChoice="stationary", phaseBaseChoice="stationary", badDataSetChoice="set00", 
        test_forChoice="carbonSource", test_contrastChoice="lactate",test_baseChoice="glucose")


RUN13=c(dataTypeChoice="mrna", phaseChoice="exponential", phaseBaseChoice="exponential", badDataSetChoice="set00", 
        test_forChoice="Mg_mM_Levels", test_contrastChoice="lowMg",test_baseChoice="baseMg")
RUN14=c(dataTypeChoice="mrna", phaseChoice="exponential", phaseBaseChoice="exponential", badDataSetChoice="set00", 
        test_forChoice="Mg_mM_Levels", test_contrastChoice="highMg",test_baseChoice="baseMg")

RUN15=c(dataTypeChoice="protein", phaseChoice="exponential", phaseBaseChoice="exponential", badDataSetChoice="set00", 
        test_forChoice="Mg_mM_Levels", test_contrastChoice="lowMg",test_baseChoice="baseMg")
RUN16=c(dataTypeChoice="protein", phaseChoice="exponential", phaseBaseChoice="exponential", badDataSetChoice="set00", 
        test_forChoice="Mg_mM_Levels", test_contrastChoice="highMg",test_baseChoice="baseMg")

RUN17=c(dataTypeChoice="mrna", phaseChoice="stationary", phaseBaseChoice="stationary", badDataSetChoice="set00", 
        test_forChoice="Mg_mM_Levels", test_contrastChoice="lowMg",test_baseChoice="baseMg")
RUN18=c(dataTypeChoice="mrna", phaseChoice="stationary", phaseBaseChoice="stationary", badDataSetChoice="set00", 
        test_forChoice="Mg_mM_Levels", test_contrastChoice="highMg",test_baseChoice="baseMg")

RUN19=c(dataTypeChoice="protein", phaseChoice="stationary", phaseBaseChoice="stationary", badDataSetChoice="set00", 
        test_forChoice="Mg_mM_Levels", test_contrastChoice="lowMg",test_baseChoice="baseMg")
RUN20=c(dataTypeChoice="protein", phaseChoice="stationary", phaseBaseChoice="stationary", badDataSetChoice="set00", 
        test_forChoice="Mg_mM_Levels", test_contrastChoice="highMg",test_baseChoice="baseMg")


RUN21=c(dataTypeChoice="mrna", phaseChoice="exponential", phaseBaseChoice="exponential", badDataSetChoice="set00", 
        test_forChoice="Na_mM_Levels", test_contrastChoice="highNa",test_baseChoice="baseNa")
RUN22=c(dataTypeChoice="protein", phaseChoice="exponential", phaseBaseChoice="exponential", badDataSetChoice="set00", 
        test_forChoice="Na_mM_Levels", test_contrastChoice="highNa",test_baseChoice="baseNa")
RUN23=c(dataTypeChoice="mrna", phaseChoice="stationary", phaseBaseChoice="stationary", badDataSetChoice="set00", 
        test_forChoice="Na_mM_Levels", test_contrastChoice="highNa",test_baseChoice="baseNa")
RUN24=c(dataTypeChoice="protein", phaseChoice="stationary", phaseBaseChoice="stationary", badDataSetChoice="set00", 
        test_forChoice="Na_mM_Levels", test_contrastChoice="highNa",test_baseChoice="baseNa")

differentCompartisonsList=list(RUN01,RUN02,RUN03,RUN04,RUN05,RUN06,
                               RUN07,RUN08,RUN09,RUN10,RUN11,RUN12,
                               RUN13,RUN14,RUN15,RUN16,RUN17,RUN18,
                               RUN19,RUN20,RUN21,RUN22,RUN23,RUN24)
###*****************************


for(counter01 in 1:length(differentCompartisonsList))
{
  print(paste0("counter :",counter01))
  
  RUN=differentCompartisonsList[[counter01]]
  
  dataTypeChoice=RUN[["dataTypeChoice"]]
  badDataSetChoice=RUN[["badDataSetChoice"]]
  phaseBaseChoice=RUN[["phaseBaseChoice"]]
  phaseChoice=RUN[["phaseChoice"]]
  test_forChoice=RUN[["test_forChoice"]]
  test_baseChoice=RUN[["test_baseChoice"]]
  test_contrastChoice=RUN[["test_contrastChoice"]]
  
  ###*****************************
  ###*****************************
  # Find the csv files that need to be imported
  dataName=name_data(initialValue=c("ez_P0.05Fold2"), # can be c("genes0.05","genes_P0.05Fold2","resDf"),
                     dataType = dataTypeChoice, # can be "rna", "mrna", "protein", "protein_wo_NA"
                     badDataSet = badDataSetChoice, # can be "set00",set01","set02", "set03"
                     # referenceParameters can be a vector like
                     # c("growthPhase", "Mg_mM_Levels", "Na_mM_Levels", "carbonSource", "experiment")
                     referenceParameters=c("growthPhase",
                                           "Mg_mM_Levels", 
                                           "Na_mM_Levels", 
                                           "carbonSource", 
                                           "experiment"),
                     # referenceLevels can be a vector like
                     # c("exponential", "baseMg", "baseNa", "glucose", "glucose_time_course")
                     referenceLevels=c(phaseBaseChoice,
                                       "baseMg", 
                                       "baseNa", 
                                       "glucose", 
                                       "glucose_time_course"),
                     experimentVector = c("allEx"), # can be "Stc","Ytc","Nas","Agr","Ngr","Mgl","Mgh" // "allEx"
                     carbonSourceVector = "SYAN", # can be any sub combination of "SYAN"
                     MgLevelVector = c("allMg"), # can be "lowMg","baseMg","highMg" // "allMg"
                     NaLevelVector = c("allNa"), # can be "baseNa","highNa" // "allNa"
                     # can be "exponential","stationary","late_stationary" // "allPhase"
                     growthPhaseVector = c(phaseChoice), 
                     filterGenes = "noMatchFilter", # can be either "noFilter", or any combination of c("meanFilter", "maxFilter", "sdFilter", "noMatchFilter") 
                     threshold=NA, # the threshold value for "meanFilter", "maxFilter", "sdFilter" can be  c(meanFilter=5,maxFilter=3,sdFilter=7)
                     roundData=TRUE,
                     sumTechnicalReplicates=TRUE,
                     deSeqSfChoice="p1Sf", # can be "regSf", "p1Sf", "noSf"
                     normalizationMethodChoice= "noNorm", # can be "vst", "rlog", "log10", "noNorm"
                     test_for = test_forChoice)  # works only if normalizationMethodChoice == noNorm
  # c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource", "noTest")
  

  objectName=as.data.frame(dataName[[1]])
  ###*****************************
  
  ###*****************************
  #Update objectName
  test_base=test_baseChoice
  test_contrast=test_contrastChoice
  
  test_for=objectName$test_for
  objectName$test_for=paste0("_batchNumberPLUS",gsub("^_","",test_for),"PLUSdoublingTimeMinutes")
  objectName$contrast=paste0("_",test_contrast,"VS",test_base)
  ###*****************************
  
  
  ###****************************
  # Generate Object name
  fileName=paste(objectName,collapse = "_")
  ###*****************************
  
  
  print(fileName)
  ###*****************************
  # Load data
  davidInputData=read.csv(file = paste0("../c_results/" ,fileName,".csv"),header = TRUE)
  davidInputData=as.vector(davidInputData[[1]])
  #davidInputData=read.csv(file="../../../Desktop/test2.csv")
  #davidInputData=as.character(as.vector(davidInputData[[1]]))
  ###*****************************
  
  if(length(davidInputData)!=0)
  {
    ###*****************************
    # Connect to david_d for analyse
    david_d<-DAVIDWebService(email="umut.caglar@utexas.edu",
                             url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
    RDAVIDWebService::setTimeOut(david_d, 90000)
    ###*****************************
    
    
    ###*****************************
    #result<-addList(david_d, davidInputData, idType="OFFICIAL_GENE_SYMBOL", listName="testList", listType="Gene")
    result<-addList(david_d, davidInputData, idType="ENTREZ_GENE_ID", listName="testList", listType="Gene")
    ###*****************************
    
    
    ###*****************************
    # Set species and backround
    selectedSpecie="Escherichia coli str. K-12 substr. MG1655"
    backgroundLocation=grep(selectedSpecie,RDAVIDWebService::getBackgroundListNames(david_d))
    specieLocation=grep(selectedSpecie,RDAVIDWebService::getSpecieNames(david_d))
    setCurrentSpecies(object=david_d, species=specieLocation);setCurrentBackgroundPosition(object=david_d,position=backgroundLocation)
    ###*****************************
    
    
    ###*****************************
    # KEGG TEST
    setAnnotationCategories(david_d, c("KEGG_PATHWAY"))
    objectName$analyzeType="kegg"
    fileName=paste(objectName,collapse = "_")
    
    keggObject<- as.data.frame(getFunctionalAnnotationChart(object=david_d,  threshold=1, count=0L))
    write.csv(x =keggObject ,file = paste0("../c_results/david_results_batchGrowth/",fileName,".csv"))
    ###*****************************
    
    
    ###*****************************
    # MF NEW TEST
    setAnnotationCategories(david_d, c("GOTERM_MF_ALL"))
    objectName$analyzeType="mf_n"
    fileName=paste(objectName,collapse = "_")
    
    mfObject<- as.data.frame(getFunctionalAnnotationChart(object=david_d,  threshold=1, count=0L))
    write.csv(x =mfObject ,file = paste0("../c_results/david_results_batchGrowth/",fileName,".csv"))
    ###*****************************
    
    
    ###*****************************
    # MF OLD
    # Connect to david for analyse
    david<-DAVIDWebService(email="umut.caglar@utexas.edu",
                           url="https://david-d.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
    RDAVIDWebService::setTimeOut(david, 90000)
    ###*****************************
    
    
    ###*****************************
    #result<-addList(david, davidInputData, idType="OFFICIAL_GENE_SYMBOL", listName="testList", listType="Gene")
    result<-addList(david, davidInputData, idType="ENTREZ_GENE_ID", listName="testList", listType="Gene")
    ###*****************************
    
    
    ###*****************************
    # Set species and backround
    selectedSpecie="Escherichia coli"
    backgroundLocation=grep(selectedSpecie,RDAVIDWebService::getBackgroundListNames(david))
    specieLocation=grep(selectedSpecie,RDAVIDWebService::getSpecieNames(david))
    setCurrentSpecies(object=david, species=specieLocation);setCurrentBackgroundPosition(object=david,position=backgroundLocation)
    ###*****************************
    
    
    ###*****************************
    # MF OLD TEST
    setAnnotationCategories(david, c("GOTERM_MF_ALL"))
    objectName$analyzeType="mf_o"
    fileName=paste(objectName,collapse = "_")
    
    mfObject<- as.data.frame(getFunctionalAnnotationChart(object=david,  threshold=1, count=0L))
    write.csv(x =mfObject ,file = paste0("../c_results/david_results_batchGrowth/",fileName,".csv"))
    ###*****************************
    
    
    ###*****************************
    remove(david,david_d)
    ###*****************************
  }
  
}



for(i in list.files(pattern="*.csv$", path = "../c_results/david_results_batchGrowth/"))  # iterate i over .csv files
{
  if(R.utils::countLines(paste0("../c_results/david_results_batchGrowth/",i))==0)
  {
    file.remove(paste0("../c_results/david_results_batchGrowth/",i)) #delete the ones with 0 lines.
  }
  if(file.info(paste0("../c_results/david_results_batchGrowth/",i))$size<=3)
  {
    file.remove(paste0("../c_results/david_results_batchGrowth/",i)) #delete the ones with 3 bytes.
  }
} 





