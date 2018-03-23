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
require(org.Hs.eg.db)
###*****************************


###*****************************
#Load Functions
source("../a_code_dataPreperation_RNA&Protein/data_naming_functions.R")
###*****************************


###*****************************
# Find the csv files that need to be imported
dataName=name_data(initialValue=c("ez_P0.05Fold2"), # can be c("genes0.05","genes_P0.05Fold2","resDf")
                   dataType = "mrna", # can be "rna", "mrna", "protein", "protein_wo_NA"
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
                   growthPhaseVector = c("exponential"), # can be "exponential","stationary","late_stationary" // "allPhase"
                   filterGenes = "noMatchFilter", # can be "noFilter", "meanFilter", "maxFilter", "sdFilter"
                   threshold=NA, # the threshold value for "meanFilter", "maxFilter", "sdFilter"
                   roundData=TRUE,
                   sumTechnicalReplicates=TRUE,
                   deSeqSfChoice="p1Sf", # can be "regSf", "p1Sf", "noSf"
                   normalizationMethodChoice= "noNorm", # can be "vst", "rlog", "log10", "noNorm"
                   test_for = "carbonSource")  # works only if normalizationMethodChoice == noNorm
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource", "noTest")

objectName=as.data.frame(dataName[[1]])
###*****************************



###*****************************
#Update objectName
test_base="glucose"
test_contrast="glycerol"

test_for=objectName$test_for
objectName$test_for=paste0("_batchNumberPLUS",gsub("_","",test_for))
objectName$contrast=paste0("_",test_contrast,"VS",test_base)
###*****************************


###****************************
# Generate Object name
fileName=paste(objectName,collapse = "_")
###*****************************


###*****************************
# Load data
davidInputData=read.csv(file = paste0("../c_results/" ,fileName,".csv"),header = TRUE)
davidInputData=as.vector(davidInputData[[1]])
#davidInputData=read.csv(file="../../../Desktop/test2.csv")
#davidInputData=as.character(as.vector(davidInputData[[1]]))
###*****************************

###*****************************
# Connect to david for analyse
david_d<-DAVIDWebService(email="umut.caglar@utexas.edu",
                       url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
#result<-addList(david_d, davidInputData, idType="OFFICIAL_GENE_SYMBOL", listName="testList", listType="Gene")
result<-addList(david_d, davidInputData, idType="ENTREZ_GENE_ID", listName="testList", listType="Gene")
###*****************************


###*****************************
# Set species and backround
selectedSpecie="Escherichia coli str. K-12 substr. MG1655"
backgroundLocation=grep(selectedSpecie,RDAVIDWebService::getBackgroundListNames(david_d))
specieLocation=grep(selectedSpecie,RDAVIDWebService::getSpecieNames(david_d))
setCurrentSpecies(object=david_d, species=specieLocation);
setCurrentBackgroundPosition(object=david_d,position=backgroundLocation)
###*****************************


###*****************************
# KEGG TEST
setAnnotationCategories(david_d, c("KEGG_PATHWAY"))
objectName$analyzeType="kegg"
fileName=paste(objectName,collapse = "_")

keggObject<- as.data.frame(getFunctionalAnnotationChart(object=david_d,  threshold=1, count=0L))
write.csv(x =keggObject ,file = paste0("../c_results/david_results_batch/",fileName,".csv"))
###*****************************


###*****************************
# MF NEW TEST
setAnnotationCategories(david_d, c("GOTERM_MF_ALL"))
objectName$analyzeType="mf_n"
fileName=paste(objectName,collapse = "_")

mfObject<- as.data.frame(getFunctionalAnnotationChart(object=david_d,  threshold=1, count=0L))
write.csv(x =mfObject ,file = paste0("../c_results/david_results_batch/",fileName,".csv"))
###*****************************
###*****************************


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
write.csv(x =mfObject ,file = paste0("../c_results/david_results_batch/",fileName,".csv"))
###*****************************
