###*****************************
# @Title pick_data
# @Param data name could be 
#  "rna" (that installs all rna data including mRNA tRNA and RF)
#  "mrna" (that installs all mrna data)
#  "protein" (that installs all protein data)
#  "protein_wo_NA" (that installs protein data with removed NA rows)
# @Description installs necessary rawData and metaData
pick_data<-function(data_name){
  if(! data_name %in% c("rna", "mrna", "protein", "protein_wo_NA"))
  {stop("data_name should be one of those rna, mrna, protein, protein_wo")}
  objectName=c(initial="unnormalzied", pick_data=data_name)
  
  if(data_name=="rna")
  {rawData=read.csv(file=paste0("../a_results/","rnaMatrix_RNA.csv"))}
  if(data_name=="mrna")
  {rawData=read.csv(file=paste0("../a_results/","rnaMatrix_mRNA.csv"))}
  if(data_name=="protein")
  {rawData=read.csv(file=paste0("../a_results/","proteinMatrix.csv"))}
  if(data_name=="protein_wo_NA")
  {rawData=read.csv(file=paste0("../a_results/","proteinMatrix_wo_NA.csv"))}
  
  if(data_name %in% c("rna","mrna"))
  {metaData=read.csv(paste0("../a_results/","metaRNA.csv"))}
  if(data_name %in% c("protein","protein_wo_NA"))
  {metaData=read.csv(paste0("../a_results/","metaProtein.csv"))}
  
  output=list(objectName=objectName, rawData=rawData, metaData=metaData)
  return(output)
}
###*****************************


###*****************************
# Remove problematic data sets
# problematic data sets are categorized in 3 different sets
#@param problematic_set can be one of the three options.
# Set 0. No filtering
# set 1 all problematic rna related files
# set 2 filtered problematic rna related files
filter_bad_quality_data<-function(dataInput,problematic_set)
{
  complete_problematicSets=c("set00","set01","set02")
  if(length(problematic_set)!=1 | !problematic_set %in% complete_problematicSets)
  {stop(paste0("there should be only one set and it should be one of ",
               paste(complete_problematicSets,collapse = " ")))}
  
  if(problematic_set=="set00")
  {badDataSet=c()}
  if(problematic_set=="set01")
  {badDataSet=c("MURI_029","MURI_067","MURI_075","MURI_084",
                "MURI_086","MURI_091","MURI_136","MURI_138")}
  if(problematic_set=="set02")
  {badDataSet=c("MURI_029","MURI_067","MURI_075","MURI_084",
                "MURI_136","MURI_138")}
  
  #find related samples
  all_colnames=as.vector(colnames(dataInput$rawData))
  badDataSet_fullName=grep(paste(badDataSet,collapse = "|"),all_colnames,value = TRUE)
  
  # Seperate data input
  objectName=dataInput$objectName
  rawData=dataInput$rawData
  metaData=dataInput$metaData
  
  # find remaining columns in raw data
  if(length(badDataSet_fullName)!=0){
    rawData %>%
      dplyr::select_(.dots=paste0("-",badDataSet_fullName))->rawData
    
    # find remaining rows in meta data
    metaData %>%
      dplyr::filter(!dataSet %in% badDataSet_fullName) ->metaData
  }
  
  # add information to file name
  objectName$bad_data_set=problematic_set
  objectName=as.data.frame(objectName)
  
  # Package the results
  dataOutput=list()
  dataOutput$objectName=objectName
  dataOutput$rawData=rawData
  dataOutput$metaData=metaData
  
  # return the results
  return(dataOutput)
  
}
###*****************************


###*****************************
# @Title pick_experiment
# @Param experimentVector Experiment vector can be composed of abbrevations of experiment names
#   "Stc"="glucose_time_course", 
#   "Ytc"="glycerol_time_course", 
#   "Nas"="NaCl_stress", 
#   "Agr"="lactate_growth", 
#   "Ngr"="gluconate_growth", 
#   "Mgl"="MgSO4_stress_low", 
#   "Mgh"="MgSO4_stress_high"
# @Description it filters the experiments from raw data and meta data 
#      also adds remainng experiments to objectName
pick_experiments<-function(dataInput,experimentVector)
{
  completeExperimentVector=c("Stc","Ytc","Nas","Agr","Ngr","Mgl","Mgh")
  completeExperimentVectorLong = c("glucose_time_course", 
                                   "glycerol_time_course", 
                                   "NaCl_stress", 
                                   "lactate_growth", 
                                   "gluconate_growth", 
                                   "MgSO4_stress_low", 
                                   "MgSO4_stress_high")
  
  if ("allEx" %in% experimentVector & length(experimentVector)!=1)
  {stop("\"allEx\" is in the experimentVector it should be the only element")}
  if(!("allEx" %in% experimentVector) & !all(experimentVector %in% completeExperimentVector))
  {stop(paste0("possible options for experiments are ",paste((completeExperimentVector),collapse = " ")))}
  
  # make the experiment vector for allEx
  if ("allEx" %in% experimentVector){experimentVector=completeExperimentVector}
  # sort the experiment vector wrt completeExperimentVector
  # and convert the experiment vector into the long form
  if (all(experimentVector %in% completeExperimentVector))
  {
    experimentVector=experimentVector[order(match(experimentVector,completeExperimentVector))]
    experimentVectorLong=completeExperimentVectorLong[match(experimentVector,completeExperimentVector)]
    
    initialConditions=as.vector(unique(dataInput$metaData$experiment))
    difference=setdiff(experimentVectorLong,initialConditions)
    if(length(difference!=0))
    {warning(paste0("the following required conditions does not exist in data ",
                    paste(difference, collapse = " ")))}
    experimentVectorLong=intersect(experimentVectorLong,initialConditions)
    experimentVector=
      completeExperimentVector[match(experimentVectorLong,completeExperimentVectorLong)]
  }
  
  # at that point I have the experiment names that I can use to filter data
  
  # Seperate data input
  objectName=dataInput$objectName
  rawData=dataInput$rawData
  metaData=dataInput$metaData
  
  # Filter meta data
  metaData %>% 
    dplyr::filter(experiment %in% experimentVectorLong) -> metaData
  
  # get column names that do not start with MURI in raw data
  colnames_rawData=as.vector(colnames(mainData$rawData))
  colnames_wo_MURI=colnames_rawData[-grep("MURI*",colnames_rawData)]
  # get the col names from already filtered meta data
  colnames_from_metaData=as.vector(metaData$dataSet)
  # combine those two to obtain new set of column names
  colnames_rawData_new=c(colnames_wo_MURI,colnames_from_metaData)
  
  # filter raw data
  rawData %>%
    dplyr::select_(.dots=colnames_rawData_new)->rawData
  
  # add information to file name
  if (all(completeExperimentVector == experimentVector))
  {objectName$experiment_names=c("allEx")}
  if (!all(completeExperimentVector == experimentVector))
  {objectName$experiment_names=paste(experimentVector, collapse = "_")}
  
  objectName=as.data.frame(objectName)
  
  # Package the results
  dataOutput=list()
  dataOutput$objectName=objectName
  dataOutput$rawData=rawData
  dataOutput$metaData=metaData
  
  # return the results
  return(dataOutput)
}
###*****************************


###*****************************
# @Title pick_carbonSource
# @param carbonSourceVector Carbon source vector can be composed of abbrevations of carbon sources
#   "S"="glucose_time_course", 
#   "Y"="glycerol_time_course", 
#   "A="lactate_growth", 
#   "N"="gluconate_growth" 
# @Description it filters the carbon sources from raw data and meta data 
#      also adds remainig carbon sources to objectName
pick_carbonSource<-function(dataInput,carbonSourceVector)
{
  completeCarbonSourceVector=c("S","Y","A","N")
  completeCarbonSourceVectorLong = c("glucose", 
                                     "glycerol", 
                                     "lactate", 
                                     "gluconate")
  carbonSourceVector=strsplit(x=carbonSourceVector, split = NULL)
  carbonSourceVector=carbonSourceVector[[1]]
  
  if(!all(carbonSourceVector %in% completeCarbonSourceVector))
  {stop(paste0("possible options for carbon sources are ",paste((completeCarbonSourceVector),collapse = " ")))}
  
  # sort the carbon source vector wrt completeCarbonSourceVector
  # and convert the carbon source vector into the long form
  if (all(carbonSourceVector %in% completeCarbonSourceVector))
  {
    carbonSourceVector=carbonSourceVector[order(match(carbonSourceVector,completeCarbonSourceVector))]
    carbonSourceVectorLong=completeCarbonSourceVectorLong[match(carbonSourceVector,completeCarbonSourceVector)]
    
    initialConditions=as.vector(unique(dataInput$metaData$carbonSource))
    difference=setdiff(carbonSourceVectorLong,initialConditions)
    if(length(difference!=0))
    {warning(paste0("the following required conditions does not exist in data ",
                    paste(difference, collapse = " ")))}
    carbonSourceVectorLong=intersect(carbonSourceVectorLong,initialConditions)
    carbonSourceVector=
      completeCarbonSourceVector[match(carbonSourceVectorLong,completeCarbonSourceVectorLong)]
  }
  
  # at that point I have the carbon source names that I can use to filter data
  
  # Seperate data input
  objectName=dataInput$objectName
  rawData=dataInput$rawData
  metaData=dataInput$metaData
  
  # Filter meta data
  metaData %>% 
    dplyr::filter(carbonSource %in% carbonSourceVectorLong) -> metaData
  
  # get column names that do not start with MURI in raw data
  colnames_rawData=as.vector(colnames(mainData$rawData))
  colnames_wo_MURI=colnames_rawData[-grep("MURI*",colnames_rawData)]
  # get the col names from already filtered meta data
  colnames_from_metaData=as.vector(metaData$dataSet)
  # combine those two to obtain new set of column names
  colnames_rawData_new=c(colnames_wo_MURI,colnames_from_metaData)
  
  # filter raw data
  rawData %>%
    dplyr::select_(.dots=colnames_rawData_new)->rawData
  
  # add information to file name
  objectName$carbonSource_names=paste(carbonSourceVector, collapse = "")
  objectName=as.data.frame(objectName)
  
  # Package the results
  dataOutput=list()
  dataOutput$objectName=objectName
  dataOutput$rawData=rawData
  dataOutput$metaData=metaData
  
  # return the results
  return(dataOutput)
}
###*****************************


###*****************************
# @Title pick_MgLevel
# @Param MgLevelVector MgLevel vector can be composed of abbrevations of MgLevel names
#   "lowMg"="low Mg levels", 
#   "baseMg"="midium Mg levels", 
#   "highMg"="high Mg levels"
# @Description it filters the Mg levels from raw data and meta data 
#      also adds remainng Mg levels to objectName
pick_MgLevel<-function(dataInput,MgLevelVector)
{
  completeMgLevelVector=c("lowMg","baseMg","highMg")
  if ("allMg" %in% MgLevelVector & length(MgLevelVector)!=1)
  {stop("\"allMg\" is in the MgLevelVector it should be the only element")}
  if(!("allMg" %in% MgLevelVector) & !all(MgLevelVector %in% completeMgLevelVector))
  {stop(paste0("possible options for mg levels are ",paste((completeMgLevelVector),collapse = " ")))}
  
  # make the experiment vector for allEx
  if ("allMg" %in% MgLevelVector){MgLevelVector=completeMgLevelVector}
  # sort the experiment vector wrt completeExperimentVector
  if (all(MgLevelVector %in% completeMgLevelVector))
  {
    MgLevelVector=MgLevelVector[order(match(MgLevelVector,completeMgLevelVector))]
    initialConditions=as.vector(unique(dataInput$metaData$Mg_mM_Levels))
    difference=setdiff(MgLevelVector,initialConditions)
    if(length(difference!=0))
    {warning(paste0("the following required conditions does not exist in data ",
                    paste(difference, collapse = " ")))}
    MgLevelVector=intersect(MgLevelVector,initialConditions)
  }
  
  # at that point I have the MgLevel names that I can use to filter data
  
  # Seperate data input
  objectName=dataInput$objectName
  rawData=dataInput$rawData
  metaData=dataInput$metaData
  
  # Filter meta data
  metaData %>% 
    dplyr::filter(Mg_mM_Levels %in% MgLevelVector) -> metaData
  
  # get column names that do not start with MURI in raw data
  colnames_rawData=as.vector(colnames(mainData$rawData))
  colnames_wo_MURI=colnames_rawData[-grep("MURI*",colnames_rawData)]
  # get the col names from already filtered meta data
  colnames_from_metaData=as.vector(metaData$dataSet)
  # combine those two to obtain new set of column names
  colnames_rawData_new=c(colnames_wo_MURI,colnames_from_metaData)
  
  # filter raw data
  rawData %>%
    dplyr::select_(.dots=colnames_rawData_new)->rawData
  
  # add information to file name
  if (all(completeMgLevelVector == MgLevelVector))
  {objectName$MgLevel_names=c("allMg")}
  if (!all(completeMgLevelVector == MgLevelVector))
  {objectName$MgLevel_names=paste(MgLevelVector, collapse = "_")}
  
  objectName=as.data.frame(objectName)
  
  # Package the results
  dataOutput=list()
  dataOutput$objectName=objectName
  dataOutput$rawData=rawData
  dataOutput$metaData=metaData
  
  # return the results
  return(dataOutput)
}
###*****************************


###*****************************
# @Title pick_NaLevel
# @Param NaLevelVector NaLevel vector can be composed of abbrevations of NaLevel names
#   "baseNa"="midium Na levels", 
#   "highNa"="high Na levels"
# @Description it filters the Na levels from raw data and meta data 
#      also adds remainng Na levels to objectName
pick_NaLevel<-function(dataInput,NaLevelVector)
{
  completeNaLevelVector=c("baseNa","highNa")
  if ("allNa" %in% NaLevelVector & length(NaLevelVector)!=1)
  {stop("\"allNa\" is in the NaLevelVector it should be the only element")}
  if(!("allNa" %in% NaLevelVector) & !all(NaLevelVector %in% completeNaLevelVector))
  {stop(paste0("possible options for Na levels are ",paste((completeNaLevelVector),collapse = " ")))}
  
  # make the experiment vector for allEx
  if ("allNa" %in% NaLevelVector){NaLevelVector=completeNaLevelVector}
  # sort the experiment vector wrt completeExperimentVector
  if (all(NaLevelVector %in% completeNaLevelVector))
  {
    NaLevelVector=NaLevelVector[order(match(NaLevelVector,completeNaLevelVector))]
    initialConditions=as.vector(unique(dataInput$metaData$Na_mM_Levels))
    difference=setdiff(NaLevelVector,initialConditions)
    if(length(difference!=0))
    {warning(paste0("the following required conditions does not exist in data ",
                    paste(difference, collapse = " ")))}
    NaLevelVector=intersect(NaLevelVector,initialConditions)
  }
  
  # at that point I have the NaLevel names that I can use to filter data
  
  # Seperate data input
  objectName=dataInput$objectName
  rawData=dataInput$rawData
  metaData=dataInput$metaData
  
  # Filter meta data
  metaData %>% 
    dplyr::filter(Na_mM_Levels %in% NaLevelVector) -> metaData
  
  # get column names that do not start with MURI in raw data
  colnames_rawData=as.vector(colnames(mainData$rawData))
  colnames_wo_MURI=colnames_rawData[-grep("MURI*",colnames_rawData)]
  # get the col names from already filtered meta data
  colnames_from_metaData=as.vector(metaData$dataSet)
  # combine those two to obtain new set of column names
  colnames_rawData_new=c(colnames_wo_MURI,colnames_from_metaData)
  
  # filter raw data
  rawData %>%
    dplyr::select_(.dots=colnames_rawData_new)->rawData
  
  # add information to file name
  if (all(completeNaLevelVector == NaLevelVector))
  {objectName$NaLevel_names=c("allNa")}
  if (!all(completeNaLevelVector == NaLevelVector))
  {objectName$NaLevel_names=paste(NaLevelVector, collapse = "_")}
  objectName=as.data.frame(objectName)
  
  # Package the results
  dataOutput=list()
  dataOutput$objectName=objectName
  dataOutput$rawData=rawData
  dataOutput$metaData=metaData
  
  # return the results
  return(dataOutput)
}
###*****************************


###*****************************
# @Title pick_growthPhase
# @Param growthPhaseVector growthPhase vector can be composed of abbrevations of growthPhase names
#   "exponential" = "exponential phase"
#   "stationary" = "stationary phase", 
#   "late_stationary" = "late stationary phase"
# @Description it filters the growth phase from raw data and meta data 
#      also adds remainng chosen growth phase to objectName
pick_growthPhase<-function(dataInput,growthPhaseVector)
{
  completeGrowthPhaseVector=c("exponential","stationary","late_stationary")
  if ("allPhase" %in% growthPhaseVector & length(growthPhaseVector)!=1)
  {stop("\"allPhase\" is in the growthPhaseVector it should be the only element")}
  if(!("allPhase" %in% growthPhaseVector) & !all(growthPhaseVector %in% completeGrowthPhaseVector))
  {stop(paste0("possible options for Na levels are ",paste((completeGrowthPhaseVector),collapse = " ")))}
  
  # make the experiment vector for allEx
  if ("allPhase" %in% growthPhaseVector){growthPhaseVector=completeGrowthPhaseVector}
  # sort the experiment vector wrt completeExperimentVector
  if (all(growthPhaseVector %in% completeGrowthPhaseVector))
  {
    growthPhaseVector=growthPhaseVector[order(match(growthPhaseVector,completeGrowthPhaseVector))]
    initialConditions=as.vector(unique(dataInput$metaData$growthPhase))
    difference=setdiff(growthPhaseVector,initialConditions)
    if(length(difference!=0))
    {warning(paste0("the following required conditions does not exist in data ",
                    paste(difference, collapse = " ")))}
    growthPhaseVector=intersect(growthPhaseVector,initialConditions)
  }
  
  # at that point I have the growth phase names that I can use to filter data
  
  # Seperate data input
  objectName=dataInput$objectName
  rawData=dataInput$rawData
  metaData=dataInput$metaData
  
  # Filter meta data
  metaData %>% 
    dplyr::filter(growthPhase %in% growthPhaseVector) -> metaData
  
  # get column names that do not start with MURI in raw data
  colnames_rawData=as.vector(colnames(mainData$rawData))
  colnames_wo_MURI=colnames_rawData[-grep("MURI*",colnames_rawData)]
  # get the col names from already filtered meta data
  colnames_from_metaData=as.vector(metaData$dataSet)
  # combine those two to obtain new set of column names
  colnames_rawData_new=c(colnames_wo_MURI,colnames_from_metaData)
  
  # filter raw data
  rawData %>%
    dplyr::select_(.dots=colnames_rawData_new)->rawData
  
  # add information to file name
  if (all(completeGrowthPhaseVector == growthPhaseVector))
  {objectName$growthPhase_names=c("allPhase")}
  if (!all(completeGrowthPhaseVector == growthPhaseVector))
  {objectName$growthPhase_names=paste(growthPhaseVector, collapse = "_")}
  
  objectName=as.data.frame(objectName)
  
  
  # Package the results
  dataOutput=list()
  dataOutput$objectName=objectName
  dataOutput$rawData=rawData
  dataOutput$metaData=metaData
  
  # return the results
  return(dataOutput)
}
###*****************************



###*****************************
# @Title filter_rows
# @Param fitering method is the way tha we filter out the mon-significant rows
#   methods are
#   "noFilter" (basically keeps all the rows)
#   "meanFilter" (filter out the rows whose mean is below threshold)
#   "maxFilter" (filter out the rows whose max is below threshold)
#   "sdFilter" (filter out the rows whose standard deviation is below threshold)
# @param threshold the threshold related with relevant model. Default is 0.
# @Description The function filter out the meaningless rows from data
filter_rows<-function(dataInput, filtering_method, threshold=NA)
{
  # Seperate data input
  objectName=dataInput$objectName
  rawData=dataInput$rawData
  metaData=dataInput$metaData
  
  all_colnames=as.vector(colnames(rawData))
  selected_colnames=grep(pattern = "MURI*",x = all_colnames,value = TRUE)
  unselected_colnames=all_colnames[-grep(pattern = "MURI*",x = all_colnames)]
  
  rawData %>%
    tidyr::gather_(key = "dataSet", 
                   value = "numRead", 
                   selected_colnames)%>% 
    dplyr::group_by_(unselected_colnames) %>% 
    dplyr::summarize(meanValue=mean(numRead),
                     sdValue=sd(numRead),
                     maxValue=max(numRead))->rawData_summarize
  
  dplyr::left_join(rawData,rawData_summarize)->rawData
  
  
  
  # Do the filtering
  if(filtering_method=="noFilter"){rawData = rawData }
  if(filtering_method=="meanFilter")
  {rawData %>%dplyr::filter(meanValue>threshold)->rawData}
  
  if(filtering_method=="sdFilter")
  {rawData %>% dplyr::filter(sdValue>threshold)->rawData}
  
  if(filtering_method=="maxFilter")
  {rawData %>% dplyr::filter(maxValue>threshold)->rawData}
  
  # get rid of additional columns
  rawData %>%
    dplyr::select(-meanValue, -maxValue, -sdValue)->rawData
  
  # add information to file name
  if(filtering_method=="noFilter")
  {objectName$filter_Name=paste(filtering_method)}
  if(!filtering_method=="noFilter")
  {objectName$filter_Name=paste0(filtering_method, "_", threshold)}
  
  objectName=as.data.frame(objectName)
  
  # Package the results
  dataOutput=list()
  dataOutput$objectName=objectName
  dataOutput$rawData=rawData
  dataOutput$metaData=metaData
  
  # return the results
  return(dataOutput)
}
###*****************************