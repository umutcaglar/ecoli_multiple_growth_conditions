###*****************************
# @Title filter_data function
# The data filtering function that controls sub functions.
name_data<-function(initialValue, # can be "genes0.05"
                    dataType, # can be "rna", "mrna", "protein", "protein_wo_NA"
                    referenceParameters=NA, # c("growthPhase", "Mg_mM_Levels", "Na_mM_Levels", "carbonSource", "experiment")
                    referenceLevels=NA, # c("exponential", "baseMg", "baseNa", "glucose", "glucose_time_course")
                    badDataSet, # can be "set00","set01","set02", "set03"
                    experimentVector, # can be "Stc","Ytc","Nas","Agr","Ngr","Mgl","Mgh" // "allEx"
                    carbonSourceVector, # can be any sub combination of "SYAN"
                    MgLevelVector, # can be "lowMg","baseMg","highMg" // "allMg"
                    NaLevelVector, # can be "baseNa","highNa" // "allNa"
                    growthPhaseVector, # can be "exponential","stationary","late_stationary" // "allPhase"
                    filterGenes, # can be "noFilter", "meanFilter", "maxFilter", "sdFilter" 
                    threshold=NA, # the threshold value for "meanFilter", "maxFilter", "sdFilter"
                    roundData,
                    sumTechnicalReplicates,
                    deSeqSfChoice, # can be "regSf", "p1Sf"
                    normalizationMethodChoice, # can be "vst", "rlog", "log10", "noNorm")
                    test_for)  # works only if normalizationMethodChoice == noNorm
                               # c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
{
  mainData_internal=pick_data_n(dataType=dataType,initialValue)
  mainData_internal=prepeare_data_n(dataInput = mainData_internal,
                                    roundData,
                                    sumTechnicalReplicates)
  for(counter01 in 1:2)
  {
    mainData_internal=select_reference_levels_n(dataInput = mainData_internal, 
                                                referenceParameters,referenceLevels)
    mainData_internal=filter_bad_quality_data_n(dataInput = mainData_internal, 
                                                badDataSet = badDataSet)
    mainData_internal=pick_experiments_n(dataInput = mainData_internal, 
                                         experimentVector = experimentVector)
    mainData_internal=pick_carbonSource_n(dataInput = mainData_internal, 
                                          carbonSourceVector=carbonSourceVector)
    mainData_internal=pick_MgLevel_n(dataInput = mainData_internal, MgLevelVector=MgLevelVector)
    mainData_internal=pick_NaLevel_n(dataInput = mainData_internal, NaLevelVector=NaLevelVector)
    mainData_internal=pick_growthPhase_n(dataInput = mainData_internal, 
                                         growthPhaseVector = growthPhaseVector)
    mainData_internal=filter_rows_n(dataInput=mainData_internal, filterGenes=filterGenes)
  }
  mainData_internal=sizefactors_deseq_n(dataInput=mainData_internal, deSeqSfChoice)
  mainData_internal=normalizeData_n(dataInput=mainData_internal, normalizationMethodChoice)
  if(normalizationMethodChoice=="noNorm")
  {mainData_internal=test_data_function_n(dataInput=mainData_internal, test_for)}
  for(counter01 in 1:length(mainData_internal$objectName))
      {mainData_internal$objectName[[counter01]]=as.character(mainData_internal$objectName[[counter01]])}
  
  return(mainData_internal)
}
###*****************************


###*****************************
# @Title pick_data
# @Param data name could be 
#  "rna" (that installs all rna data including mRNA tRNA and RF)
#  "mrna" (that installs all mrna data)
#  "protein" (that installs all protein data)
#  "protein_wo_NA" (that installs protein data with removed NA rows)
# @Description installs necessary rawData and metaData
pick_data_n<-function(dataType,initialValue){
  if(! dataType %in% c("rna", "mrna", "protein", "protein_wo_NA"))
  {stop("dataType should be one of those rna, mrna, protein, protein_wo")}
  objectName=c(initial=initialValue, pick_data=dataType)
  
  
  if(dataType %in% c("rna","mrna"))
  {metaData=read.csv(paste0("../a_results/","metaRNA.csv"))}
  if(dataType %in% c("protein","protein_wo_NA"))
  {metaData=read.csv(paste0("../a_results/","metaProtein.csv"))}
  
  output=list(objectName=objectName, metaData=metaData)
  return(output)
}
###*****************************

###*****************************
# Prepeare data and sum technical replicates
prepeare_data_n<-function(dataInput=dataType,
                          roundData=TRUE,
                          sumTechnicalReplicates=TRUE)
{
  
  # Seperate data input
  objectName=dataInput$objectName
  metaData=dataInput$metaData
  
  
  # Modify objectName
  if(sumTechnicalReplicates==TRUE){objectName$sumTechRep="trT"}
  if(sumTechnicalReplicates==FALSE){objectName$sumTechRep="trF"}
  
  # Package the results
  dataOutput=list()
  dataOutput$objectName=objectName
  dataOutput$metaData=metaData
  
  
  return(dataOutput)
}
###*****************************

###*****************************
# Remove problematic data sets
# problematic data sets are categorized in 3 different sets
#@param badDataSet can be one of the three options.
# Set 0. No filtering
# set 1 all problematic rna related files
# set 2 filtered problematic rna related files
# set 3 test file to filter our random genes 
filter_bad_quality_data_n<-function(dataInput,badDataSet)
{
  complete_problematicSets=c("set00","set01","set02", "set03")
  if(length(badDataSet)!=1 | !badDataSet %in% complete_problematicSets)
  {stop(paste0("there should be only one set and it should be one of ",
               paste(complete_problematicSets,collapse = " ")))}
  
  if(badDataSet=="set00")
  {badDataVector=c()}
  if(badDataSet=="set01")
  {badDataVector=c("MURI_029","MURI_067","MURI_075","MURI_084",
                   "MURI_086","MURI_091","MURI_136","MURI_138")}
  if(badDataSet=="set02")
  {badDataVector=c("MURI_029","MURI_067","MURI_075","MURI_084",
                   "MURI_136","MURI_138")}
  if(badDataSet=="set03")
  {badDataVector=c("MURI_029","MURI_067","MURI_075","MURI_084",
                   "MURI_136","MURI_138","MURI_131","MURI_119","MURI_107")}
  
  #find related samples
  all_colnames=as.vector(colnames(dataInput$rawData))
  if(badDataSet!="set00")
  {badDataVector_fullName=grep(paste(badDataVector,collapse = "|"),all_colnames,value = TRUE)}
  if(badDataSet=="set00"){badDataVector_fullName=vector()}
  
  
  # Seperate data input
  objectName=dataInput$objectName
  metaData=dataInput$metaData
  
  
  # add information to file name
  objectName$bad_data_set=badDataSet
  objectName=as.data.frame(objectName)
  
  # Package the results
  dataOutput=list()
  dataOutput$objectName=objectName
  dataOutput$metaData=metaData
  
  # return the results
  return(dataOutput)
  
}
###*****************************


###*****************************
# select reference levels
# select the reference levels in metaData by using relevel
select_reference_levels_n<-function(dataInput,
                                    referenceParameters=NA,
                                    referenceLevels=NA)
{
  # if no reference related parameter is mentioned the output is same as input
  if(any(is.na(referenceParameters))&any(is.na(referenceLevels)))
  {
    dataOutput=dataInput
    return(dataOutput)
  }
  
  # controls
  if(length(referenceParameters)!=length(referenceLevels))
  {stop("length of referece parameters should be equal to length of reference levels")}
  if(!all(referenceParameters %in% colnames(dataInput$metaData)))
  {stop("not all reference parameters are in metaData")}
  
  # Seperate data input
  objectName=dataInput$objectName
  metaData=dataInput$metaData
  
  # check for each reference parameter
  for(counter02 in 1: length(referenceParameters))
  {
    referenceParameter=referenceParameters[counter02]
    referenceLevel=referenceLevels[counter02]
    possibleReferenceLevels=as.vector(unique(metaData[[referenceParameter]]))
    if(!referenceLevel %in% possibleReferenceLevels)
    {
      print(referenceParameter)
      print(referenceLevel)
      stop("mentioned reference level is not in possible reference levels")
    }
    
    metaData[[referenceParameter]]=factor(metaData[[referenceParameter]])
    metaData[[referenceParameter]] <- relevel(metaData[[referenceParameter]],
                                              referenceLevel) 
  }
  
  # Package the results
  dataOutput=list()
  dataOutput$objectName=objectName
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
pick_experiments_n<-function(dataInput,experimentVector)
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
  {stop(paste0("possible options for experiments are ",
               paste((completeExperimentVector),collapse = " ")))}
  
  # make the experiment vector for allEx
  if ("allEx" %in% experimentVector){experimentVector=completeExperimentVector}
  # sort the experiment vector wrt completeExperimentVector
  # and convert the experiment vector into the long form
  if (all(experimentVector %in% completeExperimentVector))
  {
    experimentVector=experimentVector[order(match(experimentVector,completeExperimentVector))]
    experimentVectorLong=completeExperimentVectorLong[match(experimentVector,
                                                            completeExperimentVector)]
    
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
  metaData=dataInput$metaData
  
  # Filter meta data
  metaData %>% 
    dplyr::filter(experiment %in% experimentVectorLong) -> metaData
  
  # get the col names from already filtered meta data
  colnames_from_metaData=as.vector(metaData$dataSet)
  
  # add information to file name
  if (all(completeExperimentVector == experimentVector))
  {
    refLevLong=levels(metaData$experiment)[1]
    refLev=completeExperimentVector[which(refLevLong==completeExperimentVectorLong)]
    objectName$experiment_names=paste0(refLev,"AllEx")
  }
  if (!all(completeExperimentVector == experimentVector))
  {
    refLevLong=levels(metaData$experiment)[1]
    refLev=completeExperimentVector[which(refLevLong==completeExperimentVectorLong)]
    experimentVector=c(experimentVector[which(experimentVector==refLev)],
                       experimentVector[-which(experimentVector==refLev)])
    objectName$experiment_names=paste(experimentVector, collapse = "")
  }
  
  objectName=as.data.frame(objectName)
  
  # Package the results
  dataOutput=list()
  dataOutput$objectName=objectName
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
pick_carbonSource_n<-function(dataInput,carbonSourceVector)
{
  completeCarbonSourceVector=c("S","Y","A","N")
  completeCarbonSourceVectorLong = c("glucose", 
                                     "glycerol", 
                                     "lactate", 
                                     "gluconate")
  carbonSourceVector=strsplit(x=carbonSourceVector, split = NULL)
  carbonSourceVector=carbonSourceVector[[1]]
  
  if(!all(carbonSourceVector %in% completeCarbonSourceVector))
  {stop(paste0("possible options for carbon sources are ",
               paste((completeCarbonSourceVector),collapse = " ")))}
  
  # sort the carbon source vector wrt completeCarbonSourceVector
  # and convert the carbon source vector into the long form
  if (all(carbonSourceVector %in% completeCarbonSourceVector))
  {
    carbonSourceVector=carbonSourceVector[order(match(carbonSourceVector,
                                                      completeCarbonSourceVector))]
    carbonSourceVectorLong=completeCarbonSourceVectorLong[match(carbonSourceVector,
                                                                completeCarbonSourceVector)]
    
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
  metaData=dataInput$metaData
  
  # Filter meta data
  metaData %>% 
    dplyr::filter(carbonSource %in% carbonSourceVectorLong) -> metaData
  
  # add information to file name
  refLevLong=levels(metaData$carbonSource)[1]
  refLev=completeCarbonSourceVector[which(refLevLong==completeCarbonSourceVectorLong)]
  carbonSourceVector=c(carbonSourceVector[which(carbonSourceVector==refLev)],
                       carbonSourceVector[-which(carbonSourceVector==refLev)])
  objectName$carbonSource_names=paste(carbonSourceVector, collapse = "")
  objectName=as.data.frame(objectName)
  
  # Package the results
  dataOutput=list()
  dataOutput$objectName=objectName
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
pick_MgLevel_n<-function(dataInput,MgLevelVector)
{
  completeMgLevelVector=c("lowMg","baseMg","highMg")
  if ("allMg" %in% MgLevelVector & length(MgLevelVector)!=1)
  {stop("\"allMg\" is in the MgLevelVector it should be the only element")}
  if(!("allMg" %in% MgLevelVector) & !all(MgLevelVector %in% completeMgLevelVector))
  {stop(paste0("possible options for mg levels are ",
               paste((completeMgLevelVector),collapse = " ")))}
  
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
  metaData=dataInput$metaData
  
  # Filter meta data
  metaData %>% 
    dplyr::filter(Mg_mM_Levels %in% MgLevelVector) -> metaData
  
  # add information to file name
  if (all(completeMgLevelVector == MgLevelVector))
  {
    refLev=levels(metaData$Mg_mM_Levels)[1]
    objectName$MgLevel_names=paste0(refLev,"AllMg")
  }
  if (!all(completeMgLevelVector == MgLevelVector))
  {
    refLev=levels(metaData$Mg_mM_Levels)[1]
    MgLevelVector=c(MgLevelVector[which(MgLevelVector==refLev)],
                    MgLevelVector[-which(MgLevelVector==refLev)])
    objectName$MgLevel_names=paste(MgLevelVector, collapse = "")
  }
  
  objectName=as.data.frame(objectName)
  
  # Package the results
  dataOutput=list()
  dataOutput$objectName=objectName
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
pick_NaLevel_n<-function(dataInput,NaLevelVector)
{
  completeNaLevelVector=c("baseNa","highNa")
  if ("allNa" %in% NaLevelVector & length(NaLevelVector)!=1)
  {stop("\"allNa\" is in the NaLevelVector it should be the only element")}
  if(!("allNa" %in% NaLevelVector) & !all(NaLevelVector %in% completeNaLevelVector))
  {stop(paste0("possible options for Na levels are ",
               paste((completeNaLevelVector),collapse = " ")))}
  
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
  metaData=dataInput$metaData
  
  # Filter meta data
  metaData %>% 
    dplyr::filter(Na_mM_Levels %in% NaLevelVector) -> metaData
  
  # add information to file name
  if (all(completeNaLevelVector == NaLevelVector))
  {
    refLev=levels(metaData$Na_mM_Levels)[1]
    objectName$NaLevel_names=paste0(refLev,"AllNa")
  }
  if (!all(completeNaLevelVector == NaLevelVector))
  {
    refLev=levels(metaData$Na_mM_Levels)[1]
    NaLevelVector=c(NaLevelVector[which(NaLevelVector==refLev)],
                    NaLevelVector[-which(NaLevelVector==refLev)])
    objectName$NaLevel_names=paste(NaLevelVector, collapse = "")
  }
  objectName=as.data.frame(objectName)
  
  # Package the results
  dataOutput=list()
  dataOutput$objectName=objectName
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
pick_growthPhase_n<-function(dataInput,growthPhaseVector)
{
  completeGrowthPhaseVector=c("exponential","stationary","late_stationary")
  if ("allPhase" %in% growthPhaseVector & length(growthPhaseVector)!=1)
  {stop("\"allPhase\" is in the growthPhaseVector it should be the only element")}
  if(!("allPhase" %in% growthPhaseVector) & !all(growthPhaseVector %in% completeGrowthPhaseVector))
  {stop(paste0("possible options for Na levels are ",paste((completeGrowthPhaseVector),collapse = " ")))}
  
  # make the growth phase vector for allPhase
  if ("allPhase" %in% growthPhaseVector){growthPhaseVector=completeGrowthPhaseVector}
  # sort the growth phase vector wrt completeGrowthPhaseVector
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
  metaData=dataInput$metaData
  
  # Filter meta data
  metaData %>% 
    dplyr::filter(growthPhase %in% growthPhaseVector) -> metaData
  
  # add information to file name
  if (all(completeGrowthPhaseVector == growthPhaseVector))
  {
    refLev=levels(metaData$growthPhase)[1]
    if(refLev=="exponential"){refLev="Exp"}
    if(refLev=="stationary"){refLev="Sta"}
    if(refLev=="late_stationary"){refLev="Ltsta"}
    objectName$growthPhase_names=paste0(refLev,"AllPhase")
  }
  if (!all(completeGrowthPhaseVector == growthPhaseVector))
  {
    refLev=levels(metaData$growthPhase)[1]
    growthPhaseVector=c(growthPhaseVector[which(growthPhaseVector==refLev)],
                        growthPhaseVector[-which(growthPhaseVector==refLev)])
    growthPhaseVector=replace(growthPhaseVector, growthPhaseVector=="exponential", "Exp")
    growthPhaseVector=replace(growthPhaseVector, growthPhaseVector=="stationary", "Sta")
    growthPhaseVector=replace(growthPhaseVector, growthPhaseVector=="late_stationary", "Ltsta")
    objectName$growthPhase_names=paste(growthPhaseVector, collapse = "")
  }
  
  objectName=as.data.frame(objectName)
  
  
  # Package the results
  dataOutput=list()
  dataOutput$objectName=objectName
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
filter_rows_n<-function(dataInput, filterGenes, threshold=NA)
{
  # Seperate data input
  objectName=dataInput$objectName
  metaData=dataInput$metaData
  
  # add information to file name
  if(filterGenes=="noFilter")
  {objectName$filter_Name=paste(filterGenes)}
  if(!filterGenes=="noFilter")
  {objectName$filter_Name=paste0(filterGenes, "_", threshold)}
  
  objectName=as.data.frame(objectName)
  
  # Package the results
  dataOutput=list()
  dataOutput$objectName=objectName
  dataOutput$metaData=metaData
  
  # return the results
  return(dataOutput)
}
###*****************************


###*****************************
# @Title sizefactors_deseq
# @Description calculates size factors for data by using DeSeq2
# @Param deSeqSfChoice picks the method for DeSeqSf Calculation
sizefactors_deseq_n<-function(dataInput=mainData_internal, deSeqSfChoice)
{
  if(!deSeqSfChoice %in% c("regSf", "p1Sf"))
  {stop("deSeqSfChoice should be one of regSf p1Sf")}
  
  # Seperate data input
  objectName=dataInput$objectName
  metaData=dataInput$metaData
  
  
  if(deSeqSfChoice=="regSf")
  {
    # modify object name
    objectName$deSeqSfChoice="regSf"
  }
  
  
  if(deSeqSfChoice=="p1Sf")
  {
    
    # modify object name
    objectName$deSeqSfChoice="p1Sf"
  }
  
  deseq_Data_Container=list(objectName=objectName, metaData=metaData)
  return(deseq_Data_Container)
}
###*****************************


###*****************************
# @Title normalizeData
# @Description The function normalize the data by 
# @Param deSeqSfChoice picks the method for DeSeqSf Calculation
normalizeData_n<-function(dataInput=mainData_internal, normalizationMethodChoice)
{
  # decompose the container
  objectName=dataInput[["objectName"]]
  metaData=dataInput[["metaData"]]
  
  if(!length(normalizationMethodChoice) == 1)
  {stop("normalizationMethodChoice should have length 1")}
  
  if(!normalizationMethodChoice %in% c("vst", "rlog" , "log10", "noNorm"))
  {stop("normalizationMethodChoice should be one of vst, rlog, log10, noNorm")}
  
  if(normalizationMethodChoice == "vst")
  {objectName$normalizationMethodChoice="vst"}
  
  if(normalizationMethodChoice == "rlog")
  {objectName$normalizationMethodChoice="rlog"}
  
  if(normalizationMethodChoice == "noNorm")
  {objectName$normalizationMethodChoice="noNorm"}
  
  if(normalizationMethodChoice == "log10")
  {objectName$normalizationMethodChoice="log10"}
  
  # recompose the container
  deseq_Data_Container=list(objectName=objectName, metaData=metaData)
  return(deseq_Data_Container)
}
###*****************************



###*****************************
# @Title testing function
# name the test variables
test_data_function_n<-function(dataInput=mainData_internal, test_for)
{
  # decompose the container
  objectName=dataInput[["objectName"]]
  metaData=dataInput[["metaData"]]
  
  # control test_for
  if(!test_for %in% c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource"))
    {stop("test_for should be in one of \"Mg_mM_Levels\", \"Na_mM_Levels\", \"growthPhase\", \"carbonSource\"")}
  
  # modify object name
  objectName$test_for=paste0("_",test_for)
  
  # recompose the container
  deseq_Data_Container=list(objectName=objectName)
  return(deseq_Data_Container)
}