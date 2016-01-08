###*****************************
# CALCULATION OF SIZE FACTORS (BY USING P1 METHOD) AND CALCULATION OF VARIANCE STABILIZING TRANSFORMATIONS

# FUNCTIONS

# Generate DeSeq Object (With P1) function
generate_p1_object=function(unnormalized_sampleData, meta_condition){
  
  # Generate "p1 Data Object"
  unnormalized_sampleData_p1=unnormalized_sampleData+1;
  rnaObject_p1=DESeqDataSetFromMatrix(countData = unnormalized_sampleData_p1,
                                      colData = meta_condition,design = ~ uniqueCondition)
  levelList<-as.character(unique(meta_condition)[["uniqueCondition"]]);
  colData(rnaObject_p1)$uniqueCondition=factor(colData(rnaObject_p1)$uniqueCondition, levels=levelList)
  
  # Calculate Size Factors
  rnaObject_p1=estimateSizeFactors(rnaObject_p1)
  sizeFactors_p1=sizeFactors(rnaObject_p1)
  
  # Generate "regular Data Object"
  rnaObject=DESeqDataSetFromMatrix(countData = unnormalized_sampleData, colData = meta_condition,design = ~ uniqueCondition)
  levelList<-as.character(unique(meta_condition)[["uniqueCondition"]]);
  colData(rnaObject)$uniqueCondition=factor(colData(rnaObject)$uniqueCondition, levels=levelList)
  
  # Import Size Factors From "p1 Data Object"
  sizeFactors(rnaObject) <- sizeFactors_p1
  
  return(rnaObject)
}


# Generate DeSeq Object (With P6) function
generate_p6_object=function(unnormalized_sampleData, meta_condition){
  
  # Generate "p1 Data Object"
  unnormalized_sampleData_p6=unnormalized_sampleData+6;
  rnaObject_p6=DESeqDataSetFromMatrix(countData = unnormalized_sampleData_p6,colData = meta_condition,design = ~ uniqueCondition)
  levelList<-as.character(unique(meta_condition)[["uniqueCondition"]]);
  colData(rnaObject_p6)$uniqueCondition=factor(colData(rnaObject_p6)$uniqueCondition, levels=levelList)
  
  # Calculate Size Factors
  rnaObject_p6=estimateSizeFactors(rnaObject_p6)
  sizeFactors_p6=sizeFactors(rnaObject_p6)
  
  # Generate "regular Data Object"
  rnaObject=DESeqDataSetFromMatrix(countData = unnormalized_sampleData, colData = meta_condition,design = ~ uniqueCondition)
  levelList<-as.character(unique(meta_condition)[["uniqueCondition"]]);
  colData(rnaObject)$uniqueCondition=factor(colData(rnaObject)$uniqueCondition, levels=levelList)
  
  # Import Size Factors From "p6 Data Object"
  sizeFactors(rnaObject) <- sizeFactors_p6
  
  return(rnaObject)
}

# Generate DeSeq Object (Without P1) function
generate_object=function(unnormalized_sampleData, meta_condition){
  
  # Generate "regular Data Object"
  rnaObject=DESeqDataSetFromMatrix(countData = unnormalized_sampleData, colData = meta_condition,design = ~ uniqueCondition)
  levelList<-as.character(unique(meta_condition)[["uniqueCondition"]]);
  colData(rnaObject)$uniqueCondition=factor(colData(rnaObject)$uniqueCondition, levels=levelList)
  
  # Calculate Size Factors
  rnaObject=estimateSizeFactors(rnaObject)
  
  return(rnaObject)
}

# Generate vsd varianceStabilizingTransformation function
generate_vsd=function(rnaObject){
  
  #Calculate Variance Stabilizing Transformations
  Vst=varianceStabilizingTransformation(rnaObject)
  vst=assay(Vst)
  
  # Return Result
  return(vst)
}

# Generate vsd varianceStabilizingTransformation function
generate_log2=function(rnaObject){
  
  #Normalize "regular Data Object" by using DeSeq size factors
  normalized_rna_data=counts(rnaObject, normalized=TRUE)
  
  # Calculate log2 transformation
  log2Norm=log2(normalized_rna_data+1)
  
  # Return Result
  return(log2Norm)
}


generateNormalizedData=function(unnormalized_rna_Input, 
                                unnormalized_mrna_Input,
                                unnormalized_protein_Input,
                                unnormalized_protein_Input_wo_NA,
                                unnormalized_protein_Input_wo_NAx6,
                                #meta_Input,
                                conditionNumberChoice, # "uniqueCondition", "uniqueCondition02"
                                badDataFilterSetChoice, # "set01" , "set02"
                                dataTypeChoice,  # can be "rna", "mrna","protein","protein_wo_NA","protein_wo_NAx6"
                                dataTimeChoice, # can be "wholeSet", "exponential", "stationary", "late_stationary"
                                MgLevelChoice,  # can be "lowMg", "midMg", "highMg", "allMg"
                                NaLevelChoice,  # can be "lowNa","highNa", "allNa"
                                carbonTypeChoice, # a letter combination from the list "SYAN"  S (glucose), Y (glycerol), A (lactate), N (gluconate) 
                                filterTypeChoice,  # can be "noFilter", "mean", "sd", "max", "threshold"
                                deSeqNormChoice, # can be "p1", "p6", reg"
                                normalizationMethodChoice, # can be "vst", "log2"
                                experimentChoice) # can be "allEx", 
                                                  # or a combination of below
                                                  # "Stc" for "glucose_time_course", 
                                                  # "Ytc" for "glycerol_time_course", 
                                                  # "Nas" for "NaCl_stress", 
                                                  # "Agr" for "lactate_growth", 
                                                  # "Ngr" for "gluconate_growth", 
                                                  # "Mgl" for "MgSO4_stress_low", 
                                                  # "Mgh" for "MgSO4_stress_high"
                                {
  ###*****************************
  

  ###*****************************
  # Generating Meta Input
  meta_rna=meta_rna
  meta_protein=meta_protein
  if(grepl(pattern = "rna",x = dataTypeChoice)){meta_Input=meta_rna}
  if(grepl(pattern = "protein",x = dataTypeChoice)){meta_Input=meta_protein}
  ###*****************************
  
  ###*****************************
  # Order Experiment List
  experimentListOrder=c("allEx","Stc","Ytc","Nas","Agr","Ngr","Mgl","Mgh")
  currentOrder=experimentChoice
  experimentChoice=currentOrder[order(match(currentOrder,experimentListOrder))]
  #*******************************
  

  ###*****************************
  # FILE NAME GENERATE
  
  experimentList=paste(experimentChoice, collapse = '')
  longNameList=paste0(dataTypeChoice,"_",dataTimeChoice,
                      "_",MgLevelChoice,"_",NaLevelChoice,
                      "_",carbonTypeChoice,"_",badDataFilterSetChoice,
                      "_",experimentList)
  step01=paste0("unnormalized_",longNameList)
  step02=paste0("unnormalized_",filterTypeChoice,"_",longNameList)
  step03=paste0("normalized_",normalizationMethodChoice,"_",deSeqNormChoice,"_",
                filterTypeChoice,"_",longNameList)  
  conditionName=paste0("condition_",normalizationMethodChoice,"_",deSeqNormChoice,"_",
                       filterTypeChoice,"_",longNameList) 
  parametersName=paste0("parameters_",normalizationMethodChoice,"_",deSeqNormChoice,"_",
                       filterTypeChoice,"_",longNameList) 
  
  print(step03)
  ###*****************************
  
  
  
  ###*****************************
  # SELECT DATA // rna - mrna - protein - protein_wo_NA
  if(dataTypeChoice=="rna"){mainDataFrame=unnormalized_rna_Input}
  if(dataTypeChoice=="mrna"){mainDataFrame=unnormalized_mrna_Input}
  if(dataTypeChoice=="protein"){mainDataFrame=unnormalized_protein_Input}
  if(dataTypeChoice=="protein_wo_NA"){mainDataFrame=unnormalized_protein_Input_wo_NA}
  if(dataTypeChoice=="protein_wo_NAx6"){mainDataFrame=unnormalized_protein_Input_wo_NAx6}
  ###*****************************
  
  
  
  ###*****************************
  # SELECT GOOD SAMPLES FILTER OUT BAD ONES
  badDataFilterOutDataSet_set00=c()
  badDataFilterOutDataSet_set01=c("MURI_029","MURI_067","MURI_075","MURI_084",
                                  "MURI_086","MURI_091","MURI_136","MURI_138")
  badDataFilterOutDataSet_set02=c("MURI_029","MURI_067","MURI_075","MURI_084",
                                  "MURI_136","MURI_138")
  
  badDataFilterOutDataSets=get(grep(badDataFilterSetChoice,ls(), value = TRUE))
  badDataRemainingDataSets=base::setdiff(as.vector(meta_Input$dataSet),badDataFilterOutDataSets)
  

  mainDataFrame=mainDataFrame[,badDataRemainingDataSets]
  meta_Input %>%
    dplyr::filter(dataSet %in% badDataRemainingDataSets)->meta_Input
  ###*****************************
  
  
  
  ###*****************************
  # select dataSets
  
  if(dataTimeChoice!="wholeSet"){dataTimeVector=dataTimeChoice} 
  if(dataTimeChoice=="wholeSet"){dataTimeVector=as.vector(unique(meta_Input$growthPhase))}
  
  if(MgLevelChoice!="allMg"){MgLevelVector=MgLevelChoice} 
  if(MgLevelChoice=="allMg"){MgLevelVector=as.vector(unique(meta_Input$Mg_mM_Levels))}
  
  if(NaLevelChoice!="allNa"){NaLevelVector=NaLevelChoice} 
  if(NaLevelChoice=="allNa"){NaLevelVector=as.vector(unique(meta_Input$Na_mM_Levels))}
  
  carbonSourceVector = c("glucose", "glycerol", "lactate", "gluconate")
  carbonSourceConditionVector = c(grepl("S",carbonTypeChoice),grepl("Y",carbonTypeChoice),
                                  grepl("A",carbonTypeChoice),grepl("N",carbonTypeChoice))
  carbonSourceVector = carbonSourceVector[carbonSourceConditionVector]
  
  experimentVector = c("glucose_time_course", "glycerol_time_course", "NaCl_stress", "lactate_growth", 
                       "gluconate_growth", "MgSO4_stress_low", "MgSO4_stress_high")
  if(experimentChoice=="allEx"){experimentChoice=c("Stc","Ytc","Nas","Agr","Ngr","Mgl","Mgh")}
  experimentConditionVector = c(any(experimentChoice=="Stc"),
                                any(experimentChoice=="Ytc"),
                                any(experimentChoice=="Nas"),
                                any(experimentChoice=="Agr"),
                                any(experimentChoice=="Ngr"),
                                any(experimentChoice=="Mgl"),
                                any(experimentChoice=="Mgh"))
  experimentVector = experimentVector[experimentConditionVector]
  
  meta_Input %>%
    dplyr::filter(growthPhase %in% dataTimeVector) %>%
    dplyr::filter(Mg_mM_Levels %in% MgLevelVector) %>%
    dplyr::filter(Na_mM_Levels %in% NaLevelVector) %>%
    dplyr::filter(carbonSource %in% carbonSourceVector) %>%
    dplyr::filter(experiment %in% experimentVector)->meta_sub
  ###*****************************
  
  
  ###*****************************
  # !!! This part is only for GO analyze !!!
  if (all(experimentChoice==c("Stc","Mgl","Mgh"))){
  meta_sub %>% 
    dplyr::filter((experiment!="glucose_time_course") | 
                    (experiment=="glucose_time_course" & growthTime_hr %in% c(5,6,24)))->meta_sub
  }
  ###*****************************
  
  
  ###*****************************
  meta_sub %>%
    select_(conditionNumberChoice)->conditionDf
  row.names(conditionDf)<-meta_sub[["dataSet"]]
  conditionDf <- droplevels(conditionDf)
  colnames(conditionDf)<-"uniqueCondition"
  
  mainDataFrame=mainDataFrame[,as.vector(meta_sub$dataSet)]
  assign(value = mainDataFrame, x = step01)
  ###*****************************
  
  
  
  ###*****************************
  # FILTERING OUT GENES
  # I am going to try several methods for filtering 
  # a) No filtering at all
  # b) Filtering with mean >1
  # c) Filtering with sd >10
  # d) Filtering with max>10
  
  if(filterTypeChoice=="noFilter"){mainDataFrame = mainDataFrame }
  if(filterTypeChoice=="mean"){mainDataFrame = mainDataFrame[(apply(mainDataFrame,1,mean)>1),]}
  if(filterTypeChoice=="sd"){mainDataFrame = mainDataFrame[(apply(mainDataFrame,1,sd)>10),]}
  if(filterTypeChoice=="max"){mainDataFrame = mainDataFrame[(apply(mainDataFrame,1,max)>10),]}
  if(filterTypeChoice=="threshold"){mainDataFrame[mainDataFrame<3]<-0}
  assign(value = mainDataFrame, x = step02)
  ###*****************************
  
  
  
  ###*****************************
  # GENERATE P1 or Non P1 DeSeq OBJECT
  if (deSeqNormChoice=="p1")
    {deseqObject=generate_p1_object(mainDataFrame, conditionDf)}
  if (deSeqNormChoice=="p6")
    {deseqObject=generate_p6_object(mainDataFrame, conditionDf)
    sizeFactors(deseqObject)<-sizeFactors(deseqObject)*6}
  if (deSeqNormChoice=="reg")
    {deseqObject=generate_object(mainDataFrame, conditionDf)}
  
  if(normalizationMethodChoice=="vst")
    {mainDataFrame=generate_vsd(deseqObject)}
  if(normalizationMethodChoice=="log2")
  {mainDataFrame=generate_log2(deseqObject)}
  if(normalizationMethodChoice=="none")
  {mainDataFrame=counts(deseqObject, normalized=TRUE)}
  
  
  assign(value = mainDataFrame, x = step03)
  assign(value = meta_sub, x = conditionName)
  ###*****************************
  
  
  ###*****************************
  # GENERATE PARAMETER NAMES
  parameters<-list()
  
  parameters$conditionNumber=conditionNumberChoice
  parameters$badDataFilterSet=badDataFilterSetChoice
  parameters$dataType=dataTypeChoice
  parameters$dataTime=dataTimeChoice
  parameters$MgLevel=MgLevelChoice
  parameters$NaLevel=NaLevelChoice
  parameters$carbonType=carbonTypeChoice
  parameters$filterType=filterTypeChoice
  parameters$deSeqNorm=deSeqNormChoice
  parameters$normalizationMethod=normalizationMethodChoice
  parameters$experimet=experimentChoice
  
  assign(x=parametersName, value = parameters)
  ###*****************************
  
  
  ###*****************************
  #SAVE FILES
  savedFilename=paste0("../a_results/",step03,".RData")
  save(list = c(step02,step03,conditionName, parametersName),
       file=savedFilename, compress = "xz")
  
  return(mainDataFrame)
}