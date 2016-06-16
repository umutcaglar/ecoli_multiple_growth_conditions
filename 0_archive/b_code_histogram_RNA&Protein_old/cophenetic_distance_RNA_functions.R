###*****************************
# Distance data frame calculation
distanceDataFrame<-function(dataMatrix, metaData, condition){
  a=condition
  metaData %>% dplyr::filter(condition == a) -> subMeta
  categoryVar=unique(subMeta$category)
  subMeta=as.vector(subMeta$dataSet)->subMeta
  dataMatrix[subMeta,subMeta]->subDataMatrix
  subDataMatrix[lower.tri(subDataMatrix,diag = T)]<-NA
  numberOfSamples=nrow(subDataMatrix)
  
  subDataMatrix=as.data.frame(subDataMatrix)
  colNameList=colnames(subDataMatrix)
  subDataMatrix %>% 
    dplyr::mutate(rowNames=row.names(subDataMatrix))%>%
    tidyr::gather(key = colNames, value = cophenetic_distance, seq(1:length(colNameList)))%>%
    dplyr::filter(!is.na(cophenetic_distance)) %>%
    dplyr::mutate(condition=condition,
                  category=categoryVar)->subDataMatrix
  
  output<-list()
  output[[1]]=subDataMatrix
  output[[2]]=numberOfSamples
  return(output)
}
###*****************************


###*****************************
# mean distance for a spesific category
mainDistanceCategory<-function(dataMatrix, metaData, category){
  a=category
  metaData %>% dplyr::filter(category == a) -> subMeta
  conditions=unique(subMeta$condition)
  
  for(counter01 in 1: length(conditions)){
    condition=conditions[counter01]
    temp=distanceDataFrame(dataMatrix, metaData, condition)
    temp=temp[[1]]
    if(counter01==1){combinedDistance=temp}
    if(counter01!=1){combinedDistance=rbind(combinedDistance,temp)}
  }
  
  combinedDistance %>%
    dplyr::group_by(condition) %>%
    dplyr::summarize(lengthVar=length(cophenetic_distance),
                     meanVal=mean(cophenetic_distance),
                     medianVal=median(cophenetic_distance),
                     numVar=(0.5)*(1+sqrt(1+8*lengthVar)),
                     category=unique(category)) %>%
    dplyr::mutate(weighted_Mean=sum(meanVal*numVar)/sum(numVar),
                  weighted_Median=sum(medianVal*numVar)/sum(numVar),
                  overal_Mean=sum(meanVal*lengthVar)/sum(lengthVar))->combinedDistance
  return(combinedDistance)
}
###*****************************


###*****************************
# fake mean distance for a spesific category
mainDistanceCategoryF<-function(dataMatrix, metaData, category){
  a=category
  metaData %>% dplyr::filter(category == a) -> subMeta
  
  subMetaF=subMeta
  subMetaF %>% 
    dplyr::mutate(condition=sample(condition)) -> subMetaF
  
  fakeMean=mainDistanceCategory(dataMatrix, subMetaF, category)
  
  return(fakeMean)
}