calculateClusterMeans<-function(variableInput,metaInput,mainDataFrame)
{
  metaInput %>%
    dplyr::select_("dataSet","variable"=variableInput)->meta_variable
  
  
  mainDataFrame %>%
    dplyr::add_rownames("rowNames")%>%
    tidyr::gather(key = colNames,
                  value = cophenetic_distance, 
                  ...=1:ncol(mainDataFrame)+1) ->mainDataFrame_tidy
  

  mainDataFrame_tidy %>%
    mutate(cophenetic_distanceNa=ifelse(colNames==rowNames,NA,cophenetic_distance))%>%
    left_join(.,meta_variable, by=c("rowNames"="dataSet"))%>%
    dplyr::rename(variableR=variable)%>%
    left_join(.,meta_variable, by=c("colNames"="dataSet"))%>%
    dplyr::rename(variableC=variable) %>%
    dplyr::mutate(cophenetic_distanceNa=ifelse(variableC!=variableR,
                                               NA,cophenetic_distanceNa))->mainDataFrame_tidy
  

  mainDataFrame_tidy %>%
    dplyr::group_by(variableC)%>%
    dplyr::summarise(meanValue=mean(cophenetic_distanceNa,na.rm = TRUE))->mainDataFrame_summary
  mainDataFrame_summary$variableC<-as.character(mainDataFrame_summary$variableC)
  
  mainDataFrame_summary<-rbind(mainDataFrame_summary,
                               c(variableC="allData",
                                 meanValue=mean(mainDataFrame_tidy$cophenetic_distanceNa,na.rm = TRUE)))
  mainDataFrame_summary$meanValue=as.numeric(mainDataFrame_summary$meanValue)
  
  row.names(mainDataFrame_summary)=mainDataFrame_summary[["variableC"]]
  mainDataFrame_summary[["variableC"]]<-NULL
  mainDataFrame_summary<-as.data.frame(t(mainDataFrame_summary))
  
  return(mainDataFrame_summary)
}


reShapeZScores<-function(variableInput,metaInput,z_data)
{
  old_z<-z_data
  z_data<-as.data.frame(t(z_data))
  z_data%>%
    dplyr::add_rownames(var="Condition")->z_data
  colnames(z_data)<-c("Condition","Z_score")
  z_data[["Z_score"]]<-as.numeric(z_data[["Z_score"]])
  
  Overall_Z_score<-z_data[nrow(z_data),]
  z_data<-z_data[-c(nrow(z_data)), ] 
  
  metaInput %>%
    dplyr::group_by_(variableInput)%>%
    dplyr::summarise(num_elements=length(dataSet))->metaInput_summary
  metaInput_summary[[variableInput]]<-as.character(metaInput_summary[[variableInput]])
  
  
  z_data2<-cbind(Variable=variableInput,
                 Overall_Z_score=Overall_Z_score[[2]],
                 z_data)
  z_data3<-dplyr::left_join(z_data2,metaInput_summary,by=c("Condition"=variableInput))
  return(z_data3)
}