reShapeZScores<-function(variableInput,metaInput,z_data)
{
  old_z<-z_data
  z_data<-as.data.frame(t(z_data))
  z_data%>%
    dplyr::add_rownames(var="Condition")->z_data
  colnames(z_data)<-c("Condition","Z_score")
  z_data[["Z_score"]]<-as.numeric(z_data[["Z_score"]])
  
  Overall_Z_score<-data[nrow(z_data),]
  z_data<-z_data[-c(nrow(z_data)), ] 
  
  metaInput %>%
    dplyr::group_by_(variableInput)%>%
    dplyr::summarise(num_elements=length(dataSet))->metaInput_summary
  metaInput_summary[[variableInput]]<-as.character(metaInput_summary[[variableInput]])
  
  
  z_data<-cbind(variableInput,Overall_Z_score,z_data)
  z_data<-dplyr::left_join(z_data,metaInput_summary,by=c("Condition"=variableInput))
  return(z_data)
}



z_batch_re<-reShapeZScores(variableInput="batchNumber",
                           metaInput=conditionSummary,
                           z_data=z_batch)




data<-t(data)
colnames(data)<-"Z-score"

Overall_Z_score<-data[nrow(data)]