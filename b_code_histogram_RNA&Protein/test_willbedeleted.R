mainDataFrame %>%
  dplyr::add_rownames("rowNames")%>%
  tidyr::gather(key = colNames,
                value = cophenetic_distance, 
                ...=1:ncol(mainDataFrame)+1) ->q1

q1 %>%
  mutate(cophenetic_distanceNa=ifelse(colNames==rowNames,NA,cophenetic_distance))%>%
  left_join(.,meta_variable, by=c("rowNames"="dataSet"))%>%
  dplyr::rename(variableR=variable)%>%
  left_join(.,meta_variable, by=c("colNames"="dataSet"))%>%
  dplyr::rename(variableC=variable) %>%
  dplyr::mutate(cophenetic_distanceNa=ifelse(variableC!=variableR,
                                             NA,cophenetic_distanceNa))->q2


q2%>%
  dplyr::select(rowNames, variableR, colNames,cophenetic_distanceNa)%>%
  tidyr::spread(key = colNames,
                value = cophenetic_distanceNa)%>%
  dplyr::arrange(variableR)->q3

m<-c("rowNames","variableR",as.vector(q3[["rowNames"]]))
q3[,m]->q4











distance_df1=real_distance_df
distance_df2=cophenetic_distance_df_real
