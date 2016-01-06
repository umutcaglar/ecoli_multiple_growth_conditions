# Go functions // Calculation of slope and inverse p

spearman_function<-function(main_data,x_axis_input,y_axis_input){
  
  
  var_inp=x_axis_input
  var_out="x_axis"
  mutateFunction <-function(name){u=name; return(u)}
  dots <- list(interp(~mutateFunction(var), var = as.name(var_inp)))
  main_data %>%  dplyr::mutate_(.dots = setNames(dots, c(var_out))) -> main_data

    
  var_inp=y_axis_input
  var_out="y_axis"
  mutateFunction <-function(name){u=name; return(u)}
  dots <- list(interp(~mutateFunction(var), var = as.name(var_inp)))
  main_data %>%  dplyr::mutate_(.dots = setNames(dots, c(var_out))) -> main_data
  
  main_data %>% 
    dplyr::group_by(dataSet) %>% 
    dplyr::summarize(x_axis=unique(x_axis)) %>% 
    dplyr::group_by(x_axis) %>% 
    dplyr::summarise(numSamples=length(x_axis))->repeatTable
  
  for (counter01 in 1:nrow(repeatTable)){
    setSize=repeatTable$numSamples[counter01]
    b=rep(counter01,setSize)
    
    if(counter01==1){
      bb=b}
    if(counter01!=1){
      bb=c(bb,b)}
  }
  
  aa= seq(1,sum(repeatTable$numSamples))
  
  main_data %>%
    dplyr::group_by(ID) %>%
    dplyr::arrange(y_axis) %>%
    dplyr::summarize(p_value_wrep=cor.test(x_axis,y_axis, 
                                           method = "spearman", 
                                           alternative = "two.sided")[["p.value"]],
                     p_value_worep=cor.test(aa,x_axis, 
                                            method = "spearman", 
                                            alternative = "two.sided")[["p.value"]],
                     rho_values=cor.test(x_axis,y_axis, 
                                              method = "spearman", 
                                              alternative = "two.sided")[["estimate"]],
                     p_max_wrep=cor.test(bb,y_axis, 
                                         method = "spearman", 
                                         alternative = "two.sided")[["p.value"]],
                     p_max_worep=cor.test(aa,bb, 
                                          method = "spearman", 
                                          alternative = "two.sided")[["p.value"]],
                     p_ratio_wrep=-log10(p_value_wrep)/(-log10(p_max_wrep)),
                     p_ratio_worep=-log10(p_value_worep)/(-log10(p_max_worep)),
                     gene_name=unique(gene_name)) %>%
    dplyr::mutate(p_related_wrep=-log10(p_value_wrep)) %>%
    dplyr::mutate(p.adj_wrep=p.adjust(p=p_value_wrep, method = "fdr")) %>%
    dplyr::mutate(p.adj_related_wrep = -log10(p.adj_wrep)) %>%
    dplyr::mutate(p_related_sign=sign(rho_values)) %>%
    dplyr::arrange(desc(p_related_wrep))->summaryDataFrame
  
  return(summaryDataFrame)
}
###************************************


###************************************
ICC_function<-function(main_data,x_axis_input,y_axis_input){
  
  
  var_inp=x_axis_input
  var_out="x_axis"
  mutateFunction <-function(name){u=name; return(u)}
  dots <- list(interp(~mutateFunction(var), var = as.name(var_inp)))
  main_data %>%  dplyr::mutate_(.dots = setNames(dots, c(var_out))) -> main_data
  
  
  var_inp=y_axis_input
  var_out="y_axis"
  mutateFunction <-function(name){u=name; return(u)}
  dots <- list(interp(~mutateFunction(var), var = as.name(var_inp)))
  main_data %>%  dplyr::mutate_(.dots = setNames(dots, c(var_out))) -> main_data
  
  
  fitF<-function(dataInput,x_axiss,y_axiss){
    result=list() 
 
    vec=as.vector(unique(dataInput[[x_axiss]]))
    q=gtools::permutations(n = length(vec), 
                           r = length(vec), 
                           v = paste0("x",c(1:length(vec))))
    
    
    for(counter01 in 1:nrow(q)){
      dataInput %>%
        dplyr::mutate(x_axis2=x_axis)->dataInput
      
      for(counter02 in 1: length(vec)){
        vec2=as.vector(dataInput$x_axis2)
        
        oldName=vec[counter02]
        newName=paste0(q[counter01,counter02],vec[counter02])
        
        vec2=replace(vec2,vec2==oldName,newName)
        dataInput$x_axis2<-vec2
      }
      
      dataInput$x_axis2<-factor(dataInput$x_axis2)

      ICC_val=abs(ICCbare(x_axis2, get(y_axiss), dataInput))
      if(counter01==1){maxIccVal=ICC_val}
      if(counter01!=1){maxIccVal=max(ICC_val,maxIccVal)}
    }
    
    result$ICC_val=maxIccVal
    
    
    return(result)
  }
  
  
  main_data %>%
    dplyr::group_by(ID,gene_name) %>%
    dplyr::do(fit_result=fitF(.,"x_axis","y_axis"))->summaryDataFrame
  
  summaryDataFrame %>% 
    group_by(ID,gene_name) %>%
    do(., do.call(data.frame, .$fit_result))->summaryDataFrame

  
  return(summaryDataFrame)
}
###************************************


###************************************
spearman_function_v2<-function(main_data,x_axis_input,y_axis_input){
  
  # Generate a column named x-axis
  var_inp=x_axis_input
  var_out="x_axis"
  mutateFunction <-function(name){u=name; return(u)}
  dots <- list(interp(~mutateFunction(var), var = as.name(var_inp)))
  main_data %>%  dplyr::mutate_(.dots = setNames(dots, c(var_out))) -> main_data
  
  # generate a column named y-axis
  var_inp=y_axis_input
  var_out="y_axis"
  mutateFunction <-function(name){u=name; return(u)}
  dots <- list(interp(~mutateFunction(var), var = as.name(var_inp)))
  main_data %>%  dplyr::mutate_(.dots = setNames(dots, c(var_out))) -> main_data
  
  main_data %>% 
    dplyr::group_by(dataSet) %>% 
    dplyr::summarize(x_axis=unique(x_axis)) %>% 
    dplyr::group_by(x_axis) %>% 
    dplyr::summarise(numSamples=length(x_axis))->repeatTable
  
  for (counter01 in 1:nrow(repeatTable)){
    setSize=repeatTable$numSamples[counter01]
    b=rep(counter01,setSize)
    
    if(counter01==1){
      bb=b}
    if(counter01!=1){
      bb=c(bb,b)}
  }
  
  aa= seq(1,sum(repeatTable$numSamples))
  
  main_data %>%
    dplyr::group_by(ID) %>%
    dplyr::arrange(y_axis) %>%
    dplyr::summarize(p_value=cor.test(x_axis,y_axis, 
                                      method = "spearman", 
                                      alternative = "two.sided")[["p.value"]],
                     rho_values=cor.test(x_axis,y_axis, 
                                         method = "spearman", 
                                         alternative = "two.sided")[["estimate"]],
                     gene_name=unique(gene_name)) %>%
    dplyr::mutate(p_related=-log10(p_value)) %>%
    dplyr::mutate(p.adj=p.adjust(p=p_value, method = "fdr")) %>%
    dplyr::mutate(p.adj_related= -log10(p.adj)) %>%
    dplyr::mutate(p_related_sign=sign(rho_values)) %>%
    dplyr::arrange(desc(p_related))->summaryDataFrame
  
  
  randomResultsForz<-function(dataInput){
    #print(dataInput$x_axis)
    
    scoreVector=rep(0,20)
    for(counter01 in 1:20){
      fake_x_axis=sample(dataInput$x_axis)
      fake_y_axis=dataInput$y_axis
      fake_testResult=cor.test(fake_x_axis,fake_y_axis, 
                               method = "spearman", 
                               alternative = "two.sided")
      scoreVector[counter01]=sign(fake_testResult[["estimate"]])*(-log10(fake_testResult[["p.value"]]))
    }
    outPut=list(meanVal=mean(scoreVector),stdVal=sd(scoreVector))
    return(outPut)
  }
  
  main_data %>%
    dplyr::group_by(ID) %>% 
    do(ststOut=randomResultsForz(.))->statisticalTable
  
  statisticalTable<-do(statisticalTable, do.call(data.frame,list(ID=.$ID,
                                                                  meanVal=.$ststOut$meanVal,
                                                                  stdVal=.$ststOut$stdVal)))
  summaryDataFrame<-dplyr::left_join(summaryDataFrame,statisticalTable)
  summaryDataFrame %>%
    dplyr::mutate(z_score=((p_related*p_related_sign)-meanVal)/stdVal)->summaryDataFrame
  
  
  
  return(summaryDataFrame)
}