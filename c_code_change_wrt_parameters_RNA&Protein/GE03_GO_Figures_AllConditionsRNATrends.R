# Best seperated genes in Glucose/ Glycerol / Lactate / Gluconate


# trend-GO analysis of MgSO4 

###*****************************
# INITIAL COMMANDS TO RESET THE SYSTEM

rm(list = ls())
if (is.integer(dev.list())){dev.off()}
cat("\014")

# Set Working Directory
setwd('/Users/umut/GitHub/AG3C_Analyze/initialPaper03c_change_wrt_parameters_RNA&Protein/') # mac computer
###*****************************




###*****************************
# INSTALL LIBRARIES

library("dplyr")
library("tidyr")
library("lazyeval")
require("flashClust")

# For Plotting
library("ggplot2")
library("RColorBrewer")
library("gridExtra")
library("cowplot")
require("ggdendro")
require("scales")
library("VennDiagram")
require("gtable")
###*****************************


#*********************************
# Not In Operator 
"%!in%" <- function(x,table) match(x,table, nomatch = 0) == 0
#*********************************


#*********************************
# ordering Function
orderFunction<-function(df,column,target){
  
  temp=df[[column]]
  temp2<-factor(temp, levels = target)
  df[column]<-temp2
  
  df %>% 
    dplyr::select_(column) %>%
    dplyr::mutate(order=seq(1:nrow(df)))%>%
    dplyr::group_by_(column) %>%
    dplyr::summarize(order_length=length(order))%>%
    dplyr::mutate(setNo=seq(1:length(target)))->dfs
  
  correctOrderDf<-dplyr::left_join(df,dfs)
  correctOrderDf %>% dplyr::arrange(setNo)->correctOrderDf
  
  correctOrderDf %>% dplyr::select(-c(setNo,order_length))->df
  df %>% dplyr::group_by()->df
  return(df)
}
#*********************************


###*****************************
# PARAMETERS FOR DATA
conditionNumberChoice="uniqueCondition"
dataTypeChoice="mrna" # can be "rna", "mrna","protein","protein_wo_NA","protein_wo_NAx6"
badDataFilterSetChoice="set02" # "set00", "set01", "set02"
dataTimeChoice="wholeSet" # exponential / stationary / late_stationary / wholeSet
MgLevelChoice="allMg" # allMg highMg midMg lowMg
NaLevelChoice="allNa" # allNa lowNa highNa
carbonTypeChoice="SYAN" # a letter combination from the list "SYAN"  
#S (glucose), Y (glycerol), A (lactate), N (gluconate)
filterTypeChoice="mean" # mean/sd/ max/ noFilter
deSeqNormChoice="p1"
normalizationMethodChoice="vst" # "vst" , "log2" 
experimentChoice=c("allEx") # can be "allEx", 
# or a combination of below
# "Stc" for "glucose_time_course", 
# "Ytc" for "glycerol_time_course", 
# "Nas" for "NaCl_stress", 
# "Agr" for "lactate_growth", 
# "Ngr" for "gluconate_growth", 
# "Mgh" for "MgSO4_stress_low", 
# "Mgl" for "MgSO4_stress_high"

test_type="spearman" # "spearman", "ICC"

dictionaryPick="BL"
# there are 4 alternatves for dictionary
# barricks own dictionary -Barrick Lab (BL)
# what I found on web (WEB)
# a combination of 2 -Combined (COM)
# Protein dictionary (PD)
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
step02=paste0("unnormalized_",filterTypeChoice,"_",longNameList)
step03=paste0("goDetailedInput_",normalizationMethodChoice,"_",deSeqNormChoice,"_",
              filterTypeChoice,"_",longNameList,"_",test_type) 
step04=paste0(normalizationMethodChoice,"_",deSeqNormChoice,"_",
              filterTypeChoice,"_",longNameList,"_",test_type) 
generalFigureName=paste0(step04,"_dic",dictionaryPick)
###*****************************


###*****************************
# Load data
name=paste0(step03,"_test_glucose_dic",dictionaryPick)
test_glucose_result=read.csv(file = paste0("../initialPaper03r/",name,".csv"))
test_glucose_result %>% dplyr::mutate(test_condition="glucose",
                                      overal_test_condition="carbonSource")->test_glucose_result

name=paste0(step03,"_test_glycerol_dic",dictionaryPick)
test_glycerol_result=read.csv(file = paste0("../initialPaper03r/",name,".csv"))
test_glycerol_result %>% dplyr::mutate(test_condition="glycerol",
                                       overal_test_condition="carbonSource")->test_glycerol_result

name=paste0(step03,"_test_lactate_dic",dictionaryPick)
test_lactate_result=read.csv(file = paste0("../initialPaper03r/",name,".csv"))
test_lactate_result %>% dplyr::mutate(test_condition="lactate",
                                      overal_test_condition="carbonSource")->test_lactate_result

name=paste0(step03,"_test_gluconate_dic",dictionaryPick)
test_gluconate_result=read.csv(file = paste0("../initialPaper03r/",name,".csv"))
test_gluconate_result %>% 
  dplyr::mutate(test_condition="gluconate",
                overal_test_condition="carbonSource")->test_gluconate_result

name=paste0(step03,"_test_Na_mM_dic",dictionaryPick)
test_Na_mM_result=read.csv(file = paste0("../initialPaper03r/",name,".csv"))
test_Na_mM_result %>% 
  dplyr::mutate(test_condition="Na_mM",
                overal_test_condition="Na")->test_Na_mM_result

name=paste0(step03,"_test_Mg_mM_dic",dictionaryPick)
test_Mg_mM_result=read.csv(file = paste0("../initialPaper03r/",name,".csv"))
test_Mg_mM_result %>% 
  dplyr::mutate(test_condition="Mg_mM",
                overal_test_condition="Mg")->test_Mg_mM_result

name=paste0(step03,"_test_growthTime_hr_dic",dictionaryPick)
test_growthTime_hr_result=read.csv(file = paste0("../initialPaper03r/",name,".csv"))
test_growthTime_hr_result %>% 
  dplyr::mutate(test_condition="growthTime_hr",
                overal_test_condition="time")->test_growthTime_hr_result

name=paste0(step03,"_test_exponential_dic",dictionaryPick)
test_exponential_result=read.csv(file = paste0("../initialPaper03r/",name,".csv"))
test_exponential_result %>% 
  dplyr::mutate(test_condition="exponential",
                overal_test_condition="time")->test_exponential_result

name=paste0(step03,"_test_stationary_dic",dictionaryPick)
test_stationary_result=read.csv(file = paste0("../initialPaper03r/",name,".csv"))
test_stationary_result %>% 
  dplyr::mutate(test_condition="stationary",
                overal_test_condition="time")->test_stationary_result

name=paste0(step03,"_test_late_stationary_dic",dictionaryPick)
test_late_stationary_result=read.csv(file = paste0("../initialPaper03r/",name,".csv"))
test_late_stationary_result %>% 
  dplyr::mutate(test_condition="late_stationary",
                overal_test_condition="time")->test_late_stationary_result
###*****************************


###*****************************
# Combining results
test_combined_result=dplyr::rbind_all(list(test_glucose_result,
                                           test_glycerol_result,
                                           test_gluconate_result,
                                           test_lactate_result,
                                           test_Na_mM_result,
                                           test_Mg_mM_result,
                                           test_growthTime_hr_result,
                                           test_exponential_result,
                                           test_stationary_result,
                                           test_late_stationary_result))
###*****************************


###*****************************
# add new variable "carbon source"
test_combined_result %>%
  dplyr::group_by() %>%
  dplyr::filter(test_condition %in% c("glucose","gluconate","lactate","glycerol"))%>%
  dplyr::select(ID, gene_name, p_related, test_condition)%>%
  tidyr::spread(test_condition, p_related)%>%
  dplyr::mutate(p_related=glucose+gluconate+lactate+glycerol,
                test_condition="carbonSource",
                overal_test_condition="carbonSource",
                p_value=10^-(p_related))%>%
  dplyr::select(-glucose,-glycerol,-lactate,-gluconate) %>%
  dplyr::arrange(desc(p_related))-> test_carbon_source_result

test_combined_result %>%
  dplyr::group_by() %>%
  dplyr::filter(test_condition %in% c("glucose","gluconate","lactate","glycerol"))%>%
  dplyr::select(ID, gene_name, p.adj_related, test_condition)%>%
  tidyr::spread(test_condition, p.adj_related)%>%
  dplyr::mutate(p.adj_related=glucose+gluconate+lactate+glycerol,
                test_condition="carbonSource",
                overal_test_condition="carbonSource",
                p.adj=10^-(p.adj_related))%>%
  dplyr::select(-glucose,-glycerol,-lactate,-gluconate) %>%
  dplyr::arrange(desc(p.adj_related))-> test_carbon_source_result2

test_carbon_source_result=dplyr::left_join(test_carbon_source_result,test_carbon_source_result2)
dplyr::rbind_all(list(test_combined_result,test_carbon_source_result))->test_combined_result
###*****************************


###*****************************
# Reordering Results
test_combined_result$test_condition=factor(test_combined_result$test_condition, 
                                           levels = c("glucose",
                                                      "glycerol",
                                                      "lactate",
                                                      "gluconate",
                                                      "carbonSource",
                                                      "Na_mM",
                                                      "Mg_mM",
                                                      "growthTime_hr",
                                                      "exponential",
                                                      "stationary",
                                                      "late_stationary"))
###*****************************


###*****************************
listColors=c("glucose"="#bcbddc",
             "glycerol"="#9e9ac8",
             "lactate"="#807dba",
             "gluconate"="#6a51a3",
             "carbonSource"="#4a1486", # carbon sources
             "Na_mM"="#fd8d3c", # Na
             "Mg_mM"="#6baed6", # Mg
             "growthTime_hr"="#005a32",
             "exponential"="#bae4b3",
             "stationary"="#74c476",
             "late_stationary"="#238b45") # time

# Draw figure 1 violin
fig01a<-ggplot(test_combined_result,aes(x=test_condition,y=p_related))+
  geom_violin(scale = "area",aes(fill=test_condition))+
  scale_fill_manual(values=listColors)+
  guides(fill = guide_legend(override.aes = list(colour = NULL)))+
  scale_y_continuous(expand = c(0,0))+
  ylab("-log10(p_value)")+
  xlab("test condition")+
  ggtitle("-log10(p value) distribution")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

fig01a

#***********************************


#***********************************
# generate the same figure with negative values and without carbon source
test_combined_result %>%
  dplyr::filter(test_condition %!in% c("carbonSource","exponential","stationary","late_stationary")) %>%
  dplyr::mutate(p_with_sign=p_related_sign*p.adj_related) ->test_combined_result2

listColors=c("glucose"="#bcbddc",
             "glycerol"="#9e9ac8",
             "lactate"="#807dba",
             "gluconate"="#6a51a3", # carbon sources
             "Na_mM"="#fd8d3c", # Na
             "Mg_mM"="#6baed6", # Mg
             "growthTime_hr"="#005a32") # time

fig01b<-ggplot(test_combined_result2,aes(x=test_condition,y=p_with_sign))+
  geom_violin(scale = "width",aes(fill=test_condition))+
  scale_fill_manual(values=listColors)+
  geom_hline(aes(yintercept=c(-log10(0.05))),colour="red",linetype="longdash")+
  geom_hline(aes(yintercept=c(log10(0.05))),colour="red",linetype="longdash")+
  geom_hline(aes(yintercept=c(-log10(0.01))),colour="orange",linetype="longdash")+
  geom_hline(aes(yintercept=c(log10(0.01))),colour="orange",linetype="longdash")+
  geom_hline(aes(yintercept=c(-log10(0.001))),colour="green",linetype="longdash")+
  geom_hline(aes(yintercept=c(log10(0.001))),colour="green",linetype="longdash")+
  #scale_y_log10()+
  guides(fill = guide_legend(override.aes = list(colour = NULL)))+
  ylab("Gene Score")+
  xlab("test condition")+
  #ggtitle("Gene Score Distribution")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

print(fig01b)
#***********************************



###*****************************


###*****************************
# parameters related with sets
test_combined_result %>%
  dplyr::filter(test_condition %in% c("carbonSource","Na_mM","Mg_mM","growthTime_hr"))%>%
  dplyr::filter(p.adj<0.05)%>%
  dplyr::select(ID, gene_name,test_condition)%>%
  dplyr::mutate(exist=1)%>%
  tidyr::spread(key = test_condition, value = exist, fill = 0)%>%
  dplyr::mutate(Na_mM=0)%>%
  dplyr::group_by(carbonSource,Mg_mM,Na_mM,growthTime_hr)%>%
  dplyr::summarise(numElements=length(carbonSource))->test_combined_resultS

test_combined_resultS %>% 
  dplyr::group_by()%>%
  dplyr::filter(carbonSource==1)%>%
  dplyr::summarize(sum(numElements))->carbonSource_data

test_combined_resultS %>% 
  dplyr::group_by()%>%
  dplyr::filter(Mg_mM==1)%>%
  dplyr::summarize(sum(numElements))->Mg_mM

test_combined_resultS %>% 
  dplyr::group_by()%>%
  dplyr::filter(Na_mM==1)%>%
  dplyr::summarize(sum(numElements))->Na_mM

test_combined_resultS %>% 
  dplyr::group_by()%>%
  dplyr::filter(growthTime_hr==1)%>%
  dplyr::summarize(sum(numElements))->growthTime_hr

test_combined_resultS %>% 
  dplyr::group_by()%>%
  dplyr::filter(carbonSource==1 & Mg_mM==1)%>%
  dplyr::summarize(sum(numElements))->carbon_Mg_data

test_combined_resultS %>% 
  dplyr::group_by()%>%
  dplyr::filter(carbonSource==1 & Na_mM==1)%>%
  dplyr::summarize(sum(numElements))->carbon_Na_data

test_combined_resultS %>% 
  dplyr::group_by()%>%
  dplyr::filter(carbonSource==1 & growthTime_hr==1)%>%
  dplyr::summarize(sum(numElements))->carbon_growth_data

test_combined_resultS %>% 
  dplyr::group_by()%>%
  dplyr::filter(Mg_mM==1 & Na_mM==1)%>%
  dplyr::summarize(sum(numElements))->Mg_Na_data

test_combined_resultS %>% 
  dplyr::group_by()%>%
  dplyr::filter(Mg_mM==1 & growthTime_hr==1)%>%
  dplyr::summarize(sum(numElements))->Mg_growth_data

test_combined_resultS %>% 
  dplyr::group_by()%>%
  dplyr::filter(Na_mM==1 & growthTime_hr==1)%>%
  dplyr::summarize(sum(numElements))->Na_growth_data

test_combined_resultS %>% 
  dplyr::group_by()%>%
  dplyr::filter(carbonSource==1 & Mg_mM==1 & Na_mM==1)%>%
  dplyr::summarize(sum(numElements))->carbon_Mg_Na_data

test_combined_resultS %>% 
  dplyr::group_by()%>%
  dplyr::filter(carbonSource==1 & Mg_mM==1 & growthTime_hr==1)%>%
  dplyr::summarize(sum(numElements))->carbon_Mg_growth_data

test_combined_resultS %>% 
  dplyr::group_by()%>%
  dplyr::filter(carbonSource==1 & Na_mM==1 & growthTime_hr==1)%>%
  dplyr::summarize(sum(numElements))->carbon_Na_growth_data

test_combined_resultS %>% 
  dplyr::group_by()%>%
  dplyr::filter(Mg_mM==1 & Na_mM==1 & growthTime_hr==1)%>%
  dplyr::summarize(sum(numElements))->Mg_Na_growth_data

test_combined_resultS %>% 
  dplyr::group_by()%>%
  dplyr::filter(carbonSource==1 & Mg_mM==1 & Na_mM==1 & growthTime_hr==1)%>%
  dplyr::summarize(sum(numElements))->carbon_Mg_Na_growth_data

png(filename = paste0("../initialPaper03f/venn_",generalFigureName,".png")
    ,width = 1600,height = 800, units = "px")

fig02<-VennDiagram::draw.quad.venn(area1 = carbonSource_data[[1]],
                                   area2 = Mg_mM[[1]],
                                   area3 = Na_mM[[1]],
                                   area4 = growthTime_hr[[1]],
                                   n12 = carbon_Mg_data[[1]],
                                   n13 = carbon_Na_data[[1]],
                                   n14 = carbon_growth_data[[1]],
                                   n23 = Mg_Na_data[[1]],
                                   n24 = Mg_growth_data[[1]],
                                   n34 = Na_growth_data[[1]],
                                   n123 = carbon_Mg_Na_data[[1]],
                                   n124 = carbon_Mg_growth_data[[1]],
                                   n134 = carbon_Na_growth_data[[1]],
                                   n234 = Mg_Na_growth_data[[1]],
                                   n1234 = carbon_Mg_Na_growth_data[[1]],
                                   category=c("Carbon source","Mg Strees","Na Stress","Growth phase"),
                                   fill=c("#decbe4","#b3cde3","#fdcdac","#ccebc5"),
                                   cex=rep(3,15),
                                   cat.cex =rep(3.5,4)
)
dev.off()

png(paste0("../initialPaper03f/vennTriplet_",generalFigureName,".png")
    ,width = 1000,height = 800, units = "px")

VennDiagram::draw.triple.venn(area1 = carbonSource_data[[1]],
                              area2 = Mg_mM[[1]],
                              area3 = growthTime_hr[[1]],
                              n12 = carbon_Mg_data[[1]],
                              n23 = Mg_growth_data[[1]],
                              n13 = carbon_growth_data[[1]],
                              n123 = carbon_Mg_growth_data[[1]],
                              category = c("Carbon Source","Mg Strees","Growth phase"),
                              fill = c("#decbe4","#b3cde3","#ccebc5"),
                              cex = rep(3,7),
                              cat.cex = rep(3.5,3),
                              margin=0.05)
dev.off()
###*****************************


###*****************************
# Gene names for Venn Diagram
test_combined_result %>%
  dplyr::filter(test_condition %in% c("carbonSource","Na_mM","Mg_mM","growthTime_hr"))%>%
  dplyr::filter(p_value<0.01)%>%
  dplyr::select(ID, gene_name,test_condition)%>%
  dplyr::mutate(exist=1)%>%
  tidyr::spread(key = test_condition, value = exist, fill = 0)%>%
  dplyr::group_by(carbonSource,Na_mM,Mg_mM,growthTime_hr)->test_combined_result_tidy

test_combined_result %>%
  dplyr::filter(test_condition %in% c("carbonSource","Na_mM","Mg_mM","growthTime_hr"))%>%
  dplyr::select(ID, gene_name,test_condition,p_related)%>%
  dplyr::mutate(score=sprintf("%.3f",p_related))%>%
  dplyr::select(-p_related)%>%
  dplyr::mutate(test_condition=paste0(test_condition,"_p"))%>%
  tidyr::spread(key = test_condition, value = score, fill = 0)%>%
  dplyr::select(ID,gene_name,
                carbonSource_p, Na_mM_p, Mg_mM_p, growthTime_hr_p)->test_combined_result_tidy2

test_combined_result_tidy %>%
  dplyr::filter(carbonSource==1 & Na_mM==1 & Mg_mM==1 & growthTime_hr==1)%>%
  dplyr::left_join(.,test_combined_result_tidy2)->CNMG_list

write.csv(x = CNMG_list, file = "../Tables/CNMG_list.csv",row.names = F,quote = T)

test_combined_result_tidy %>%
  dplyr::filter(carbonSource==1 & Na_mM==1 & Mg_mM==1 & growthTime_hr==0)%>%
  dplyr::left_join(.,test_combined_result_tidy2)->CNM_list

write.csv(x = CNM_list, file = "../Tables/CNM_list.csv",row.names = F,quote = T)

test_combined_result_tidy %>%
  dplyr::filter(carbonSource==1 & Na_mM==1 & Mg_mM==0 & growthTime_hr==1)%>%
  dplyr::left_join(.,test_combined_result_tidy2)->CNG_list

write.csv(x = CNG_list, file = "../Tables/CNG_list.csv",row.names = F,quote = T)

test_combined_result_tidy %>%
  dplyr::filter(carbonSource==1 & Na_mM==0 & Mg_mM==1 & growthTime_hr==1)%>%
  dplyr::left_join(.,test_combined_result_tidy2)->CMG_list
write.csv(x = CMG_list, file = "../Tables/CMG_list.csv",row.names = F,quote = T)

test_combined_result_tidy %>%
  dplyr::filter(carbonSource==0 & Na_mM==1 & Mg_mM==1 & growthTime_hr==1)%>%
  dplyr::left_join(.,test_combined_result_tidy2)->NMG_list

write.csv(x = NMG_list, file = "../Tables/NMG_list.csv",row.names = F,quote = T)
###*****************************


###*****************************
# Prepeare data frame for heatmap
test_combined_result %>%
  dplyr::mutate(p.adj_related_ws=p.adj_related*p_related_sign)%>%
  dplyr::select(ID,gene_name,test_condition,p.adj_related_ws)%>%
  tidyr::spread(test_condition, p.adj_related_ws)%>%
  dplyr::mutate(carbonSource=glucose+gluconate+lactate+glycerol)->combinedHeatMap
###*****************************


###*****************************
#remove columns with na
listOfRowsWithNA=which(rowSums(is.na(combinedHeatMap)) > 0,)
NumRowsWithNa=length(listOfRowsWithNA)
if (NumRowsWithNa!=0){
combinedHeatMap<-combinedHeatMap[-as.vector(listOfRowsWithNA), ]}
###*****************************

# combinedHeatMap %>% dplyr::select(ID,gene_name)->combinedHeatMap2
# for(counter01 in 3: ncol(combinedHeatMap)){
#   inputVector=t(as.vector(combinedHeatMap[,counter01]))
#   x=match(inputVector,sort(inputVector))
#   
#   combinedHeatMap2=cbind(combinedHeatMap2,as.data.frame(x))
# }
# colnames(combinedHeatMap2)<-colnames(combinedHeatMap)
# combinedHeatMap=combinedHeatMap2

row.names(combinedHeatMap)<-combinedHeatMap$gene_name
combinedHeatMap%>%dplyr::select(-ID,-gene_name,
                                -carbonSource, 
                                -exponential, 
                                -stationary, 
                                -late_stationary)->combinedHeatMap
###*****************************


###*****************************
# Do Dendogram for data sets
dist_x=dist(t(combinedHeatMap), 
            method = "manhattan", 
            diag = FALSE, 
            upper = FALSE, 
            p = 2)

METree = flashClust(dist_x, method = "complete");
METree = as.dendrogram(METree)

#Reordering
METree.reorder=rev(reorder(METree, c(1,2,3,4,5,8,9,10,11,6,7)))

ddata <- dendro_data(METree.reorder, type = "rectangle")
realOrderText.x<-as.vector(ddata$labels$label) # this real order is important for everything

fig04a<-ggplot(segment(ddata),aes(x = x, y = y, xend = xend, yend = yend)) + 
  geom_segment() + 
  theme_classic() +
  scale_x_discrete(labels=ddata$labels$label, expand = c(0,.5))+
  scale_y_log10(expand = c(0,0))+
  xlab("")+
  ylab("Height")+
  #ggtitle("Clustering of Data Sets")+
  theme( axis.text.x=element_text(size=10,angle = 90, hjust = 1),
         axis.text.y=element_blank(),
         axis.title.x=element_text(size=16),
         axis.title.y=element_text(size=16),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14),
         axis.line.y=element_blank(),
         axis.ticks.y=element_blank())
print(fig04a)
###*****************************


###*****************************
# Do Dendogram for genes
dist_x=dist(combinedHeatMap, 
            method = "euclidian", 
            diag = FALSE, 
            upper = FALSE, 
            p = 2)


METree = flashClust(dist_x, method = "complete");
METree = as.dendrogram(METree)

#Reordering
METree.reorder=rev(METree)

ddata <- dendro_data(METree.reorder, type = "rectangle")
realOrderText.y<-as.vector(ddata$labels$label) # this real order is important for everything

fig04b<-ggplot(segment(ddata),aes(x = y, y = x, xend = yend, yend = xend)) + 
  geom_segment() + 
  theme_classic() +
  scale_x_reverse(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  coord_cartesian(xlim = c(5, 60))+  #Changed from 25 (like in the mRNA) to 15, indicates that the proteins are better clustered
  ylab("")+
  xlab("")+
  ggtitle("Clustering of Proteins")+
  theme( axis.text.y=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14),
         axis.line.y=element_blank(),
         axis.ticks.y=element_blank())
print(fig04b)
###*****************************



###*****************************
# Reordering Genes
temp=rownames(combinedHeatMap)
combinedHeatMap=combinedHeatMap[,realOrderText.x]
combinedHeatMap %>%
  dplyr::mutate(gene_names=temp)->combinedHeatMap
combinedHeatMap<-orderFunction(combinedHeatMap,
                               column="gene_names",
                               target=realOrderText.y)
combinedHeatMap %>% 
  dplyr::mutate(order_y=seq(1,nrow(combinedHeatMap))) %>%
  tidyr::gather_("test_condition", "p_related", realOrderText.x)%>%
  dplyr::mutate(p_rescale=sign(p_related)*log10(abs(p_related)))->combinedHeatMapTidy


combinedHeatMapTidy$test_condition=factor(combinedHeatMapTidy$test_condition, levels = realOrderText.x)
###*****************************


###*****************************
fig04c=ggplot(combinedHeatMapTidy, aes( y=order_y,x=test_condition )) +
  geom_tile(aes(fill=p_related))+
  scale_fill_gradientn(colours=rev(c("#a50026",
                                     "#d73027",
                                     "#f46d43",
                                     "#fdae61",
                                     "#fee090",
                                     "#e0f3f8",
                                     "#abd9e9",
                                     "#74add1",
                                     "#4575b4",
                                     "#313695")),
                       values=rescale(c(-20,-5,-3,-1,-.1,.1,1,3,5,20)),
                       limits=c(-35,35),
                       guide = guide_colorbar(title = "Score"))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0),
                   labels = c("glucose" = "glucose",
                              "glycerol" = "glycerol",
                              "lactate" = "lactate",
                              "gluconate" = "gluconate",
                              "Na_mM" = "Na",
                              "Mg_mM"="Mg",
                              "growthTime_hr"="growth time"))+
  ylab("Genes")+
  xlab("Test Condition")+
  #ggtitle("Clustering between conditions")+
  theme_classic()+
  theme( axis.text.x=element_text(size=16,angle = 90, hjust = 1),
         axis.text.y=element_blank(),legend.title=element_text(size=18),
         axis.title.x=element_text(size=18),
         axis.title.y=element_text(size=18),
         legend.text=element_text(size=16),
         axis.line.x=element_blank(),
         axis.line.y=element_blank(),
         axis.ticks.x=element_blank(),
         axis.ticks.y=element_blank(),
         plot.title = element_text(size =18, face = "bold"))
print(fig04c)



# Combine figure4 
figComb=fig04c
dendogram.x <- gtable_filter(ggplotGrob(fig04a), "panel")
dendogram.y <- gtable_filter(ggplotGrob(fig04b), "panel")

g.main <- ggplotGrob(figComb) # convert the main plot into gtable

# add dendogram.x
index <- subset(g.main$layout, name == "panel") # locate where we want to insert the images
# add a row, one as spacer and one to take the images
g.main <- gtable_add_rows(g.main, unit.c(unit(0.35, "null")), index$b-1)
# add the grob that holds the images
g.main <- gtable_add_grob(g.main, dendogram.x, t = index$b, l = index$l, 
                          b = index$b, r = index$r, name="dendogram-x")
# add dendogram.y
# add a row, one as spacer and one to take the images
g.main <- gtable_add_cols(g.main, unit.c(unit(0.1, "null")), index$l-1)
index <- subset(g.main$layout, name == "panel")
# add the grob that holds the images
g.main <- gtable_add_grob(g.main, dendogram.y, 
                          t = index$b, 
                          l = index$l-1, 
                          b = index$b,
                          r = index$r-1, 
                          name="dendogram-y")


figComb04=ggdraw(g.main)
print(figComb04)

g.main$layout
###****************************


###****************************
# Draw bar graphs
test_combined_result %>%
  dplyr::filter(test_condition %in% c("glucose","gluconate","lactate","glycerol",
                                      "growthTime_hr","Na_mM","Mg_mM")) %>%
  dplyr::filter(p.adj < 0.05) -> test_combined_result_filtered

test_combined_result_filtered%>%
  dplyr::group_by(test_condition)%>%
  dplyr::summarize(data_length=length(test_condition))%>%
  dplyr::filter(data_length<=2)->temp01

test_combined_result_filtered%>%
  dplyr::group_by(test_condition, p_related_sign)%>%
  dplyr::summarize(data_length=length(test_condition))

if(nrow(temp01)>0){
test_combined_result_filtered%>%
  filter(test_condition!=as.vector(temp01$test_condition))->test_combined_result_filtered}

fig05 <- ggplot(test_combined_result_filtered, 
                       aes(x = test_condition, 
                           fill = as.factor(p_related_sign) )) + 
  geom_bar( position = position_dodge()) +
  xlab("Test Condition") + ylab("Gene Count") +
  scale_x_discrete(labels = c("glucose" = "Glucose",
                              "glycerol" = "Glycerol",
                              "lactate" = "Lactate",
                              "gluconate" = "Gluconate",
                              "Na_mM" = "Na",
                              "Mg_mM"="Mg",
                              "growthTime_hr"="Growth Time"))+
  annotation_logticks(sides = "l",long = unit(0.2, "cm")) + 
  scale_y_log10(expand = c(0, 0),limits= c(1,5000))+
  #ggtitle("Number of significantly changed genes") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.minor.x = element_blank())+
  theme_bw()+
  scale_fill_manual(values = c("blue","red"), 
                    name="Gene Regulation",
                    breaks=c("-1", "1"),
                    labels=c("Down Regulated", "Up Regulated"))

fig05


###****************************
# Combine plots
figComb<-plot_grid(fig05, figComb04)
###****************************


###****************************
save_plot(paste0("../initialPaper03f/violin_",generalFigureName,".png")
          ,fig01b, ncol = 1, nrow=3)
save_plot(paste0("../initialPaper03f/barGraph_",generalFigureName,".png")
          ,fig05, ncol = 1.2)
save_plot(paste0("../initialPaper03f/histogram_",generalFigureName,".png")
          ,figComb04, ncol = 1.3,nrow=1.5)
###****************************