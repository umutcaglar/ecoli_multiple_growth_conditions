# Calculate significant changes

###*****************************
# INITIAL COMMANDS TO RESET THE SYSTEM

rm(list = ls())
if (is.integer(dev.list())){dev.off()}
cat("\014")

# Set Working Directory
setwd('/Users/umut/GitHub/ecoli_multiple_growth_conditions/b_code_histogram_RNA&Protein/') # mac computer
###*****************************


###*****************************
# INSTALL LIBRARIES
library("dplyr")
library("tidyr")
require("reshape2")
library("lazyeval")
require("flashClust")

# For Plotting
library("ggplot2")
library("RColorBrewer")
library("grid")
library("gridExtra")
library("cowplot")
require("ggdendro")
require("scales")
require("gtable")
###*****************************


###*****************************
# PARAMETERS FOR DATA
conditionNumberChoice="uniqueCondition"
dataTypeChoice="protein_wo_NA" # can be "rna", "mrna","protein","protein_wo_NA","protein_wo_NAx6"
badDataFilterSetChoice="set00" # "set00", "set01", "set02"
dataTimeChoice="wholeSet" # exponential / stationary / late_stationary / wholeSet
MgLevelChoice="allMg" # allMg highMg midMg lowMg
NaLevelChoice="allNa" # allNa lowNa highNa
carbonTypeChoice="SY" # a letter combination from the list "SYAN"  
#S (glucose), Y (glycerol), A (lactate), N (gluconate)
filterTypeChoice="noFilter" # mean/sd/ max/ noFilter
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
step03=paste0("normalized_",normalizationMethodChoice,"_",deSeqNormChoice,"_",
              filterTypeChoice,"_",longNameList) 
conditionName=paste0("condition_",normalizationMethodChoice,"_",deSeqNormChoice,"_",
                     filterTypeChoice,"_",longNameList) 
###*****************************


###*****************************
# Load data
load(file = paste0("../a_results/",step03,".RData")) #Used for parameters data
assign(x = "mainDataFrame",value = get(step03))
###*****************************


###*****************************
# Calculate the average of repeated conditions
mainDataFrame<-as.data.frame(mainDataFrame)
mainDataFrame %>% mutate(gene_id=row.names(mainDataFrame))->mainDataFrameTidy
endVal=ncol(mainDataFrame)
mainDataFrameTidy %>%
  tidyr::gather(condition,read,1:endVal)%>%
  dplyr::mutate(conditionShort=gsub("_[[:digit:]]$","",condition))->mainDataFrameTidy
mainDataFrameTidy %>% 
  group_by(gene_id,conditionShort)%>%
  dplyr::summarise(read=mean(read))->mainDataFrameSummary

mainDataFrameSummary %>%
  dplyr::group_by()%>%  
  tidyr::spread(key = conditionShort, value = read)->mainDataFrameSummary

row.names(mainDataFrameSummary)<-mainDataFrameSummary$gene_id
mainDataFrameSummary %>%
  dplyr::select(-gene_id)->mainDataFrameSummary
mainDataFrame=as.matrix(mainDataFrameSummary)
###*****************************


###*****************************
#Reduce number of rows in condition and get rid of repeats
condition<-get(conditionName)
condition=condition[which(substr(as.vector(condition$dataSet),10,10) %in% c(0)),]
condition$dataSet<-gsub("_[[:digit:]]$","",as.vector(condition$dataSet))
condition$sampleNum<-gsub("_[[:digit:]]$","",as.vector(condition$sampleNum))
###*****************************
browser()

###*****************************
# Generate a new column to condition df that uniquely defines the condition named dataSet2
condition %>%
  dplyr::mutate(dataSet2d=ifelse(carbonSource=="glucose","S",NA)) %>%
  dplyr::mutate(dataSet2d=ifelse(carbonSource=="glycerol","Y",dataSet2d)) %>%
  dplyr::mutate(dataSet2d=ifelse(carbonSource=="lactate","A",dataSet2d)) %>%
  dplyr::mutate(dataSet2d=ifelse(carbonSource=="gluconate","N",dataSet2d))->condition

condition %>%
  dplyr::mutate(dataSet2a=ifelse(growthPhase=="exponential","exp",NA)) %>%
  dplyr::mutate(dataSet2a=ifelse(growthPhase=="stationary","sta",dataSet2a)) %>%
  dplyr::mutate(dataSet2a=ifelse(growthPhase=="late_stationary","lat",dataSet2a)) ->condition

condition %>%
  dplyr::mutate(dataSet2e=ifelse(experiment=="glucose_time_course","Stc",NA)) %>%
  dplyr::mutate(dataSet2e=ifelse(experiment=="glycerol_time_course","Ytc",dataSet2e)) %>%
  dplyr::mutate(dataSet2e=ifelse(experiment=="NaCl_stress","Nas",dataSet2e)) %>%
  dplyr::mutate(dataSet2e=ifelse(experiment=="lactate_growth","Agr",dataSet2e)) %>%
  dplyr::mutate(dataSet2e=ifelse(experiment=="gluconate_growth","Ngr",dataSet2e)) %>%
  dplyr::mutate(dataSet2e=ifelse(experiment=="MgSO4_stress_low","Mgl",dataSet2e)) %>%
  dplyr::mutate(dataSet2e=ifelse(experiment=="MgSO4_stress_high","Mgh",dataSet2e)) -> condition


# "Stc" for "glucose_time_course", 
# "Ytc" for "glycerol_time_course", 
# "Nas" for "NaCl_stress", 
# "Agr" for "lactate_growth", 
# "Ngr" for "gluconate_growth", 
# "Mgl" for "MgSO4_stress_low", 
# "Mgh" for "MgSO4_stress_high" 


condition %>%
  mutate(dataSet2=paste0(dataSet,"_",
                         dataSet2a,"_",
                         Mg_mM_Levels,"_",
                         Na_mM_Levels,"_",
                         dataSet2d,"_",
                         dataSet2e))->condition
###*****************************


###*****************************??
# Change Colnames of main data frame
#colnames(mainDataFrame)<-condition$dataSet2
###*****************************


###*****************************
# Find Ideal Order for dataSet Dendogram
condition %>%
  dplyr::select(dataSet,
                growthPhase, 
                carbonSource, 
                Mg_mM_Levels, 
                Na_mM_Levels,
                batchNumber) %>%
  dplyr::mutate(orderNo=1:nrow(condition))->conditionSummary


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


conditionSummary %>%
  dplyr::group_by()->conditionSummary
conditionSummaryF<-conditionSummary

# there is an order in the application of the 4 commands below
# the order in here is 
# Mg, growth phase, NA, carbons source from less ordered to most ordered
# Determined by code "cophenetic_distance_main_final.R"


Mg_mM_LevelsTarget <- c("lowMg", "baseMg", "highMg")
conditionSummary<-orderFunction(conditionSummary,
                                column="Mg_mM_Levels",
                                target=Mg_mM_LevelsTarget)


growthPhaseTarget <- c("exponential", "stationary", "late_stationary")
conditionSummary<-orderFunction(conditionSummary,
                                column="growthPhase",
                                target=growthPhaseTarget)

Na_mM_LevelsTarget <- c("baseNa", "highNa")
conditionSummary<-orderFunction(conditionSummary,
                                column="Na_mM_Levels",
                                target=Na_mM_LevelsTarget)

carbonSourceTarget <- c("glucose", "glycerol") #, "lactate", "gluconate")
conditionSummary<-orderFunction(conditionSummary,
                                column="carbonSource",
                                target=carbonSourceTarget)



desiredOrderNo.x=as.vector(conditionSummary$orderNo)
desiredOrderText.x=as.vector(conditionSummary$dataSet)
###*****************************


###*****************************
# Do Dendogram for data sets
dist_x=dist(t(mainDataFrame), 
            method = "euclidean", 
            diag = FALSE, 
            upper = FALSE, 
            p = 2)

METree_x = flashClust(dist_x, method = "complete");
METree_x_d = as.dendrogram(METree_x)

#Reordering
METree.reorder_x_d=rev(reorder(METree_x_d, desiredOrderNo.x))
save(list = c("METree_x","METree.reorder_x_d","dist_x","conditionSummary"),
     file = paste0("../b_results/treeFile_",step03,".RData"))


ddata <- dendro_data(METree.reorder_x_d, type = "rectangle")
realOrderText.x<-as.vector(ddata$labels$label) # this real order is important for everything

fig01a<-ggplot(segment(ddata),aes(x = x, y = y, xend = xend, yend = yend)) + 
  geom_segment() + 
  theme_classic() +
  scale_x_discrete(labels=ddata$labels$label, expand = c(0,.5))+
  scale_y_continuous(expand = c(0,0))+
  xlab("")+
  ylab("Height")+
  ggtitle("Clustering of Data Sets")+
  theme( axis.text.x=element_text(size=10,angle = 90, hjust = 1,vjust=0.5),
#         axis.text.y=element_blank(),
         axis.title.x=element_text(size=16),
         axis.title.y=element_text(size=16),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14),
         axis.line.y=element_blank(),
         axis.ticks.y=element_blank())
print(fig01a)
###*****************************


###*****************************
# Do Color Code Main X Axis 
dataSetTarget <- realOrderText.x
conditionSummary<-orderFunction(conditionSummary,
                                column="dataSet",
                                target=dataSetTarget)
conditionSummary %>% 
  dplyr::mutate(orderNoCurrent=seq(1,nrow(conditionSummary)))->conditionSummary

conditionSummary %>%
  tidyr::gather(columnName,condition,growthPhase:Na_mM_Levels)->conditionSummaryTidy

conditionSummaryTidy$condition <- factor(conditionSummaryTidy$condition, levels = 
                                           c("glucose","glycerol",
                                             "baseNa","highNa",
                                             "exponential","stationary","late_stationary",
                                             "lowMg","baseMg","highMg"))
conditionSummaryTidy$columnName <- factor(conditionSummaryTidy$columnName, levels = 
                                            rev(c("carbonSource", "Na_mM_Levels", "growthPhase","Mg_mM_Levels")))

listColors=c("#bcbddc","#9e9ac8",
             "#fdbe85","#fd8d3c",
             "#bae4b3","#74c476","#238b45",
             "#bdd7e7","#6baed6","#2171b5")


fig02a<-ggplot(conditionSummaryTidy, aes( y=columnName,x=factor(orderNoCurrent)))+
  geom_tile(aes(fill=condition), color="black")+
  #geom_text(aes(label=orderNo,angle = 90))+
  scale_fill_manual(values = listColors)+
  scale_y_discrete(expand = c(0,0), 
                   labels = rev(c("Carbon source","Na levels","Growth Phase","Mg levels"))) +
  scale_x_discrete(labels=as.vector(conditionSummary$dataSet),expand = c(0,0))+
  guides(fill = guide_legend(override.aes = list(colour = NULL),
                             nrow=2,byrow = TRUE))+
  theme(axis.text.x = element_text(angle = 90, hjust = -1, size = 6),
        axis.text.y = element_text(size=14, face = "bold"),
        legend.title=element_text(size=18),
        legend.text=element_text(size=18),
        axis.title=element_blank(),
        axis.line = element_blank(),
        axis.ticks.y=element_blank(),
        legend.position="bottom",
        legend.key.size= unit(.6,"cm"))

print(fig02a)
###*****************************


###*****************************
# Do Dendogram for proteins
dist_y=dist(mainDataFrame, 
            method = "euclidean", 
            diag = FALSE, 
            upper = FALSE, 
            p = 2)

METree_y = flashClust(dist_y, method = "complete");
METree_y_d = rev(as.dendrogram(METree_y))
ddata <- dendro_data(METree_y_d, type = "rectangle")
realOrderText.y<-as.vector(ddata$labels$label) # this real order is important for everything

fig03a<-ggplot(segment(ddata),aes(x = y, y = x, xend = yend, yend = xend)) + 
  geom_segment() + 
  theme_classic() +
  scale_x_reverse(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  coord_cartesian(xlim = c(15, 175))+  #Changed from 25 (like in the mRNA) to 15, indicates that the proteins are better clustered
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
print(fig03a)

###*****************************


###*****************************
# Do Main Data
mainDataFrame<-mainDataFrame[realOrderText.y,realOrderText.x]
mainDataFrame<-as.data.frame(mainDataFrame)
dataSetList=colnames(mainDataFrame)
mainDataFrame %>% 
  dplyr::mutate(order_y=seq(1,nrow(mainDataFrame))) %>%
  tidyr::gather_("dataSet", "read", dataSetList)->mainDataFrameTidy

mainDataFrameTidy$dataSet=factor(mainDataFrameTidy$dataSet, levels = dataSetList)



fig04a=ggplot(mainDataFrameTidy, aes( y=order_y,x=dataSet )) +
  geom_tile(aes(fill=read))+
  scale_fill_gradientn(colours=c("#4575b4","#91bfdb","#e0f3f8","#fee090","#fc8d59","#d73027"),
                       values=rescale(c(0,1.5,2.5,3,5,10)),
                       limits=c(0,12),
                       guide = guide_colorbar(title = "Log \nCount"))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))+
  ylab("Proteins")+
  xlab("Data Sets")+
  theme( axis.text.x = element_blank(),
         axis.text.y=element_blank(),legend.title=element_text(size=18),
         #axis.title.x=element_text(size=18),
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=18),
         legend.text=element_text(size=16),
         axis.line.x=element_blank(),
         axis.line.y=element_blank(),
         axis.ticks.x=element_blank(),
         axis.ticks.y=element_blank())



fig04b=ggplot(mainDataFrameTidy, aes( y=order_y,x=dataSet )) +
  geom_tile(aes(fill=read))+
  scale_fill_gradientn(colours=c("#4575b4","#91bfdb","#e0f3f8","#fee090","#fc8d59","#d73027"),
                       values=rescale(c(0,1.5,2.5,3,5,10)),
                       limits=c(0,15),
                       guide = guide_colorbar(title = "Log \nCount"))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))+
  xlab("Data Sets")+
  theme_classic()+
  theme( axis.text.x = element_text(angle = 90, hjust = -1, size = 6),
         axis.text.y=element_blank(),
         axis.title.x=element_text(size=18),
         axis.title.y=element_blank(),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14),
         axis.line.y=element_blank(),
         axis.ticks.y=element_blank())

###*****************************


###*****************************
# Combine Figures
figComb=fig04a
dendogram.x <- gtable_filter(ggplotGrob(fig01a), "panel")
color_x <- gtable_filter(ggplotGrob(fig02a), "panel")
color_x_legend <- gtable_filter(ggplotGrob(fig02a), "guide-box")
color_x_categories <- gtable_filter(ggplotGrob(fig02a), "axis-l")
dendogram.y <- gtable_filter(ggplotGrob(fig03a), "panel")
x.axis_text <- gtable_filter(ggplotGrob(fig04b), "xlab")

g.main <- ggplotGrob(figComb) # convert the main plot into gtable

# add color_x
index <- subset(g.main$layout, name == "panel") # locate where we want to insert the images
# add a row, one as spacer and one to take the images
g.main <- gtable_add_rows(g.main, unit.c(unit(0.1, "null")), index$b-1)
# add the grob that holds the images
g.main <- gtable_add_grob(g.main, color_x, t = index$b, l = index$l, 
                          b = index$b, r = index$r, name="colors")
g.main <- gtable_add_rows(g.main, unit.c(unit(0.1, "in")), index$b) #************


# add dendogram.x
# add a row, one as spacer and one to take the images
g.main <- gtable_add_rows(g.main, unit.c(unit(0.15, "null")), index$b-2)
# add the grob that holds the images
g.main <- gtable_add_grob(g.main, dendogram.x, t = index$b-1, l = index$l, 
                          b = index$b-1, r = index$r, name="dendogram-x")

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

# add color_x_categories
index <- subset(g.main$layout, name == "colors") # locate where we want to insert the images
g.main <- gtable_add_grob(g.main, color_x_categories, t = index$b, l = index$l-1, 
                          b = index$b, r = index$r-1, name="color_x_categories")

#add color_legend
#add a row, one as spacer and one to take the images
index <- subset(g.main$layout, name == "panel")
g.main <- gtable_add_rows(g.main, unit.c(unit(0.4, "in")), index$b+3)
g.main <- gtable_add_grob(g.main, color_x_legend, 
                          t = index$b+3, 
                          l = index$l, 
                          b = index$b+3, 
                          r = index$r, name="color_legend")

#add space above color legend
#add a row as spacer 
index <- subset(g.main$layout, name == "color_legend")
g.main <- gtable_add_rows(g.main, unit.c(unit(0.25, "in")), index$b-3)

#add x_axis text
#add a row, one as spacer and one to take the images
index <- subset(g.main$layout, name == "dendogram-x")

g.main <- gtable_add_rows(g.main, unit.c(unit(0.3, "in")), index$b-1)
g.main <- gtable_add_grob(g.main, x.axis_text, 
                          t = index$b-1, 
                          l = index$l, 
                          b = index$b-1, 
                          r = index$r, name="x_axis_name")

#add space above x_axis text
#add a row as spacer
index <- subset(g.main$layout, name == "dendogram-x")
g.main <- gtable_add_rows(g.main, unit.c(unit(0.3, "in")), index$b-3)



figComb=ggdraw(g.main)

figureName=paste0("../b_figures/","heatmap_",step03,".png")
cowplot::save_plot(plot = figComb, filename = figureName,ncol = 4,nrow = 3, dpi=300)

g.main$layout

