# Generate Heatmap by hand for mrna

###*****************************
# INITIAL COMMANDS TO RESET THE SYSTEM
rm(list = ls())
if (is.integer(dev.list())){dev.off()}
cat("\014")
###*****************************

###*****************************
# Set Working Directory
# One needs to arrange the correct pathway if this is not umut's computer ;)
if(as.vector(Sys.info()["effective_user"]=="umut"))
{setwd(paste0("/Users/umut/GitHub/ecoli_multiple_growth_conditions/",
              "b_code_histogram_RNA&Protein/"))} # mac computer
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
#Load Functions
source("../a_code_dataPreperation_RNA&Protein/data_naming_functions.R")
###*****************************


###*****************************
# Find the csv files that need to be imported
dataName=name_data(initialValue=c("resDf"), # can be c("genes0.05","genes_P0.05Fold2","resDf")
                   dataType = "mrna", # can be "rna", "mrna", "protein", "protein_wo_NA"
                   badDataSet = "set02", # can be "set00",set01","set02", "set03"
                   # referenceParameters can be a vector like
                   # c("growthPhase", "Mg_mM_Levels", "Na_mM_Levels", "carbonSource", "experiment")
                   referenceParameters=c("growthPhase",
                                         "Mg_mM_Levels", 
                                         "Na_mM_Levels", 
                                         "carbonSource", 
                                         "experiment"),
                   # referenceLevels can be a vector like
                   # c("exponential", "baseMg", "baseNa", "glucose", "glucose_time_course")
                   referenceLevels=c("exponential",
                                     "baseMg", 
                                     "baseNa", 
                                     "glucose", 
                                     "glucose_time_course"),
                   experimentVector = c("allEx"), # can be "Stc","Ytc","Nas","Agr","Ngr","Mgl","Mgh" // "allEx"
                   carbonSourceVector = "SYAN", # can be any sub combination of "SYAN"
                   MgLevelVector = c("allMg"), # can be "lowMg","baseMg","highMg" // "allMg"
                   NaLevelVector = c("allNa"), # can be "baseNa","highNa" // "allNa"
                   growthPhaseVector = c("allPhase"), # can be "exponential","stationary","late_stationary" // "allPhase"
                   filterGenes = "noFilter", # can be "noFilter", "meanFilter", "maxFilter", "sdFilter" 
                   threshold=NA, # the threshold value for "meanFilter", "maxFilter", "sdFilter"
                   roundData=TRUE,
                   sumTechnicalReplicates=TRUE,
                   deSeqSfChoice="p1Sf", # can be "regSf", "p1Sf"
                   normalizationMethodChoice= "vst", # can be "vst", "rlog", "log10", "noNorm"
                   test_for = "noTest")  # works only if normalizationMethodChoice == noNorm
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource", "noTest")

dataName=as.data.frame(dataName[1])

metaDataName=dataName
metaDataName$objectName.initial="metaData"
treeDataName=dataName
treeDataName$objectName.initial="treeData"
heatMapName=dataName
heatMapName$objectName.initial="heatMap"

dataName=paste(dataName,collapse = "_")
metaDataName=paste(metaDataName,collapse = "_")
treeDataName=paste(treeDataName,collapse = "_")
heatMapName=paste(heatMapName,collapse = "_")

mainDataFrame=read.csv(file = paste0("../a_results/",dataName,".csv"),header = TRUE,row.names = 1)
condition=read.csv(file = paste0("../a_results/",metaDataName,".csv"),header = TRUE)
###*****************************


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
# carbon source, Na, Mg, growth phase from less ordered to miost ordered

carbonSourceTarget <- c("glucose", "glycerol", "lactate", "gluconate")
conditionSummary<-orderFunction(conditionSummary,
                                column="carbonSource",
                                target=carbonSourceTarget)

Na_mM_LevelsTarget <- c("baseNa", "highNa")
conditionSummary<-orderFunction(conditionSummary,
                                column="Na_mM_Levels",
                                target=Na_mM_LevelsTarget)


Mg_mM_LevelsTarget <- c("lowMg", "baseMg", "highMg")
conditionSummary<-orderFunction(conditionSummary,
                                column="Mg_mM_Levels",
                                target=Mg_mM_LevelsTarget)


growthPhaseTarget <- c("exponential", "stationary", "late_stationary")
conditionSummary<-orderFunction(conditionSummary,
                                column="growthPhase",
                                target=growthPhaseTarget)



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
     file = paste0("../b_results/",treeDataName,".RData"))

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
         #axis.text.y=element_blank(),
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
                                           c("exponential","stationary","late_stationary",
                                             "lowMg","baseMg","highMg",
                                             "baseNa","highNa",
                                             "glucose","glycerol","lactate","gluconate"))
conditionSummaryTidy$columnName <- factor(conditionSummaryTidy$columnName, levels = 
                                            rev(c("growthPhase","Mg_mM_Levels","Na_mM_Levels","carbonSource")))

listColors=c("#bae4b3","#74c476","#238b45",
             "#bdd7e7","#6baed6","#2171b5",
             "#fdbe85","#fd8d3c",
             "#bcbddc","#9e9ac8","#807dba","#6a51a3")

fig02a<-ggplot(conditionSummaryTidy, aes( y=columnName,x= factor(orderNoCurrent)))+
  geom_tile(aes(fill=condition), color="black")+
  #geom_text(aes(label=orderNo,angle = 90))+
  scale_fill_manual(values = listColors)+
  scale_y_discrete(expand = c(0,0), 
                   labels = rev(c("Growth Phase","Mg levels","Na levels","Carbon source"))) +
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

fig02a
###*****************************


###*****************************
# Do Dendogram for genes
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
  coord_cartesian(xlim = c(25, 175))+
  ylab("")+
  xlab("")+
  ggtitle("Clustering of Genes")+
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
                       values=rescale(c(-1,0,.5,1.5,3,10)),
                       limits=c(-3,18),
                       guide = guide_colorbar(title = "Log \nCount"))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))+
  ylab("Genes")+
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
                       values=rescale(c(-1,0,.5,1.5,3,10)),
                       limits=c(-3,18),
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

figureName=paste0("../b_figures/",heatMapName,".png")
cowplot::save_plot(plot = figComb, filename = figureName,ncol = 4,nrow = 3, dpi=300)

g.main$layout
