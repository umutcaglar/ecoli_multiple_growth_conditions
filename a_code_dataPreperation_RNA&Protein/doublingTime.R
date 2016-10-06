# Doubling time associated figures

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
{setwd(paste0("/Users/umut/GitHub/ecoli_multiple_growth_conditions/a_code_dataPreperation_RNA&Protein/"))} # mac computer
###*****************************


###*****************************
# DOWNLOAD LIBRARIES
require("dplyr")
require("tidyr")
require("stringr")

require("ggplot2")
require("RColorBrewer")
require("scales")
require("cowplot")
require("ggrepel")
###*****************************


###*****************************
# Source functions
source("../e_Flux_Analyze/replaceFunction.R")
###*****************************


###*****************************
# install data
doublingTimeData<-read.csv(file = "../DoublingTimeRawData/AG3C_doubling_times_raw.csv")
###*****************************


###*****************************
# add a column defining base and several tests
doublingTimeData %>% dplyr::mutate(experiment=ifelse(name %in% c("Gluconate.tab", "Glucose.tab", "Glycerol.tab", "Lactate.tab"),
                                          "carbonSource", "undecided" ),
                           experiment=ifelse(grepl(pattern = "MgSO4",x = name),
                                          "Mg_concentration", experiment),
                           experiment=ifelse(grepl(pattern = "NaCl",x = name),
                                          "Na_concentration", experiment),
                           experiment=ifelse(name %in% c("Glucose.tab", "MgSO4_000.800_mM.tab", "NaCl_005_mM.tab"),
                                          "base" ,experiment))->doublingTimeData
###*****************************


###*****************************
#Add relevant concentrations 
doublingTimeData %>% dplyr::mutate(name=gsub(pattern = ".tab",replacement = "",x = name))->doublingTimeData
doublingTimeData %>% dplyr::mutate(carbonSource=ifelse(experiment!="carbonSource","glucose","not known"),
                           carbonSource=ifelse(name=="Gluconate","gluconate", carbonSource),
                           carbonSource=ifelse(name=="Glycerol","glycerol", carbonSource),
                           carbonSource=ifelse(name=="Lactate","lactate", carbonSource))->doublingTimeData

doublingTimeData %>% dplyr::mutate(Mg_mM=ifelse(experiment!="Mg_concentration",0.8,NA),
                           Mg_mM=ifelse(name=="MgSO4_000.080_mM" , .08, Mg_mM),
                           Mg_mM=ifelse(name=="MgSO4_000.800_mM" , .8, Mg_mM),
                           Mg_mM=ifelse(name=="MgSO4_008.000_mM" , 8, Mg_mM),
                           Mg_mM=ifelse(name=="MgSO4_050.000_mM" , 50, Mg_mM),
                           Mg_mM=ifelse(name=="MgSO4_200.000_mM" , 200, Mg_mM),
                           Mg_mM=ifelse(name=="MgSO4_400.000_mM" , 400, Mg_mM),
                           Mg_mM=ifelse(name=="MgSO4-2_000.005_mM" , .005, Mg_mM),
                           Mg_mM=ifelse(name=="MgSO4-2_000.010_mM" , .01, Mg_mM),
                           Mg_mM=ifelse(name=="MgSO4-2_000.020_mM" , .02, Mg_mM),
                           Mg_mM=ifelse(name=="MgSO4-2_000.040_mM" , .04, Mg_mM),
                           Mg_mM=ifelse(name=="MgSO4-2_000.080_mM" , .08, Mg_mM))->doublingTimeData


doublingTimeData %>% dplyr::mutate(Na_mM=ifelse(experiment!="Na_concentration",5,NA),
                           Na_mM=ifelse(name=="NaCl_100_mM" , 100, Na_mM),
                           Na_mM=ifelse(name=="NaCl_200_mM" , 200, Na_mM),
                           Na_mM=ifelse(name=="NaCl_300_mM" , 300, Na_mM))->doublingTimeData
###*****************************


###*****************************
doublingTimeData %>%
  dplyr::mutate(name=gsub("-2","",name))->doublingTimeData

# summarize data
doublingTimeData %>% 
  dplyr::mutate(name_base=name)%>%
  dplyr::mutate(name_base=ifelse(experiment=="base","base",as.character(name))) %>%
  dplyr::group_by(name_base) %>%
  dplyr::summarise(experiment=unique(experiment),
                   carbonSource=unique(carbonSource),
                   Mg_mM=unique(Mg_mM),
                   Na_mM=unique(Na_mM),
                   meanDT=mean(doubling.time.minutes),
                   stdDT=sd(doubling.time.minutes),
                   rep=n(),
                   stderrDT=sd(doubling.time.minutes)/sqrt(n()))->doublingTimeData_sum
###*****************************


###*****************************
# seperate the data frame into 3 data frames
doublingTimeData_sum %>% dplyr::filter(experiment %in% c("carbonSource", "base"))->carbonSourceData
doublingTimeData_sum %>% dplyr::filter(experiment %in% c("Mg_concentration", "base"))->mgStressData
doublingTimeData_sum %>% dplyr::filter(experiment %in% c("Na_concentration", "base"))->naStressData
###*****************************


###*****************************
#Figures associated with carbon source
carbonSourceData$carbonSource<-factor(carbonSourceData$carbonSource,
                                      levels=c("glucose","glycerol","lactate","gluconate"))

carbonSourceData%>%dplyr::filter(experiment=="base")%>%.$meanDT->base_y

sizeValue=2.5

fig01<-ggplot(carbonSourceData,aes(x=carbonSource,y=as.numeric(meanDT), colour=experiment))+
  #geom_vline(xintercept = as.numeric(1), color="orange", linetype = "longdash")+
  geom_hline(yintercept = base_y, 
             color="orange", linetype = "longdash")+
  scale_x_discrete()+
  geom_point(size=sizeValue)+
  scale_color_manual(values=c("Red", "Black"))+
  geom_errorbar(aes(ymax = meanDT+stderrDT, ymin=meanDT-stderrDT),width=0)+
  expand_limits(x = 0, y = 0)+
  ylim(0,120)+
  theme_bw()+
  xlab("Carbon sources")+
  ylab("Doubling time (min)")+
  theme(axis.line.y = element_blank(),
        legend.position="none",
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14))

print(fig01)


fig02<-ggplot(mgStressData,aes(x=Mg_mM,y=meanDT, colour=experiment))+
  #geom_vline(xintercept = .8, color="orange", linetype = "longdash")+
  geom_hline(yintercept = base_y, color="orange", linetype = "longdash")+
  geom_point(size=sizeValue)+
  scale_color_manual(values=c("Red", "Black"))+
  geom_errorbar(aes(ymax = meanDT+stderrDT, ymin=meanDT-stderrDT),width=0)+
  expand_limits(x = 0, y = 0)+
  ylim(0,120)+
  scale_x_log10(breaks=c(0.01,0.1,1,10,100), labels = comma)+
  theme_bw()+
  xlab("Mg Concentration mM")+
  ylab("Doubling time (min)")+
  theme(axis.line.y = element_blank(),
        legend.position="none",
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14))

print(fig02)



fig03<-ggplot(naStressData,aes(x=Na_mM,y=meanDT, colour=experiment))+
  #geom_vline(xintercept = 5, color="orange", linetype = "longdash")+
  geom_hline(yintercept = base_y, color="orange", linetype = "longdash")+
  geom_point(size=sizeValue)+
  scale_color_manual(values=c("Red", "Black"))+
  geom_errorbar(aes(ymax = meanDT+stderrDT, ymin=meanDT-stderrDT),width=0)+
  expand_limits(x = 0, y = 0)+
  ylim(0,120)+
  theme_bw()+
  xlab("Na Concentration mM")+
  ylab("Doubling time (min)")+
  theme(axis.line.y = element_blank(),
        legend.position="none",
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14))

print(fig03)
###*****************************


###*****************************
figAll<-cowplot::plot_grid(plotlist = list(fig01, fig02, fig03), 
                           align = "v", ncol = 1, scale = .95,labels = c("A", "B", "C"))

print(figAll)

cowplot::save_plot(filename = "../a_figures/duplicationTime.pdf",plot = figAll,ncol = 2,nrow = 3)
###*****************************