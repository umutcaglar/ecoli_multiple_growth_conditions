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
metaData<-read.csv(file = "../a_results/metaData.csv")
###*****************************


###*****************************
#Figures associated with carbon source
metaData %>%
  dplyr::filter(experiment==c("glucose_time_course","glycerol_time_course",
                              "lactate_growth","gluconate_growth")) %>%
  dplyr::group_by(experiment) %>%
  dplyr::summarize(doublingTimeMinutes=unique(doublingTimeMinutes), 
                   doublingTimeMinutes.95m=unique(doublingTimeMinutes.95m), 
                   doublingTimeMinutes_95p=unique(doublingTimeMinutes_95p)) %>%
  dplyr::mutate(experiment2=replace_fun(input_vector = as.vector(experiment),
                                        initialVal = c("glucose_time_course","glycerol_time_course",
                                                       "lactate_growth","gluconate_growth"),
                                        finalVal = c("glucose","glycerol","lactate","gluconate")))-> glucoseTimeCourseData
glucoseTimeCourseData$experiment2<-factor(glucoseTimeCourseData$experiment2,
                                          levels=c("glucose","glycerol","lactate","gluconate"))

fig01<-ggplot(glucoseTimeCourseData,aes(x=experiment2,y=doublingTimeMinutes))+
  geom_point()+
  geom_errorbar(aes(ymax = doublingTimeMinutes_95p, ymin=doublingTimeMinutes.95m),width=0)+
  geom_vline(xintercept = 1, color="orange", linetype = "longdash")+
  geom_point()+
  geom_errorbar(aes(ymax = doublingTimeMinutes_95p, ymin=doublingTimeMinutes.95m),width=0)+
  expand_limits(x = 0, y = 0)+
  ylim(0,120)+
  theme_bw()+
  xlab("Carbon sources")+
  theme(axis.line.y = element_blank(),
        legend.position="bottom",
        axis.title.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14))

###*****************************



###*****************************
#Figures associated with Mg levels
metaData %>%
  dplyr::filter(experiment %in% c("MgSO4_stress_low", "MgSO4_stress_high")) %>%
  dplyr::group_by(experiment,Mg_mM) %>%
  dplyr::summarize(doublingTimeMinutes=unique(doublingTimeMinutes), 
                   doublingTimeMinutes.95m=unique(doublingTimeMinutes.95m), 
                   doublingTimeMinutes_95p=unique(doublingTimeMinutes_95p))->mgStressData

fig02<-ggplot(mgStressData,aes(x=Mg_mM,y=doublingTimeMinutes))+
  geom_point()+
  geom_errorbar(aes(ymax = doublingTimeMinutes_95p, ymin=doublingTimeMinutes.95m),width=0)+
  geom_vline(xintercept = .8, color="orange", linetype = "longdash")+
  geom_point()+
  geom_errorbar(aes(ymax = doublingTimeMinutes_95p, ymin=doublingTimeMinutes.95m),width=0)+
  expand_limits(x = 0, y = 0)+
  ylim(0,120)+
  scale_x_log10(breaks=c(0.01,0.1,1,10,100))+
  theme_bw()+
  xlab("Mg Concentration mM")+
  theme(axis.line.y = element_blank(),
        legend.position="bottom",
        axis.title.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14))
###*****************************


###*****************************
#Figures associated with Na levels
metaData %>%
  dplyr::filter(experiment %in% c("NaCl_stress")) %>%
  dplyr::group_by(Na_mM) %>%
  dplyr::summarize(doublingTimeMinutes=unique(doublingTimeMinutes), 
                   doublingTimeMinutes.95m=unique(doublingTimeMinutes.95m), 
                   doublingTimeMinutes_95p=unique(doublingTimeMinutes_95p))->naStressData

fig03<-ggplot(naStressData,aes(x=Na_mM,y=doublingTimeMinutes))+
  geom_point()+
  geom_errorbar(aes(ymax = doublingTimeMinutes_95p, ymin=doublingTimeMinutes.95m),width=0)+
  geom_vline(xintercept = 5, color="orange", linetype = "longdash")+
  geom_point()+
  geom_errorbar(aes(ymax = doublingTimeMinutes_95p, ymin=doublingTimeMinutes.95m),width=0)+
  expand_limits(x = 0, y = 0)+
  ylim(0,120)+
  theme_bw()+
  xlab("Na Concentration mM")+
  theme(axis.line.y = element_blank(),
        legend.position="bottom",
        axis.title.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14))
###*****************************

print(fig01)
print(fig02)
print(fig03)


figAll<-cowplot::plot_grid(plotlist = list(fig01, fig02, fig03), 
                           align = "v", ncol = 1, scale = .95,labels = c("A", "B", "C"))

print(figAll)

cowplot::save_plot(filename = "../a_figures/duplicationTime.pdf",plot = figAll,ncol = 2,nrow = 3)