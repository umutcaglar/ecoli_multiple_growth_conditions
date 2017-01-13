# the difference between controling on growth and not controlling on growth

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
{setwd('/Users/umut/GitHub/ecoli_multiple_growth_conditions/c_code_change_wrt_variables_RNA&Protein/')} # mac computer
###*****************************


###*****************************
# DOWNLOAD LIBRARIES
require("dplyr")
require("tidyr")

require("ggplot2")
require("cowplot")
require("VennDiagram")
require("gtable")
require("gridExtra")
###*****************************


###*****************************
# Read significantly altered gene list
data<-read.csv(file = "../c_results/combinedDifferentiallyExpressedGenes_DeSeq.csv")
data %>% 
  dplyr::mutate(test_for=ifelse( test = grepl(pattern = "gluconate", x = testVSbase) , yes =  "Glc", no = NA ),
                test_for=ifelse( test = grepl(pattern = "glycerol", x = testVSbase) , yes = "Gly", no = test_for), 
                test_for=ifelse( test = grepl(pattern = "lactate", x = testVSbase) , yes = "Lac", no = test_for),
                test_for=ifelse( test = grepl(pattern = "highMg", x = testVSbase) , yes = "High Mg", no = test_for),
                test_for=ifelse( test = grepl(pattern = "lowMg", x = testVSbase) , yes = "Low Mg", no = test_for),
                test_for=ifelse( test = grepl(pattern = "highNa", x = testVSbase) , yes = "High Na", no = test_for)) %>%
  dplyr::mutate(controlledForGrowth = ifelse (test = grepl(pattern = "doublingTimeMinutes", x = investigatedEffect), yes = 1 , no = 0)) -> data 
###*****************************

###*****************************
# Divide data into 2 pieces
data %>%
  dplyr::filter(controlledForGrowth==0) %>%
  dplyr::select(genes, dataType, growthPhase, test_for) %>%
  dplyr::mutate(uncontrolled = 1)->uncontrolledGrowth

data %>%
  dplyr::filter(controlledForGrowth==1) %>%
  dplyr::select(genes, dataType, growthPhase, test_for) %>%
  dplyr::mutate(controlled = 1)->controlledGrowth
###*****************************


###*****************************
# join them
dplyr::full_join(uncontrolledGrowth, controlledGrowth) -> joined_data

joined_data$uncontrolled[is.na(joined_data$uncontrolled)]<-0
joined_data$controlled[is.na(joined_data$controlled)]<-0

joined_data %>%
  dplyr::mutate(exist_in = ifelse(test = (controlled==1 & uncontrolled==1) ,yes = "both", NA),
                exist_in = ifelse(test = (controlled==1 & uncontrolled==0) ,yes = "controlled for Growth", exist_in),
                exist_in = ifelse(test = (controlled==0 & uncontrolled==1) ,yes = "not controlled for Growth", exist_in))->joined_data 

joined_data$test_for <- factor(joined_data$test_for, levels=c("Low Mg", "High Mg", "High Na", "Gly", "Glc", "Lac"))
joined_data$exist_in <- factor(joined_data$exist_in, levels=c("controlled for Growth", "both", "not controlled for Growth"))
###*****************************


###*****************************
# Figure
fig01<-ggplot(joined_data, aes(x = test_for)) +
  geom_bar(aes(fill = exist_in), position = "fill") +
  facet_grid(growthPhase ~ dataType) +
  #scale_y_continuous(expand = c(0, 0),limits = c(0,2000))+
  theme_bw()+
  xlab("Test Condition") + 
  ylab("Number of differentially expressed genes") +
  theme(panel.margin.y = unit(2, "lines"),
        panel.grid.minor.x = element_blank(),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14),
        legend.background = element_rect(fill=alpha(colour = "white",
                                                    alpha = .4)))

print(fig01)
###*****************************


###*****************************
# save figure
cowplot::save_plot(filename = paste0("../c_figures/difference_rtw_GrowthControl",".pdf"),
                   plot = fig01, 
                   ncol = 2, nrow = 2, 
                   base_height = 3.075,base_width = 5.5,
                   units = "in",useDingbats=FALSE)
###*****************************