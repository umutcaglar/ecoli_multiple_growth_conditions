# generate figures for flux data
# WRITTEN BY: Viswanadham.Sridhara

# ---- initials ----

###*****************************
# INITIAL COMMANDS TO RESET THE SYSTEM
rm(list = ls())
if (is.integer(dev.list())){dev.off()}
cat("\014")
set.seed(14159)
###*****************************


###*****************************
# Set Working Directory
# One needs to arrange the correct pathway if this is not umut's computer ;)
if(as.vector(Sys.info()["effective_user"]=="umut"))
{setwd(paste0("/Users/umut/GitHub/ecoli_multiple_growth_conditions/e_Flux_Analyze/"))} # mac computer
###*****************************


###*****************************
#Load Libraries
library(reader) # for readLines()
library(dplyr)  # for mutate()
library(tidyr)  # for unnest()
library(purrr)  # for map(), reduce()
library(ggplot2) # for ggplot()
library(cowplot) # for plot_grid() making publication ready figures
library(stringr) # for str_extract()
library(broom) # for doing lm fit in parallel to calculate slopes
###*****************************


###*****************************
# Source functions
source("replaceFunction.R")
###*****************************


###*****************************
# Load data

#hard-code 14 lines for outputs and 36th line for FLUXRATIOS for now.
data_path="../Flux_reads_06_21_2016/"
files <- dir(data_path, pattern = "*ffr.TXT")


data <- data_frame(Filename = files) %>% # create a data frame
  # holding the file names
  mutate(file_contents = map(Filename,          # read files into
                             ~ n.readLines(file.path(data_path, .), 
                                           n=14, skip=36)), 
         Run = as.numeric(str_extract(str_extract(files,"(_)\\d+(_)"), "\\d+")), 
         sampleID = as.numeric(str_extract(str_extract(files,"(sample)\\d+"),"\\d+"))) # a new data column

data_unnest<-unnest(data)
###*****************************


###*****************************
# Add more information about data from From AG3C wiki pages
data_unnest$Salt<-"NULL"
data_unnest$Salt[which(data_unnest$sampleID>=106 & data_unnest$sampleID<=141)]<-"MGSO4"
data_unnest$Salt[which(data_unnest$sampleID>=16 & data_unnest$sampleID<=21)]<-"OTHER"
data_unnest$Salt[which(data_unnest$sampleID>=61 & data_unnest$sampleID<=84)]<-"NACL"

data_unnest$Conc<-NA
data_unnest$Conc[(data_unnest$sampleID) %in% seq(61,84,4)]<-5
data_unnest$Conc[(data_unnest$sampleID) %in% seq(62,84,4)]<-100
data_unnest$Conc[(data_unnest$sampleID) %in% seq(63,84,4)]<-200
data_unnest$Conc[(data_unnest$sampleID) %in% seq(64,84,4)]<-300
data_unnest$Conc[(data_unnest$sampleID) %in% seq(106,140,6)]<-0.08
data_unnest$Conc[(data_unnest$sampleID) %in% seq(107,140,6)]<-0.8
data_unnest$Conc[(data_unnest$sampleID) %in% seq(108,140,6)]<-8
data_unnest$Conc[(data_unnest$sampleID) %in% seq(109,127,6)]<-50
data_unnest$Conc[(data_unnest$sampleID) %in% seq(139,140,6)]<-50
data_unnest$Conc[(data_unnest$sampleID) %in% seq(110,128,6)]<-200
data_unnest$Conc[(data_unnest$sampleID) %in% seq(140,140,6)]<-200
data_unnest$Conc[(data_unnest$sampleID) %in% seq(111,129,6)]<-400

data_unnest$Phase<-NA
data_unnest$Phase[(data_unnest$sampleID) %in% c(seq(61,64,1),seq(69,72,1),seq(77,80,1),
                                                seq(106,111,1),seq(118,123,1),seq(130,132,1))]<-"EXP"
data_unnest$Phase[(data_unnest$sampleID) %in% c(seq(65,68,1),seq(73,76,1),seq(81,84,1),
                                                seq(112,117,1),seq(124,129,1),seq(136,140,1))]<-"STA"

data_unnest$Replicate<-NA
data_unnest$Replicate[(data_unnest$sampleID) %in% seq(61,68,1)]<-1
data_unnest$Replicate[(data_unnest$sampleID) %in% seq(69,76,1)]<-2
data_unnest$Replicate[(data_unnest$sampleID) %in% seq(77,84,1)]<-3
data_unnest$Replicate[(data_unnest$sampleID) %in% seq(106,117,1)]<-1
data_unnest$Replicate[(data_unnest$sampleID) %in% seq(118,129,1)]<-2
#samples 133,134,135,141 are not good
data_unnest$Replicate[(data_unnest$sampleID) %in% seq(130,132,1)]<-3
data_unnest$Replicate[(data_unnest$sampleID) %in% seq(136,140,1)]<-3
###*****************************


###*****************************
# Extracting flux ratios and branch information from ffr files
data_unnest$file_contents<-gsub("Labeled CO2", "Labeled CO2 EC", data_unnest$file_contents)
data_unnest$file_contents<-gsub("PYR from MAL \\(ub\\)", "PYR from MAL_UB", data_unnest$file_contents)
data_unnest$file_contents<-gsub("PYR from MAL \\(lb\\)", "PYR from MAL_LB", data_unnest$file_contents)
A<-str_extract(data_unnest$file_contents, "(\\d+\\.\\d+\\s+)")
A_<-str_extract(data_unnest$file_contents, "(-\\d+\\.\\d+\\s+)")
A[which(A_<0)]<-as.numeric(A[which(A_<0)])*-1
C<-str_extract(data_unnest$file_contents, "(\\S+\\ \\S+\\ \\S+)")
#names(A)<-C
data_unnest$Branch<-C
data_unnest$FluxRatio<-as.numeric(A)
#close(con)
data_unnest<-data_unnest[-which(data_unnest$Replicate %in% NA),]
#data_unnest<-data_unnest[,order(data_unnest$sampleID)]
###*****************************


###*****************************
# Averaging over the runs
data_unnest_aggregate<-aggregate(FluxRatio ~ sampleID + Branch + Conc + FluxRatio + Replicate + Salt + Phase, FUN = mean, data=data_unnest)
# values that fall above 1 are replaced with 1 and those that below 1e-5 are floored at 1e-5.
data_unnest_aggregate$FluxRatio[which(data_unnest_aggregate$FluxRatio<=1e-5)]=1e-5
data_unnest_aggregate$FluxRatio[which(data_unnest_aggregate$FluxRatio>=1)]=1
###*****************************


###*****************************
# Figure generation

# calculate Mean and std error of mean
data_unnest_Mean<-aggregate(FluxRatio ~  Phase + Salt + Conc + Branch, FUN = mean, data=data_unnest_aggregate)
sde <- function(x) sd(x)/sqrt(length(x))
data_unnest_SDE<-aggregate(FluxRatio ~  Phase + Salt + Conc + Branch, FUN = sde, data=data_unnest_aggregate)
colnames(data_unnest_Mean)[5]<-"MeanFluxRatio"
colnames(data_unnest_SDE)[5]<-"SDEFluxRatio"

# data frame that has the mean, sde for plotting
data_plots<-merge(data_unnest_Mean,data_unnest_SDE, 
                  by=intersect(names(data_unnest_SDE), 
                               names(data_unnest_Mean)))


# subsetting the data to simplify
dt1_NACL_EXP<-subset(data_plots, Salt=="NACL" & Phase=="EXP" )

dt1_NACL_STA<-subset(data_plots, Salt=="NACL" & Phase=="STA" )

dt1_MGSO4_EXP<-subset(data_plots, Salt=="MGSO4" & Phase=="EXP" )

dt1_MGSO4_STA<-subset(data_plots, Salt=="MGSO4" & Phase=="STA" )

# plotting
dt1NACLEXP<-ggplot(dt1_NACL_EXP,aes(x=Conc, y=MeanFluxRatio, colour=Branch)) + 
  geom_point(size=1.5) + 
  labs(y=" Flux Ratio") + 
  facet_grid(. ~Branch)+
  geom_errorbar(aes(y = MeanFluxRatio, 
                    ymin = MeanFluxRatio-SDEFluxRatio, 
                    ymax = MeanFluxRatio+SDEFluxRatio, 
                    color = Branch), 
                width = 5, 
                data = dt1_NACL_EXP)

dt1NACLSTA<-ggplot(dt1_NACL_STA,aes(x=Conc, y=MeanFluxRatio, colour=Branch)) + 
  geom_point(size=1.5) + 
  labs(y=" Flux Ratio") + 
  facet_grid(. ~Branch)+
  geom_errorbar(aes(y = MeanFluxRatio, 
                    ymin = MeanFluxRatio-SDEFluxRatio, 
                    ymax = MeanFluxRatio+SDEFluxRatio, 
                    color = Branch), 
                width = 5, 
                data = dt1_NACL_STA)

dt1MGSO4EXP<-ggplot(dt1_MGSO4_EXP,aes(x=Conc, y=MeanFluxRatio, colour=Branch)) + 
  geom_point(size=1.5) + 
  labs(y=" Flux Ratio") + 
  facet_grid(. ~Branch)+
  geom_errorbar(aes(y = MeanFluxRatio, 
                    ymin = MeanFluxRatio-SDEFluxRatio, 
                    ymax = MeanFluxRatio+SDEFluxRatio, 
                    color = Branch), 
                width = 5, 
                data = dt1_MGSO4_EXP)+
  scale_x_log10(limits=c(1e-2,400),
                expand = c(0, 0), 
                breaks=c(0.1,1,10,400),
                labels=c(0.1,1,10,400))

dt1MGSO4STA<-ggplot(dt1_MGSO4_STA,aes(x=Conc, y=MeanFluxRatio, colour=Branch)) + 
  geom_point(size=1.5) + 
  labs(y=" Flux Ratio") + 
  facet_grid(. ~Branch) +
  geom_errorbar(aes(y = MeanFluxRatio, 
                    ymin = MeanFluxRatio-SDEFluxRatio, 
                    ymax = MeanFluxRatio+SDEFluxRatio, 
                    color = Branch), 
                width = 5, 
                data = dt1_MGSO4_STA)+
  scale_x_log10(limits=c(1e-2,400),
                expand = c(0, 0), 
                breaks=c(0.1,1,10,400),
                labels=c(0.1,1,10,400))


dt<-plot_grid(dt1NACLEXP,dt1NACLSTA,dt1MGSO4EXP,dt1MGSO4STA, align="h",nrow=2)
print(dt)

ggsave("../e_figures/SaltStressPlots.pdf", plot=last_plot(),width = 48, height = 24)
###*****************************


###*****************************
# Figure Generation (My Turn)
data_plots %>% dplyr::mutate(Phase2=replace_fun(input_vector = Phase,
                                                initialVal = c("EXP", "STA"),
                                                finalVal = c("Exponential", "Stationary")))->data_plots


write.csv(x = data_plots, file = "../e_results/flux_data.csv")

# data_plots %>% dplyr::filter(Phase=="EXP")->data_plots_exp
# data_plots %>% dplyr::filter(Phase=="STA")->data_plots_sta
# 
# data_plots %>%dplyr::filter(Salt=="MGSO4") ->data_plots_Mg
# data_plots %>%dplyr::filter(Salt=="NACL")->data_plots_Na
###*****************************


###*****************************
# Best line calculation
data_plots %>%
  dplyr::mutate(Conc_transf = ifelse(Salt=='MGSO4', log(Conc), Conc)) %>%
  dplyr::group_by(Salt, Phase, Branch) %>%
  dplyr::do(broom::tidy(lm(MeanFluxRatio ~ Conc_transf, data=.))) %>%
  dplyr::filter(term == 'Conc_transf') %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Phase) %>%
  dplyr::mutate(p.adjusted = p.adjust(p.value, method="fdr")) -> data_plot_fits

write.csv(x = data_plot_fits,file = "../e_results/flux_p_values.csv")

###*****************************


###*****************************
# Font size
facetFontSize=12
###*****************************


###*****************************
data_plots %>%dplyr::filter(Salt=="MGSO4", Phase=="EXP") ->data_plots_Mg_exp
data_plots %>%dplyr::filter(Salt=="NACL", Phase=="EXP")->data_plots_Na_exp
data_plots %>%dplyr::filter(Salt=="MGSO4", Phase=="STA") ->data_plots_Mg_sta
data_plots %>%dplyr::filter(Salt=="NACL", Phase=="STA")->data_plots_Na_sta


fig_Na_exp<-ggplot2::ggplot(data_plots_Na_exp,aes(x=Conc, y=MeanFluxRatio, colour=Branch)) + 
  geom_point(size=1.5) +
  geom_line() +
  facet_wrap( ~ Branch,ncol = 7)+
  geom_errorbar(aes(y = MeanFluxRatio, 
                    ymin = MeanFluxRatio-SDEFluxRatio, 
                    ymax = MeanFluxRatio+SDEFluxRatio, 
                    color = Branch), 
                width = 5, 
                data = data_plots_Na_exp)+
  labs(y="Flux Ratio", x="Concentration") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text.x = element_text(size = facetFontSize, face = "bold"),
        strip.text.y = element_text(size = facetFontSize,face = "bold"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.position="none",
        panel.grid.minor = element_blank())

print(fig_Na_exp)


fig_Na_sta<-ggplot2::ggplot(data_plots_Na_sta,aes(x=Conc, y=MeanFluxRatio, colour=Branch)) + 
  geom_point(size=1.5) +
  geom_line() +
  facet_wrap( ~ Branch,ncol = 7)+
  geom_errorbar(aes(y = MeanFluxRatio, 
                    ymin = MeanFluxRatio-SDEFluxRatio, 
                    ymax = MeanFluxRatio+SDEFluxRatio, 
                    color = Branch), 
                width = 5, 
                data = data_plots_Na_sta)+
  labs(y="Flux Ratio", x="Concentration") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text.x = element_text(size = facetFontSize, face = "bold"),
        strip.text.y = element_text(size = facetFontSize, face = "bold"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.position="none",
        panel.grid.minor = element_blank())

print(fig_Na_sta)

fig_Mg_exp<-ggplot2::ggplot(data_plots_Mg_exp,aes(x=Conc, y=MeanFluxRatio, colour=Branch)) + 
  geom_point(size=1.5) +
  geom_line() +
  scale_x_log10(limits=c(1e-2,400),
                breaks=c(0.1,1,10,400),
                labels=c(0.1,1,10,400))+
  facet_wrap( ~ Branch,ncol = 7)+
  geom_errorbar(aes(y = MeanFluxRatio, 
                    ymin = MeanFluxRatio-SDEFluxRatio, 
                    ymax = MeanFluxRatio+SDEFluxRatio, 
                    color = Branch), 
                width = 5, 
                data = data_plots_Mg_exp)+
  labs(y="Flux Ratio", x="Concentration") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text.x = element_text(size = facetFontSize, face = "bold"),
        strip.text.y = element_text(size = facetFontSize, face = "bold"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.position="none",
        panel.grid.minor = element_blank())

print(fig_Mg_exp)

fig_Mg_sta<-ggplot2::ggplot(data_plots_Mg_sta,aes(x=Conc, y=MeanFluxRatio, colour=Branch)) + 
  geom_point(size=1.5) +
  geom_line() +
  scale_x_log10(limits=c(1e-2,400),
                breaks=c(0.1,1,10,400),
                labels=c(0.1,1,10,400))+
  facet_wrap( ~ Branch,ncol = 7)+
  geom_errorbar(aes(y = MeanFluxRatio, 
                    ymin = MeanFluxRatio-SDEFluxRatio, 
                    ymax = MeanFluxRatio+SDEFluxRatio, 
                    color = Branch), 
                width = 5, 
                data = data_plots_Mg_sta)+
  labs(y="Flux Ratio", x="Concentration") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text.x = element_text(size = facetFontSize, face = "bold"),
        strip.text.y = element_text(size = facetFontSize, face = "bold"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.position="none",
        panel.grid.minor = element_blank())

print(fig_Mg_sta)


# Combine plots 
fig_exp<-cowplot::plot_grid(fig_Na_exp, fig_Mg_exp, labels = c("A", "B"), nrow = 2, scale = .95)
print(fig_exp)

fig_sta<-cowplot::plot_grid(fig_Na_sta, fig_Mg_sta, labels = c("A", "B"), nrow = 2, scale = .95)
print(fig_sta)


cowplot::save_plot(filename = "../e_figures/Exp.pdf",plot = fig_exp,nrow = 1.7*2,ncol = 3)
cowplot::save_plot(filename = "../e_figures/Sta.pdf",plot = fig_sta,nrow = 1.7*2,ncol = 3)
###*****************************