# Aim of the code is to generate figure for Metabolites

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
              "d_code_Mf_Pathway_Analyze/"))} # mac computer
###*****************************


###*****************************
# DOWNLOAD LIBRARIES
require("dplyr")
require("tidyr")
require("ggplot2")
require("RColorBrewer")
require("scales")
require("cowplot")
require("Cairo")
###*****************************


###*****************************
#Load Functions
source("../b_code_histogram_RNA&Protein/replace_fun.R")
###*****************************


###*****************************
# Load files
fileList=dir("../d_results/",pattern = "*.kegg")
for(counter01 in 1 : length(fileList))
{
  if (counter01 ==1)
  {
    mainLoadedFile_kegg=read.csv(file = paste0("../d_results/",fileList[counter01]))
  }
  if (counter01 !=1)
  {
    temp=read.csv(file = paste0("../d_results/",fileList[counter01]),header = TRUE)
    mainLoadedFile_kegg=rbind(mainLoadedFile_kegg,temp)
  }
}


fileList=dir("../d_results/",pattern = "*.mf_n")

for(counter01 in 1 : length(fileList))
{
  if (counter01 ==1)
  {
    mainLoadedFile_mf=read.csv(file = paste0("../d_results/",fileList[counter01]))
  }
  if (counter01 !=1)
  {
    temp=read.csv(file = paste0("../d_results/",fileList[counter01]),header = TRUE)
    mainLoadedFile_mf=rbind(mainLoadedFile_mf,temp)
  }
}
###*****************************


###*****************************
# Generate flegellar data frame kegg
mainLoadedFile_kegg %>%
  dplyr::group_by()%>%
  dplyr::mutate(pick_data=ifelse(pick_data=="mrna","mRNA","Protein"))%>%
  dplyr::group_by(KEGG_Path, padj_gene, pick_data,growthPhase,contrast) %>%
  dplyr::mutate(grouping=paste(growthPhase,
                               contrast,sep = "-"))->kegg_metabolism_df


kegg_metabolism_df %>%
  dplyr::group_by(KEGG_Path,FDR_KEGG_Path, KEGG_Path_long, KEGG_Path_short) %>%
  dplyr::summarize(Pop.Hits=unique(Pop.Hits),
                   numSigP=unique(numSigP),
                   numSigN=unique(numSigN),
                   pick_data=unique(pick_data),
                   growthPhase=unique(growthPhase),
                   test_for=unique(test_for),
                   contrast=unique(contrast),
                   df_category=unique(df_category),
                   grouping=unique(grouping))%>%
  dplyr::group_by(pick_data, growthPhase, contrast)%>%
  dplyr::arrange(FDR_KEGG_Path)%>%
  dplyr::mutate(rank=seq(1:n()))%>%
  #dplyr::filter(rank<6) %>%
  tidyr::complete(rank) %>%
  dplyr::mutate(PlusMinusRatio=numSigP/numSigN) %>%
  dplyr::mutate(PlusMinusSign=ifelse(PlusMinusRatio>=95/5,"\u25B2\u25B2\u25B2","none"),
                PlusMinusSign=ifelse(PlusMinusRatio>=80/20 & PlusMinusRatio < 95/5, "\u25B2\u25B2",PlusMinusSign),
                PlusMinusSign=ifelse(PlusMinusRatio>=60/40 & PlusMinusRatio < 80/20, "\u25B2",PlusMinusSign),
                PlusMinusSign=ifelse(PlusMinusRatio>=40/60 & PlusMinusRatio < 60/40, "\u25B2 \u25BC",PlusMinusSign),
                PlusMinusSign=ifelse(PlusMinusRatio>=20/80 & PlusMinusRatio < 40/60, "\u25BC",PlusMinusSign),
                PlusMinusSign=ifelse(PlusMinusRatio>=5/95 & PlusMinusRatio < 20/80, "\u25BC\u25BC",PlusMinusSign),
                PlusMinusSign=ifelse(PlusMinusRatio < 5/95, "\u25BC\u25BC\u25BC",PlusMinusSign))->kegg_metabolism_summary


kegg_metabolism_summary$grouping <- factor(kegg_metabolism_summary$grouping, 
                                           levels = (c("Exp-lowMg", "Exp-highMg", "Exp-highNa",
                                                       "Exp-glycerol","Exp-gluconate","Exp-lactate",
                                                       "Sta-lowMg", "Sta-highMg", "Sta-highNa",
                                                       "Sta-glycerol","Sta-gluconate","Sta-lactate")))

kegg_metabolism_summary$contrast <- factor(kegg_metabolism_summary$contrast, 
                                     levels = (c("lowMg", "highMg", "highNa",
                                                 "glycerol","gluconate","lactate")))
###*****************************


###*****************************
# MF Vector

# Generate flegellar data frame mf
mainLoadedFile_mf %>%
  dplyr::group_by()%>%
  dplyr::mutate(pick_data=ifelse(pick_data=="mrna","mRNA","Protein"))%>%
  dplyr::group_by(MF, padj_gene, pick_data,growthPhase,contrast) %>%
  dplyr::mutate(grouping=paste(growthPhase,
                               contrast,sep = "-"))->mf_metabolism_df

mf_metabolism_df %>%
  dplyr::group_by(MF, FDR_MF, MF_Long, MF_Short) %>%
  dplyr::summarize(Pop.Hits=unique(Pop.Hits),
                   numSigP=unique(numSigP),
                   numSigN=unique(numSigN),
                   pick_data=unique(pick_data),
                   growthPhase=unique(growthPhase),
                   test_for=unique(test_for),
                   contrast=unique(contrast),
                   df_category=unique(df_category),
                   grouping=unique(grouping))%>%
  dplyr::group_by(pick_data, growthPhase, contrast)%>%
  dplyr::arrange(FDR_MF)%>%
  dplyr::mutate(rank=seq(1:n()))%>%
  #dplyr::filter(rank<6) %>%
  tidyr::complete(rank)%>%
  dplyr::mutate(PlusMinusRatio=numSigP/numSigN)%>%
  dplyr::mutate(PlusMinusSign=ifelse(PlusMinusRatio>=95/5,"\u25B2\u25B2\u25B2","none"),
                PlusMinusSign=ifelse(PlusMinusRatio>=80/20 & PlusMinusRatio < 95/5, "\u25B2\u25B2",PlusMinusSign),
                PlusMinusSign=ifelse(PlusMinusRatio>=60/40 & PlusMinusRatio < 80/20, "\u25B2",PlusMinusSign),
                PlusMinusSign=ifelse(PlusMinusRatio>=40/60 & PlusMinusRatio < 60/40, "\u25B2 \u25BC",PlusMinusSign),
                PlusMinusSign=ifelse(PlusMinusRatio>=20/80 & PlusMinusRatio < 40/60, "\u25BC",PlusMinusSign),
                PlusMinusSign=ifelse(PlusMinusRatio>=5/95 & PlusMinusRatio < 20/80, "\u25BC\u25BC",PlusMinusSign),
                PlusMinusSign=ifelse(PlusMinusRatio < 5/95, "\u25BC\u25BC\u25BC",PlusMinusSign))->mf_metabolism_summary

mf_metabolism_summary$grouping <- factor(mf_metabolism_summary$grouping, 
                                         levels = rev(c("Exp-lowMg", "Exp-highMg", "Exp-highNa",
                                                        "Exp-glycerol","Exp-gluconate","Exp-lactate",
                                                        "Sta-lowMg", "Sta-highMg", "Sta-highNa",
                                                        "Sta-glycerol","Sta-gluconate","Sta-lactate")))

mf_metabolism_summary$contrast <- factor(mf_metabolism_summary$contrast, 
                                   levels = (c("lowMg", "highMg", "highNa",
                                               "glycerol","gluconate","lactate")))

mf_metabolism_summary$pick_data <- factor(mf_metabolism_summary$pick_data, 
                                         levels = (c("mRNA","Protein")))
###*****************************


###*****************************
# Divide data frame into two pieces
kegg_metabolism_summary%>%
  dplyr::filter(growthPhase=="Exp")->kegg_metabolism_summary_exp

kegg_metabolism_summary%>%
  dplyr::filter(growthPhase=="Sta")->kegg_metabolism_summary_sta

mf_metabolism_summary%>%
  dplyr::filter(growthPhase=="Exp")->mf_metabolism_summary_exp

mf_metabolism_summary%>%
  dplyr::filter(growthPhase=="Sta")->mf_metabolism_summary_sta
###*****************************


###*****************************
# Generate Tables
fig01<-ggplot(kegg_metabolism_summary, aes( y=rank,x="condition"))+
  ylim(5.5,0.5)+
  scale_x_continuous(limits=c(0,1), expand=c(0,0)) +
  geom_text(aes(label=paste0(rank, ".", KEGG_Path_short, "  ", PlusMinusSign), x=0.02),size=3, hjust=0)+
  facet_grid(growthPhase+contrast~pick_data,drop = FALSE) +
  panel_border() +
  theme(  axis.line = element_blank(), 
          axis.text.x = element_blank(), 
          axis.text.y = element_blank(),
          axis.ticks = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank())

print(fig01)


fig02<-ggplot(mf_metabolism_summary, aes( y=rank,x="condition"))+
  ylim(5.5,0.5)+
  scale_x_continuous(limits=c(0,1), expand=c(0,0)) +
  geom_text(aes(label=paste0(rank, ".", MF_Short, "  ", PlusMinusSign), x=0.02),size=3, hjust=0)+
  facet_grid(growthPhase+contrast~pick_data,drop = FALSE)+
  panel_border() +
  theme(  axis.line = element_blank(), 
          axis.text.x = element_blank(), 
          axis.text.y = element_blank(),
          axis.ticks = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank())

print(fig02)

fig01a<-ggplot(kegg_metabolism_summary_exp, aes( y=rank,x="condition"))+
  ylim(5.5,0.5)+
  scale_x_continuous(limits=c(0,1), expand=c(0,0)) +
  geom_text(aes(label=paste0(rank, ".", KEGG_Path_short, "  ", PlusMinusSign), x=0.02),size=3, hjust=0)+
  facet_grid(contrast~pick_data,drop = FALSE) +
  panel_border() +
  theme(  axis.line = element_blank(), 
          axis.text.x = element_blank(), 
          axis.text.y = element_blank(),
          axis.ticks = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank())

fig01b<-ggplot(kegg_metabolism_summary_sta, aes( y=rank,x="condition"))+
  ylim(5.5,0.5)+
  scale_x_continuous(limits=c(0,1), expand=c(0,0)) +
  geom_text(aes(label=paste0(rank, ".", KEGG_Path_short, "  ", PlusMinusSign), x=0.02),size=3, hjust=0)+
  facet_grid(contrast~pick_data,drop = FALSE) +
  panel_border() +
  theme(  axis.line = element_blank(), 
          axis.text.x = element_blank(), 
          axis.text.y = element_blank(),
          axis.ticks = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank())

fig02a<-ggplot(mf_metabolism_summary_exp, aes( y=rank,x="condition"))+
  ylim(5.5,0.5)+
  scale_x_continuous(limits=c(0,1), expand=c(0,0)) +
  geom_text(aes(label=paste0(rank, ".", MF_Short, "  ", PlusMinusSign), x=0.02),size=3, hjust=0)+
  facet_grid(contrast~pick_data, drop = FALSE) +
  panel_border() +
  theme(  axis.line = element_blank(), 
          axis.text.x = element_blank(), 
          axis.text.y = element_blank(),
          axis.ticks = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank())



fig02b<-ggplot(mf_metabolism_summary_sta, aes( y=rank,x="condition"))+
  ylim(5.5,0.5)+
  scale_x_continuous(limits=c(0,1), expand=c(0,0)) +
  geom_text(aes(label=paste0(rank, ".", MF_Short, "  ", PlusMinusSign), x=0.02),size=3, hjust=0)+
  facet_grid(contrast~pick_data,drop = FALSE) +
  panel_border() +
  theme(  axis.line = element_blank(), 
          axis.text.x = element_blank(), 
          axis.text.y = element_blank(),
          axis.ticks = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank())
###*****************************


###*****************************
# Combine Kegg figures together
kegg_figure<-cowplot::plot_grid(fig01a,fig01b,labels = c("A","B"),ncol=1,scale = .95)
mf_figure<-cowplot::plot_grid(fig02a,fig02b,labels = c("A","B"),ncol=1,scale = .95)
###*****************************


###*****************************
# Save Figures
cowplot::save_plot(filename = paste0("../d_figures/resultTable_kegg.pdf"),
                   plot = fig01,ncol = 2,nrow=4,limitsize = FALSE, device=cairo_pdf)

cowplot::save_plot(filename = paste0("../d_figures/resultTable_mf.pdf"),
                   plot = fig02,ncol = 2,nrow=4,limitsize = FALSE, device=cairo_pdf)

cowplot::save_plot(filename = paste0("../d_figures/resultTable_kegg_exp.pdf"),
                   plot = fig01a,ncol = 2,nrow=1.7,limitsize = FALSE, device=cairo_pdf)
cowplot::save_plot(filename = paste0("../d_figures/resultTable_kegg_sta.pdf"),
                   plot = fig01b,ncol = 2,nrow= 1.7,limitsize = FALSE, device=cairo_pdf)
cowplot::save_plot(filename = paste0("../d_figures/resultTable_mf_exp.pdf"),
                   plot = fig02a,ncol = 2,nrow=1.7,limitsize = FALSE, device=cairo_pdf)
cowplot::save_plot(filename = paste0("../d_figures/resultTable_mf_sta.pdf"),
                   plot = fig02b,ncol = 2.05,nrow=1.7,limitsize = FALSE, device=cairo_pdf)


cowplot::save_plot(filename = paste0("../d_figures/resultTable_kegg.pdf"),
                   plot = kegg_figure, ncol = 2,nrow=3.4,limitsize = FALSE, device=cairo_pdf)

cowplot::save_plot(filename = paste0("../d_figures/resultTable_kegg.png"),
                   plot = kegg_figure, ncol = 2,nrow=3.4,limitsize = FALSE, dpi = 600)

cowplot::save_plot(filename = paste0("../d_figures/resultTable_mf.pdf"),
                   plot = mf_figure, ncol = 2,nrow=3.4,limitsize = FALSE, device=cairo_pdf)

cowplot::save_plot(filename = paste0("../d_figures/resultTable_mf.png"),
                   plot = mf_figure, ncol = 2,nrow=3.4,limitsize = FALSE, dpi = 600)