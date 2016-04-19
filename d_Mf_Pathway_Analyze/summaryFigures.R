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
              "d_Mf_Pathway_Analyze/"))} # mac computer
###*****************************


###*****************************
# DOWNLOAD LIBRARIES
require("dplyr")
require("tidyr")
require("ggplot2")
require("RColorBrewer")
require("scales")
require("cowplot")
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


fileList=dir("../d_results/",pattern = "*.mf")

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

metaData=read.csv(paste0(file="../c_results/",
                         "metaData_mrna_trT_set02_StcNasMgh_S_baseMghighMg_baseNa",
                         "_Sta_noFilter_p1Sf_noNorm__Mg_mM_Levels.csv"))

referanceMF=read.table(file = "ReferenceFiles/MF_david_mrna_ref_2009_tidy.txt",header = TRUE,sep="\t",quote = "")
###*****************************


###*****************************
# Generate flegellar data frame kegg
mainLoadedFile_kegg$vs=as.character(mainLoadedFile_kegg$vs)
mainLoadedFile_kegg %>%
  dplyr::group_by()%>%
  dplyr::mutate(pick_data=ifelse(pick_data=="mrna","mRNA","Protein"))%>%
  dplyr::mutate(vs=ifelse(vs=="baseNahighNa","highNa",vs),
                vs=ifelse(vs=="baseMghighMg","highMg",vs),
                vs=ifelse(vs=="baseMglowMg","lowMg",vs),
                vs=ifelse(vs=="glucoselactate","lactate",vs),
                vs=ifelse(vs=="glucoseglycerol","glycerol",vs),
                vs=ifelse(vs=="glucosegluconate","gluconate",vs))%>%
  dplyr::group_by(KEGG_Path, padj_gene, pick_data,growthPhase,vs) %>%
  dplyr::mutate(grouping=paste(growthPhase,
                               vs,sep = "-"))->kegg_metabolism_df


kegg_metabolism_df %>%
  dplyr::group_by(KEGG_Path,FDR_KEGG_Path, KEGG_Path_long, KEGG_Path_short) %>%
  dplyr::summarize(numGenesInCat=unique(numGenesInCat),
                   numSigP=unique(numSigP),
                   numSigN=unique(numSigN),
                   pick_data=unique(pick_data),
                   growthPhase=unique(growthPhase),
                   test_for=unique(test_for),
                   vs=unique(vs),
                   df_category=unique(df_category),
                   grouping=unique(grouping))%>%
  dplyr::group_by(pick_data, growthPhase, vs)%>%
  dplyr::arrange(FDR_KEGG_Path)%>%
  dplyr::mutate(rank=seq(1:n()))%>%
  dplyr::filter(rank<6) %>%
  tidyr::complete(rank)->kegg_metabolism_summary


kegg_metabolism_summary$grouping <- factor(kegg_metabolism_summary$grouping, 
                                           levels = (c("Exp-lowMg", "Exp-highMg", "Exp-highNa",
                                                       "Exp-glycerol","Exp-gluconate","Exp-lactate",
                                                       "Sta-lowMg", "Sta-highMg", "Sta-highNa",
                                                       "Sta-glycerol","Sta-gluconate","Sta-lactate")))

kegg_metabolism_summary$vs <- factor(kegg_metabolism_summary$vs, 
                                     levels = (c("lowMg", "highMg", "highNa",
                                                 "glycerol","gluconate","lactate")))
###*****************************


###*****************************
# MF Vector

# Generate flegellar data frame mf
mainLoadedFile_mf$vs=as.character(mainLoadedFile_mf$vs)
mainLoadedFile_mf %>%
  dplyr::group_by()%>%
  dplyr::mutate(pick_data=ifelse(pick_data=="mrna","mRNA","Protein"))%>%
  dplyr::mutate(vs=ifelse(vs=="baseNahighNa","highNa",vs),
                vs=ifelse(vs=="baseMghighMg","highMg",vs),
                vs=ifelse(vs=="baseMglowMg","lowMg",vs),
                vs=ifelse(vs=="glucoselactate","lactate",vs),
                vs=ifelse(vs=="glucoseglycerol","glycerol",vs),
                vs=ifelse(vs=="glucosegluconate","gluconate",vs))%>%
  dplyr::group_by(MF_Name, padj_gene, pick_data,growthPhase,vs) %>%
  dplyr::mutate(grouping=paste(growthPhase,
                               vs,sep = "-"))->mf_metabolism_df

mf_metabolism_df %>%
  dplyr::group_by(MF_Name, FDR_GoMF, MF_Name_long, MF_Name_short) %>%
  dplyr::summarize(numGenesInCat=unique(numGenesInCat),
                   numSigP=unique(numSigP),
                   numSigN=unique(numSigN),
                   pick_data=unique(pick_data),
                   growthPhase=unique(growthPhase),
                   test_for=unique(test_for),
                   vs=unique(vs),
                   df_category=unique(df_category),
                   grouping=unique(grouping))%>%
  dplyr::group_by(pick_data, growthPhase, vs)%>%
  dplyr::arrange(FDR_GoMF)%>%
  dplyr::mutate(rank=seq(1:n()))%>%
  dplyr::filter(rank<6) %>%
  tidyr::complete(rank)->mf_metabolism_summary

mf_metabolism_summary$grouping <- factor(mf_metabolism_summary$grouping, 
                                         levels = rev(c("Exp-lowMg", "Exp-highMg", "Exp-highNa",
                                                        "Exp-glycerol","Exp-gluconate","Exp-lactate",
                                                        "Sta-lowMg", "Sta-highMg", "Sta-highNa",
                                                        "Sta-glycerol","Sta-gluconate","Sta-lactate")))

mf_metabolism_summary$vs <- factor(mf_metabolism_summary$vs, 
                                   levels = (c("lowMg", "highMg", "highNa",
                                               "glycerol","gluconate","lactate")))
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
  geom_text(aes(label=paste(rank, KEGG_Path_short, sep=". "), x=0.02),size=3, hjust=0)+
  facet_grid(growthPhase+vs~pick_data) +
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
  geom_text(aes(label=paste(rank, MF_Name_short, sep=". "), x=0.02),size=3, hjust=0)+
  facet_grid(growthPhase+vs~pick_data)+
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
  geom_text(aes(label=paste(rank, KEGG_Path_short, sep=". "), x=0.02),size=3, hjust=0)+
  facet_grid(vs~pick_data) +
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
  geom_text(aes(label=paste(rank, KEGG_Path_short, sep=". "), x=0.02),size=3, hjust=0)+
  facet_grid(vs~pick_data) +
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
  geom_text(aes(label=paste(rank, MF_Name_short, sep=". "), x=0.02),size=3, hjust=0)+
  facet_grid(vs~pick_data) +
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
  geom_text(aes(label=paste(rank, MF_Name_short, sep=". "), x=0.02),size=3, hjust=0)+
  facet_grid(vs~pick_data) +
  panel_border() +
  theme(  axis.line = element_blank(), 
          axis.text.x = element_blank(), 
          axis.text.y = element_blank(),
          axis.ticks = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank())
###*****************************


###*****************************
# Save Figures
cowplot::save_plot(filename = paste0("../d_figures/resultTable_kegg.pdf"),
                   plot = fig01,ncol = 2,nrow=4,limitsize = FALSE)

cowplot::save_plot(filename = paste0("../d_figures/resultTable_mf.pdf"),
                   plot = fig02,ncol = 2,nrow=4,limitsize = FALSE)

cowplot::save_plot(filename = paste0("../d_figures/resultTable_kegg_exp.pdf"),
                   plot = fig01a,ncol = 2,nrow=1.7,limitsize = FALSE)
cowplot::save_plot(filename = paste0("../d_figures/resultTable_kegg_sta.pdf"),
                   plot = fig01b,ncol = 2,nrow= 1.7,limitsize = FALSE)
cowplot::save_plot(filename = paste0("../d_figures/resultTable_mf_exp.pdf"),
                   plot = fig02a,ncol = 2,nrow=1.7,limitsize = FALSE)
cowplot::save_plot(filename = paste0("../d_figures/resultTable_mf_sta.pdf"),
                   plot = fig02b,ncol = 2.05,nrow=1.7,limitsize = FALSE)