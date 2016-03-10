# Aim of the code is to generate figure for binding related pathways


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
###*****************************


###*****************************
# Generate flegellar data frame kegg
mainLoadedFile_kegg$vs=as.character(mainLoadedFile_kegg$vs)
mainLoadedFile_kegg %>%
  dplyr::filter(grepl("*.*binding*.*",mainLoadedFile_kegg$KEGG_Path_short)) %>%
  dplyr::mutate(pick_data=ifelse(pick_data=="mrna","mrna","protein"))%>%
  dplyr::mutate(vs=ifelse(vs=="baseNahighNa","highNa",vs),
                vs=ifelse(vs=="baseMghighMg","highMg",vs),
                vs=ifelse(vs=="baseMglowMg","lowMg",vs))%>%
  dplyr::group_by(KEGG_Path, padj_gene, pick_data,growthPhase,vs) %>%
  dplyr::filter(pick_data != "protein")%>%
  dplyr::mutate(grouping=paste(vs,sep = "-"))->kegg_binding_df

###*****************************


###*****************************
# Generate binding data frame mf
mainLoadedFile_mf$vs=as.character(mainLoadedFile_mf$vs)
mainLoadedFile_mf %>%
  dplyr::filter(grepl("*.*binding*.*",MF_Name_short)) %>%
  dplyr::mutate(pick_data=ifelse(pick_data=="mrna","mrna","protein"))%>%
  dplyr::mutate(vs=ifelse(vs=="baseNahighNa","highNa",vs),
                vs=ifelse(vs=="baseMghighMg","highMg",vs),
                vs=ifelse(vs=="baseMglowMg","lowMg",vs))%>%
  dplyr::group_by(MF_Name, padj_gene, pick_data,growthPhase,vs) %>%
  dplyr::mutate(grouping=paste(pick_data,
                               growthPhase,
                               vs,sep = "-"))->mf_binding_df

mf_binding_df %>%
  dplyr::group_by(pick_data, growthPhase, vs, grouping, signChange) %>%
  dplyr::summarize(lengthObj=length(score_gene))%>%
  dplyr::group_by() %>%
  tidyr::complete(grouping,signChange)->mf_binding_summary

mf_binding_df %>% 
  dplyr::filter(MF_Name_short %in% c("magnesium ion binding",
                                     "iron ion binding",
                                     "zinc ion binding",
                                     "manganese ion binding",
                                     "cobalt ion binding",
                                     "calcium ion binding"))%>%
  dplyr::filter(test_for !="carbonSource",
                pick_data=="mrna")%>%
  dplyr::group_by(MF_Name, padj_gene, pick_data,growthPhase,vs) %>%
  dplyr::mutate(grouping=paste(vs,
                               sep = "-"))->mf_binding_ions_df

mf_binding_ions_df$MF_Name_short=factor((mf_binding_ions_df$MF_Name_short))

mf_binding_ions_df %>%
  dplyr::group_by(MF_Name_short, pick_data, growthPhase, vs, grouping, signChange) %>%
  dplyr::summarize(lengthObj=length(score_gene))%>%
  dplyr::group_by() %>%
  tidyr::complete(grouping,signChange,MF_Name_short,pick_data)->mf_binding_ions_summary

mf_binding_ions_summary$MF_Name_short=factor((mf_binding_ions_summary$MF_Name_short))

mf_binding_ions_summary$grouping <- factor(mf_binding_ions_summary$grouping, 
                                           levels = rev(c("lowMg", "highMg", "highNa")))


###*****************************


###*****************************
# Generate binding figure mf
fig01=ggplot(mf_binding_summary, aes(x=grouping,
                                     y=lengthObj,
                                     fill=as.factor(signChange))) +
  geom_bar(position="dodge", stat="identity",width=.75)+
  scale_fill_manual(values = c("blue","red"),
                    name="Regulation",
                    breaks=c("-1", "1"),
                    labels=c("Down Regulated", "Up Regulated"))+
  xlab("conditions")+
  ylab("number of differentially expressed")+
  theme_bw()+
  coord_flip()

print(fig01)

fig02=ggplot(mf_binding_ions_summary, aes(x=MF_Name_short,
                                          y=lengthObj,
                                          fill=as.factor(signChange))) +
  #facet_grid( ~ MF_Name_short)+
  geom_bar(position="dodge", stat="identity",width=.5)+
  scale_fill_manual(values = c("blue","red"),
                    name="Regulation",
                    breaks=c("-1", "1"),
                    labels=c("Down Regulated", "Up Regulated"))+
  xlab("conditions")+
  ylab("number of differentially expressed")+
  theme_bw()

print(fig02)

cowplot::save_plot(filename = paste0("../d_figures/binding_mf.pdf"),
                   plot = fig02,ncol = 2.5,limitsize = FALSE)
