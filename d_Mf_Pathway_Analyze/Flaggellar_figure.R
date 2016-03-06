# Aim of the code is to generate figure for Flaggeler and motion


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
  dplyr::filter(KEGG_Path_short=="Flagellar assembly") %>%
  dplyr::mutate(pick_data=ifelse(pick_data=="mrna","mrna","protein"))%>%
  dplyr::mutate(vs=ifelse(vs=="baseNahighNa","highNa",vs),
                vs=ifelse(vs=="baseMghighMg","highMg",vs),
                vs=ifelse(vs=="baseMglowMg","lowMg",vs))%>%
  dplyr::group_by(KEGG_Path, padj_gene, pick_data,growthPhase,vs) %>%
  dplyr::filter(pick_data != "protein")%>%
  dplyr::mutate(grouping=paste(vs,sep = "-"))->kegg_flagellar_assembly_df

kegg_flagellar_assembly_df$grouping <- factor(kegg_flagellar_assembly_df$grouping, 
                                            levels = rev(c("lowMg", "highMg", "highNa")))
###*****************************


###*****************************
# Generate flegellar data frame mf
mainLoadedFile_mf$vs=as.character(mainLoadedFile_mf$vs)
mainLoadedFile_mf %>%
  dplyr::filter(MF_Name_short=="motor activity") %>%
  dplyr::mutate(pick_data=ifelse(pick_data=="mrna","mrna","protein"))%>%
  dplyr::mutate(vs=ifelse(vs=="baseNahighNa","highNa",vs),
                vs=ifelse(vs=="baseMghighMg","highMg",vs),
                vs=ifelse(vs=="baseMglowMg","lowMg",vs))%>%
  dplyr::group_by(MF_Name, padj_gene, pick_data,growthPhase,vs) %>%
  dplyr::mutate(grouping=paste(pick_data,
                               growthPhase,
                               vs,sep = "-"))->mf_flagellar_assembly_df
###*****************************


###*****************************
# Generate filegellar figure kegg
scaleHigh=max(abs(kegg_flagellar_assembly_df$score_gene))
scaleMid=0
scaleLow=-max(abs(kegg_flagellar_assembly_df$score_gene))

fig01=ggplot(kegg_flagellar_assembly_df, aes( x=rank,y=grouping)) +
  geom_tile(aes(fill=score_gene))+
  scale_fill_gradientn(colours=c("Blue","Grey50","Red"),
                       values=rescale(c(scaleLow,scaleMid,scaleHigh)),
                       limits=c(scaleLow,scaleHigh),
                       guide = guide_colorbar(title = "-sign(cor)*P_log10"))+
  geom_text(aes(label=ID),size=3, colour="White", fontface="bold")+
  theme_bw()+
  scale_x_continuous(breaks=min(kegg_flagellar_assembly_df$rank):max(kegg_flagellar_assembly_df$rank))+
  theme(axis.line.y = element_blank(),
        legend.position="bottom",
        axis.title.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())

print(fig01)
###*****************************


###*****************************
# Generate filegellar figure mf
scaleHigh=max(abs(mf_flagellar_assembly_df$score_gene))
scaleMid=0
scaleLow=-max(abs(mf_flagellar_assembly_df$score_gene))

fig02=ggplot(mf_flagellar_assembly_df, aes( x=rank,y=grouping)) +
  geom_tile(aes(fill=score_gene))+
  scale_fill_gradientn(colours=c("Blue","Grey50","Red"),
                       values=rescale(c(scaleLow,scaleMid,scaleHigh)),
                       limits=c(scaleLow,scaleHigh),
                       guide = guide_colorbar(title = "-sign(cor)*P_log10"))+
  geom_text(aes(label=ID),size=3, colour="White", fontface="bold")+
  theme_bw()+
  scale_x_continuous(breaks=min(mf_flagellar_assembly_df$rank):max(mf_flagellar_assembly_df$rank))+
  theme(axis.line.y = element_blank(),
        legend.position="bottom",
        axis.title.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())

print(fig02)
###*****************************


###*****************************
# Save figure 1
colWidth=0.6*ifelse(max(kegg_flagellar_assembly_df$rank)-min(kegg_flagellar_assembly_df$rank)+1<8, 
                8, 
                max(kegg_flagellar_assembly_df$rank)-min(kegg_flagellar_assembly_df$rank))
rowWidth=3.5

cowplot::save_plot(filename = paste0("../d_figures/kegg_flagellar_assembly.pdf"),
                   plot = fig01,
                   base_height = rowWidth,
                   base_width = colWidth,
                   limitsize = FALSE)


# Save figure 2
colWidth=ifelse(max(mf_flagellar_assembly_df$rank)-min(mf_flagellar_assembly_df$rank)+1<8, 
                8, 
                max(mf_flagellar_assembly_df$rank)-min(mf_flagellar_assembly_df$rank))
rowWidth=3

cowplot::save_plot(filename = paste0("../d_figures/mf_flagellar_assembly.pdf"),
                   plot = fig01,
                   base_height = rowWidth,
                   base_width = colWidth,
                   limitsize = FALSE)
###*****************************