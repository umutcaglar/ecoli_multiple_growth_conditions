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
    #print(unique((temp[,c("pick_data","growthPhase","test_for","contrast")])))
    #browser()
    mainLoadedFile_kegg=rbind(mainLoadedFile_kegg,temp)
  }
}


fileList=dir("../d_results/",pattern = "*.mf_o")

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
# Kegg Vector
kegg_vector=c("Pyruvate metabolism",
              "Oxidative phosphorylation",
              "One carbon pool by folate",
              "Pentose phosphate pathway",
              "Glycolysis / Gluconeogenesis",
              "C5−Branched dibasic acid metabolism",
              "Starch and sucrose metabolism",
              "Fructose and mannose metabolism",
              "Citrate cycle (TCA cycle)",
              "Pentose and glucuronate interconversions",
              "Butanoate metabolism",
              "Galactose metabolism",
              "Pantothenate and CoA biosynthesis",
              "Phenylalanine, tyrosine and tryptophan biosynthesis",
              "Glycine, serine and threonine metabolism",
              "Amino sugar and nucleotide sugar metabolism",
              "Arginine and proline metabolism",
              "Purine metabolism",
              "Pyrimidine metabolism",
              "Amino sugar and nucleotide sugar metabolism",
              "Fatty acid metabolism",
              "Ribosome",
              "Aminoacyl−tRNA biosynthesis",
              "Valine, leucine and isoleucine biosynthesis",
              "Cysteine and methionine metabolism",
              "Alanine, aspartate and glutamate metabolism",
              "Peptidoglycan biosynthesis",
              "Selenoamino acid metabolism",
              "Lysine biosynthesis",
              "Histidine metabolism",
              "Lipopolysaccharide biosynthesis",
              "Sulfur metabolism",
              "structural constituent of ribosome")
###*****************************


###*****************************
# write metabolism related kegg pathways
keggIntersection<-intersect(mainLoadedFile_kegg$KEGG_Path_Short , kegg_vector)
write.csv(x = keggIntersection, file = "../d_results/kegg_metabolism_related_pathways.csv")
###*****************************


###*****************************
# Generate flegellar data frame kegg
mainLoadedFile_kegg %>%
  dplyr::filter(test_for!="carbonSource")%>%
  dplyr::filter(KEGG_Path_Short %in% kegg_vector) %>%
  dplyr::mutate(pick_data=ifelse(pick_data=="mrna","mRNA","Protein"))%>%
  dplyr::group_by(KEGG_Path, padj_gene, pick_data,growthPhase,contrast) %>%
  dplyr::mutate(grouping=paste(growthPhase,
                               contrast,sep = "-"))->kegg_metabolism_df

kegg_metabolism_df %>%
  dplyr::group_by(grouping,signChange,pick_data) %>%
  dplyr::summarize(lengthObj=length(score_gene))%>%
  dplyr::group_by() %>%
  tidyr::complete(grouping,signChange,pick_data)->kegg_metabolism_summary

kegg_metabolism_summary$grouping <- factor(kegg_metabolism_summary$grouping, 
                                           levels = (c("Exp-lowMg", "Exp-highMg", "Exp-highNa",
                                                       "Sta-lowMg", "Sta-highMg", "Sta-highNa")))
###*****************************


###*****************************
# MF Vector

mf_vector=c("Pyruvate metabolism",
            "Oxidative phosphorylation",
            "One carbon pool by folate",
            "Pentose phosphate pathway",
            "Glycolysis / Gluconeogenesis",
            "C5−Branched dibasic acid metabolism",
            "Starch and sucrose metabolism",
            "Fructose and mannose metabolism",
            "Citrate cycle (TCA cycle)",
            "Pentose and glucuronate interconversions",
            "Butanoate metabolism",
            "Galactose metabolism",
            "Pantothenate and CoA biosynthesis",
            "Phenylalanine, tyrosine and tryptophan biosynthesis",
            "Glycine, serine and threonine metabolism",
            "Amino sugar and nucleotide sugar metabolism",
            "Arginine and proline metabolism",
            "Purine metabolism",
            "Pyrimidine metabolism",
            "Amino sugar and nucleotide sugar metabolism",
            "Fatty acid metabolism",
            "Ribosome",
            "Aminoacyl−tRNA biosynthesis",
            "Valine, leucine and isoleucine biosynthesis",
            "Cysteine and methionine metabolism",
            "Alanine, aspartate and glutamate metabolism",
            "Peptidoglycan biosynthesis",
            "Selenoamino acid metabolism",
            "Lysine biosynthesis",
            "Histidine metabolism",
            "Lipopolysaccharide biosynthesis",
            "Sulfur metabolism",
            "structural constituent of ribosome")

# Generate flegellar data frame mf
mainLoadedFile_mf %>%
  dplyr::filter(MF_Short %in% mf_vector) %>%
  dplyr::filter(pick_data=="mrna",test_for!="carbonSource")%>%
  dplyr::mutate(pick_data=ifelse(pick_data=="mrna","mRNA","Protein"))%>%
  dplyr::group_by(MF, padj_gene, pick_data,growthPhase,contrast) %>%
  dplyr::mutate(grouping=paste(growthPhase,
                               contrast,sep = "-"))->mf_metabolism_df

mf_metabolism_df %>%
  dplyr::group_by(grouping,signChange) %>%
  dplyr::summarize(lengthObj=length(score_gene))%>%
  dplyr::group_by() %>%
  tidyr::complete(grouping,signChange)->mf_metabolism_summary

mf_metabolism_summary$grouping <- factor(mf_metabolism_summary$grouping, 
                                         levels = rev(c("Exp-lowMg", "Exp-highMg", "Exp-highNa",
                                                        "Sta-lowMg", "Sta-highMg", "Sta-highNa")))
###*****************************


###*****************************
# rename x-axis
oldLevels=levels(kegg_metabolism_summary$grouping)
newLevels=replace_fun(input_vector = oldLevels, initialVal = oldLevels,
                      finalVal = c("Low Mg \n Exponential","High Mg \n Exponential",
                                   "High Na \n Exponential","Low Mg \n Stationary",
                                   "High Mg \n Stationary","High Na \n Stationary"))
###*****************************


###*****************************
# Generate metabolism figure kegg
fig01=ggplot(kegg_metabolism_summary, aes(x=grouping,
                                          y=lengthObj,
                                          fill=as.factor(signChange))) +
  facet_grid(pick_data~.)+
  geom_bar(position="dodge", stat="identity",width=.75)+
  scale_y_continuous(expand = c(0, 0),limits = c(0,200))+
  scale_x_discrete(labels=newLevels)+
scale_fill_manual(values = c("blue","red"),
                  name="Regulation",
                  breaks=c("-1", "1"),
                  labels=c("Down Regulated", "Up Regulated"))+
  xlab("Conditions")+
  ylab("Number of differentially expressed genes")+
  theme_bw()+
  theme(panel.margin.y = unit(2, "lines"),
        panel.grid.minor.x = element_blank(),
        legend.position=c(0.8,0.9),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        legend.background = element_rect(fill=alpha(colour = "white",
                                                    alpha = .4)))

print(fig01)


fig02=ggplot(mf_metabolism_summary, aes(x=grouping,
                                        y=lengthObj,
                                        fill=as.factor(signChange))) +
  geom_bar(position="dodge", stat="identity",width=.75)+
  scale_fill_manual(values = c("blue","red"),
                    name="Regulation",
                    breaks=c("-1", "1"),
                    labels=c("Down-regulated", "Up-regulated"))+
  xlab("conditions")+
  ylab("Number of differentially expressed genes")+
  theme_bw()

print(fig02)
###*****************************



###*****************************
# Save Figures
cowplot::save_plot(filename = paste0("../d_figures/metabolism_kegg.pdf"),
                   plot = fig01,ncol = 2,nrow=2,limitsize = FALSE)