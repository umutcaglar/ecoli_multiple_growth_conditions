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
    #print(unique((temp[,c("pick_data","growthPhase","test_for","vs")])))
    #browser()
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

# Generate flegellar data frame kegg
mainLoadedFile_kegg$vs=as.character(mainLoadedFile_kegg$vs)
mainLoadedFile_kegg %>%
  dplyr::filter(pick_data=="mrna",test_for!="carbonSource")%>%
  dplyr::filter(KEGG_Path_short %in% kegg_vector) %>%
  dplyr::mutate(pick_data=ifelse(pick_data=="mrna","mrna","protein"))%>%
  dplyr::mutate(vs=ifelse(vs=="baseNahighNa","highNa",vs),
                vs=ifelse(vs=="baseMghighMg","highMg",vs),
                vs=ifelse(vs=="baseMglowMg","lowMg",vs))%>%
  dplyr::group_by(KEGG_Path, padj_gene, pick_data,growthPhase,vs) %>%
  dplyr::mutate(grouping=paste(growthPhase,
                               vs,sep = "-"))->kegg_metabolism_df

kegg_metabolism_df %>%
  dplyr::group_by(grouping,signChange) %>%
  dplyr::summarize(lengthObj=length(score_gene))%>%
  dplyr::group_by() %>%
  tidyr::complete(grouping,signChange)->kegg_metabolism_summary

kegg_metabolism_summary$grouping <- factor(kegg_metabolism_summary$grouping, 
                                           levels = rev(c("Exp-lowMg", "Exp-highMg", "Exp-highNa",
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
mainLoadedFile_mf$vs=as.character(mainLoadedFile_mf$vs)
mainLoadedFile_mf %>%
  dplyr::filter(MF_Name_short %in% mf_vector) %>%
  dplyr::filter(pick_data=="mrna",test_for!="carbonSource")%>%
  dplyr::mutate(pick_data=ifelse(pick_data=="mrna","mrna","protein"))%>%
  dplyr::mutate(vs=ifelse(vs=="baseNahighNa","highNa",vs),
                vs=ifelse(vs=="baseMghighMg","highMg",vs),
                vs=ifelse(vs=="baseMglowMg","lowMg",vs))%>%
  dplyr::group_by(MF_Name, padj_gene, pick_data,growthPhase,vs) %>%
  dplyr::mutate(grouping=paste(growthPhase,
                               vs,sep = "-"))->mf_metabolism_df

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
# Generate metabolism figure kegg
fig01=ggplot(kegg_metabolism_summary, aes(x=grouping,
                                          y=lengthObj,
                                          fill=as.factor(signChange))) +
  geom_bar(position="dodge", stat="identity",width=.75)+
  scale_fill_manual(values = c("blue","red"),
                    name="Regulation",
                    breaks=c("-1", "1"),
                    labels=c("Down Regulated", "Up Regulated"))+
  xlab("conditions")+
  ylab("number of differentially expressed")+
  theme_bw()

print(fig01)


fig02=ggplot(mf_metabolism_summary, aes(x=grouping,
                                        y=lengthObj,
                                        fill=as.factor(signChange))) +
  geom_bar(position="dodge", stat="identity",width=.75)+
  scale_fill_manual(values = c("blue","red"),
                    name="Regulation",
                    breaks=c("-1", "1"),
                    labels=c("Down Regulated", "Up Regulated"))+
  xlab("conditions")+
  ylab("number of differentially expressed")+
  theme_bw()

print(fig02)

fig03<-ggplot(metaData,aes(x=Mg_mM_Levels, y=doublingTimeMinutes))+
  geom_point()+
  ylim(0,90)+
  theme_bw()+
  xlab("conditions")+
  ylab("doubling time")+
  theme(panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank())

print(fig03)
###*****************************



###*****************************
# Save Figures
cowplot::save_plot(filename = paste0("../d_figures/metabolism_kegg.pdf"),
                   plot = fig01,ncol = 2,limitsize = FALSE)

cowplot::save_plot(filename = paste0("../d_figures/metabolism_growthtime.pdf"),
                   plot = fig03,ncol = 1,limitsize = FALSE)