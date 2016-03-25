# generate figures from DESeq Results

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
require("Biobase") 
require("DESeq2")
require("dplyr")
require("tidyr")
require("ggplot2")
require("cowplot")
###*****************************


###*****************************
# Get Data in 4 parts

resultList=dir("../c_results/DESeq2_diffGene_Results")
resultList=grep(pattern = "resDF*", x = resultList, value = TRUE)

resultList_Exp=grep(pattern = "*_Exp_*", x = resultList, value = TRUE)
resultList_Exp_mrna=grep(pattern = "*_mrna_*", x = resultList_Exp, value = TRUE)
resultList_Exp_protein=grep(pattern = "*_protein_*", x = resultList_Exp, value = TRUE)

resultList_Sta=grep(pattern = "*_Sta_*", x = resultList, value = TRUE)
resultList_Sta_mrna=grep(pattern = "*_mrna_*", x = resultList_Sta, value = TRUE)
resultList_Sta_protein=grep(pattern = "*_protein_*", x = resultList_Sta, value = TRUE)
###*****************************


###*****************************
# a) Exp // mrna
for(counter01 in 1: length(resultList_Exp_mrna))
{
  temp=read.csv(file = paste0("../c_results/DESeq2_diffGene_Results/",
                              resultList_Exp_mrna[counter01]))
  if(counter01==1){df_combined=temp}
  if(counter01!=1){df_combined=rbind(df_combined,temp)}
}

df_combined %>%
  dplyr::filter(!is.na(padj))%>%
  dplyr::group_by(vs,signChange) %>%
  dplyr::filter(padj<0.05 & abs(log2FoldChange)>1)%>%
  dplyr::summarise(P0.05Fold2=length(padj),
                   growthPhase=unique(growthPhase),
                   data_type=unique(pick_data))->df_summary_exp_mrna

df_combined %>%
  dplyr::filter(vs %in% "baseNahighNa", 
                padj<0.05, abs(log2FoldChange)>1)->temp
exp_mrna_Na=unique(as.vector(temp$gene_name))
df_combined %>%
  dplyr::filter(vs %in% c("baseMghighMg" , "baseMglowMg"),
                padj<0.05, abs(log2FoldChange)>1)->temp
exp_mrna_Mg=unique(as.vector(temp$gene_name))
df_combined %>%
  dplyr::filter(vs %in% c("glucoselactate","glucosegluconate","glucoseglycerol"), 
                padj<0.05, abs(log2FoldChange)>1)->temp
exp_mrna_Carb=unique(as.vector(temp$gene_name))
###*****************************


###*****************************
# b) Sta // mrna
for(counter01 in 1: length(resultList_Sta_mrna))
{
  temp=read.csv(file = paste0("../c_results/DESeq2_diffGene_Results/",
                              resultList_Sta_mrna[counter01]))
  if(counter01==1){df_combined=temp}
  if(counter01!=1){df_combined=rbind(df_combined,temp)}
}


df_combined %>%
  dplyr::filter(!is.na(padj))%>%
  dplyr::group_by(vs,signChange) %>%
  dplyr::filter(padj<0.05 & abs(log2FoldChange)>1)%>%
  dplyr::summarise(P0.05Fold2=length(padj),
                   growthPhase=unique(growthPhase),
                   data_type=unique(pick_data))->df_summary_sta_mrna

df_combined %>%
  dplyr::filter(vs %in% "baseNahighNa", 
                padj<0.05, abs(log2FoldChange)>1)->temp
sta_mrna_Na=unique(as.vector(temp$gene_name))
df_combined %>%
  dplyr::filter(vs %in% c("baseMghighMg" , "baseMglowMg"),
                padj<0.05, abs(log2FoldChange)>1)->temp
sta_mrna_Mg=unique(as.vector(temp$gene_name))
df_combined %>%
  dplyr::filter(vs %in% c("glucoselactate","glucosegluconate","glucoseglycerol"), 
                padj<0.05, abs(log2FoldChange)>1)->temp
sta_mrna_Carb=unique(as.vector(temp$gene_name))
###*****************************


###*****************************
# c) Exp // protein
for(counter01 in 1: length(resultList_Exp_protein))
{
  temp=read.csv(file = paste0("../c_results/DESeq2_diffGene_Results/",
                              resultList_Exp_protein[counter01]))
  if(counter01==1){df_combined=temp}
  if(counter01!=1){df_combined=rbind(df_combined,temp)}
}

df_combined %>%
  dplyr::filter(!is.na(padj))%>%
  dplyr::group_by(vs,signChange) %>%
  dplyr::filter(padj<0.05 & abs(log2FoldChange)>1)%>%
  dplyr::summarise(P0.05Fold2=length(padj),
                   growthPhase=unique(growthPhase),
                   data_type=unique(pick_data))->df_summary_exp_protein

df_combined %>%
  dplyr::filter(vs %in% "baseNahighNa", 
                padj<0.05, abs(log2FoldChange)>1)->temp
exp_protein_Na=unique(as.vector(temp$gene_name))
df_combined %>%
  dplyr::filter(vs %in% c("baseMghighMg" , "baseMglowMg"),
                padj<0.05, abs(log2FoldChange)>1)->temp
exp_protein_Mg=unique(as.vector(temp$gene_name))
df_combined %>%
  dplyr::filter(vs %in% c("glucoselactate","glucosegluconate","glucoseglycerol"), 
                padj<0.05, abs(log2FoldChange)>1)->temp
exp_protein_Carb=unique(as.vector(temp$gene_name))
###*****************************


###*****************************
# d) Sta // protein
for(counter01 in 1: length(resultList_Sta_protein))
{
  temp=read.csv(file = paste0("../c_results/DESeq2_diffGene_Results/",
                              resultList_Sta_protein[counter01]))
  if(counter01==1){df_combined=temp}
  if(counter01!=1){df_combined=rbind(df_combined,temp)}
}

df_combined %>%
  dplyr::filter(!is.na(padj))%>%
  dplyr::group_by(vs,signChange) %>%
  dplyr::filter(padj<0.05 & abs(log2FoldChange)>1)%>%
  dplyr::summarise(P0.05Fold2=length(padj),
                   growthPhase=unique(growthPhase),
                   data_type=unique(pick_data))->df_summary_sta_protein

df_combined %>%
  dplyr::filter(vs %in% "baseNahighNa", 
                padj<0.05, abs(log2FoldChange)>1)->temp
sta_protein_Na=unique(as.vector(temp$gene_name))
df_combined %>%
  dplyr::filter(vs %in% c("baseMghighMg" , "baseMglowMg"),
                padj<0.05, abs(log2FoldChange)>1)->temp
sta_protein_Mg=unique(as.vector(temp$gene_name))
df_combined %>%
  dplyr::filter(vs %in% c("glucoselactate","glucosegluconate","glucoseglycerol"), 
                padj<0.05, abs(log2FoldChange)>1)->temp
sta_protein_Carb=unique(as.vector(temp$gene_name))
###*****************************


###*****************************
# Combine summary df's to generate figure 

df_summary<-dplyr::rbind_list(df_summary_exp_mrna,
                              df_summary_exp_protein,
                              df_summary_sta_mrna,
                              df_summary_sta_protein)
###*****************************


###*****************************
# Change Column Names
df_summary %>%
  dplyr::mutate(investigated_conditions=ifelse(vs=="glucoselactate","Lac",NA),
                investigated_conditions=ifelse(vs=="baseNahighNa","High Na",investigated_conditions),
                investigated_conditions=ifelse(vs=="baseMghighMg","High Mg",investigated_conditions),
                investigated_conditions=ifelse(vs=="baseMglowMg","Low Mg",investigated_conditions),
                investigated_conditions=ifelse(vs=="glucosegluconate","Glu",investigated_conditions),
                investigated_conditions=ifelse(vs=="glucoseglycerol" ,"Gly",investigated_conditions))->df_summary

df_summary$signChange <- factor(df_summary$signChange, levels = c("-1", "1"))
df_summary$investigated_conditions <- factor(df_summary$investigated_conditions, levels = c("Low Mg", 
                                                                                            "High Mg",
                                                                                            "High Na",
                                                                                            "Gly",
                                                                                            "Glu",
                                                                                            "Lac"))

df_summary %>% 
  dplyr::mutate(growthPhase = ifelse(growthPhase=="Exp", "Exponential", growthPhase),
                growthPhase = ifelse(growthPhase=="Sta", "Stationary", growthPhase))->df_summary
df_summary$growthPhase <- factor(df_summary$growthPhase, levels = c("Exponential", "Stationary"))

df_summary %>% 
  dplyr::mutate(data_type_abv = ifelse(data_type=="mrna", "mRNA", NA),
                data_type_abv = ifelse(data_type=="protein", "Protein", data_type_abv))->df_summary
df_summary$data_type_abv <- factor(df_summary$data_type_abv, levels = c("mRNA", "Protein"))


df_summary %>%
  dplyr::filter(data_type=="mrna") -> df_summary_mrna

df_summary %>%
  dplyr::filter(data_type=="protein") -> df_summary_protein
###*****************************


###*****************************
# Figures
figBarGraph01=ggplot(df_summary, aes(x=investigated_conditions, 
                                   y=P0.05Fold2,
                                   fill=as.factor(signChange))) +
  facet_grid(growthPhase ~ data_type_abv)+
  scale_y_continuous(expand = c(0, 0),limits = c(0,1600))+
  geom_bar(position="dodge",stat="identity",width=.75)+
  scale_fill_manual(values = c("blue","red"),
                    name="Regulation",
                    breaks=c("-1", "1"),
                    labels=c("Down Regulated", "Up Regulated"))+
  theme_bw()+
  xlab("Test Condition") + 
  ylab("Number of differentially expressed genes") +
  theme(panel.margin.y = unit(2, "lines"),
        panel.grid.minor.x = element_blank(),
        legend.position=c(0.88,0.9),
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


print(figBarGraph01)


figBarGraph02a=ggplot(df_summary_mrna, aes(x=investigated_conditions, 
                                     y=P0.05Fold2,
                                     fill=as.factor(signChange))) +
  facet_grid(. ~ growthPhase)+
  geom_bar(position="dodge",stat="identity",width=.75)+
  scale_fill_manual(values = c("blue","red"),
                    name="Regulation",
                    breaks=c("-1", "1"),
                    labels=c("Down Regulated", "Up Regulated"))+
  ylim(0,1700)+
  theme_bw()+
  xlab("Test Condition") + ylab("Count") +
  theme(panel.grid.minor.x = element_blank(),
        legend.position=c(0.88,0.8),
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
print(figBarGraph02a)


figBarGraph02b=ggplot(df_summary_protein, aes(x=investigated_conditions, 
                                           y=P0.05Fold2,
                                           fill=as.factor(signChange))) +
  facet_grid(. ~ growthPhase)+
  geom_bar(position="dodge",stat="identity",width=.75)+
  scale_fill_manual(values = c("blue","red"),
                    name="Regulation",
                    breaks=c("-1", "1"),
                    labels=c("Down Regulated", "Up Regulated"))+
  ylim(0,1700)+
  theme_bw()+
  xlab("Test Condition") + ylab("Count") +
  theme(panel.grid.minor.x = element_blank(),
        legend.position=c(0.88,0.8),
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
print(figBarGraph02b)
###*****************************


###*****************************
# Save Figure
cowplot::save_plot(filename = paste0("../c_figures/difExpressedGenes",".pdf"),
                   plot = figBarGraph01, 
                   ncol = 2, nrow = 2, 
                   base_height = 3.075,base_width = 4.5,
                   units = "in",useDingbats=FALSE)

cowplot::save_plot(filename = paste0("../c_figures/difExpressedGenes_mrna",".pdf"),
                   plot = figBarGraph02a, 
                   ncol = 2, nrow = 1, 
                   base_height = 4,base_width = 4.5,
                   units = "in",useDingbats=FALSE)

cowplot::save_plot(filename = paste0("../c_figures/difExpressedGenes_protein",".pdf"),
                   plot = figBarGraph02b, 
                   ncol = 2, nrow = 1, 
                   base_height = 4,base_width = 4.5,
                   units = "in",useDingbats=FALSE)

###*****************************
###*****************************




###*****************************
###*****************************
# VENN DIAGRAM

# Venn Diagram Count Function

VennDiagramCounter<-function(a,b,c)
{
  l_a=length(a)
  l_b=length(b)
  l_c=length(c)
  l_ab=length(intersect(a,b))
  l_ac=length(intersect(a,c))
  l_bc=length(intersect(b,c))
  l_abc=length(intersect(intersect(a,b),c))
  
  Pos_5=l_abc
  Pos_2=l_ab-l_abc
  Pos_4=l_ac-l_abc
  Pos_6=l_bc-l_abc
  Pos_1=l_a-l_ab-l_ac+l_abc
  Pos_3=l_b-l_ab-l_bc+l_abc
  Pos_7=l_c-l_ac-l_bc+l_abc
  
  Pos=c(Pos_1=Pos_1, 
        Pos_2=Pos_2, 
        Pos_3=Pos_3, 
        Pos_4=Pos_4, 
        Pos_5=Pos_5, 
        Pos_6=Pos_6, 
        Pos_7=Pos_7)
  PosPercent=c(Pos_1=Pos_1, 
               Pos_2=Pos_2, 
               Pos_3=Pos_3, 
               Pos_4=Pos_4, 
               Pos_5=Pos_5, 
               Pos_6=Pos_6, 
               Pos_7=Pos_7)/sum(Pos)
  PosPercent_a=c(Pos_1=Pos_1, 
                 Pos_2=Pos_2, 
                 Pos_4=Pos_4, 
                 Pos_5=Pos_5)/sum(c(Pos_1, Pos_2, Pos_4, Pos_5))
  PosPercent_b=c(Pos_2=Pos_2, 
                 Pos_3=Pos_3, 
                 Pos_5=Pos_5, 
                 Pos_6=Pos_6)/sum(c(Pos_2, Pos_3, Pos_5, Pos_6))
  PosPercent_c=c(Pos_4=Pos_4, 
                 Pos_5=Pos_5, 
                 Pos_6=Pos_6, 
                 Pos_7=Pos_7)/sum(c(Pos_4, Pos_5, Pos_6, Pos_7))
  output=list(Pos=Pos,
              PosPercent=PosPercent, 
              PosPercent_a = PosPercent_a, 
              PosPercent_b = PosPercent_b, 
              PosPercent_c = PosPercent_c)
  return(output)
}
###*****************************


###*****************************
venn_exp_mrna=VennDiagramCounter(b=exp_mrna_Mg, c=exp_mrna_Na, a=exp_mrna_Carb)
venn_sta_mrna=VennDiagramCounter(b=sta_mrna_Mg, c=sta_mrna_Na, a=sta_mrna_Carb)

venn_exp_protein=VennDiagramCounter(b=exp_protein_Mg, c=exp_protein_Na, a=exp_protein_Carb)
venn_sta_protein=VennDiagramCounter(b=sta_protein_Mg, c=sta_protein_Na, a=sta_protein_Carb)
###*****************************















