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
require("VennDiagram")
require("gtable")
require("gridExtra")
###*****************************


###*****************************
# Get Data in 4 parts

resultList=dir("../c_results/DeSeq2_diffGene_batch_Results/")
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
  temp=read.csv(file = paste0("../c_results/DeSeq2_diffGene_batch_Results/",
                              resultList_Exp_mrna[counter01]))
  if(counter01==1){df_combined=temp}
  if(counter01!=1){df_combined=rbind(df_combined,temp)}
}

df_combined %>%
  dplyr::filter(!is.na(padj))%>%
  dplyr::group_by(contrast,base,signChange) %>%
  dplyr::filter(padj<0.05 & abs(log2FoldChange)>1)%>%
  dplyr::summarise(P0.05Fold2=length(padj),
                   growthPhase=unique(growthPhase),
                   data_type=unique(pick_data))->df_summary_exp_mrna

df_combined %>%
  dplyr::filter(base=="baseNa", 
                padj<0.05, abs(log2FoldChange)>1)->temp
exp_mrna_Na=unique(as.vector(temp$gene_name))
df_combined %>%
  dplyr::filter(base == "baseMg",
                padj<0.05, abs(log2FoldChange)>1)->temp
exp_mrna_Mg=unique(as.vector(temp$gene_name))
df_combined %>%
  dplyr::filter(base == "glucose", 
                padj<0.05, abs(log2FoldChange)>1)->temp
exp_mrna_Carb=unique(as.vector(temp$gene_name))
###*****************************


###*****************************
# b) Sta // mrna
for(counter01 in 1: length(resultList_Sta_mrna))
{
  temp=read.csv(file = paste0("../c_results/DeSeq2_diffGene_batch_Results//",
                              resultList_Sta_mrna[counter01]))
  if(counter01==1){df_combined=temp}
  if(counter01!=1){df_combined=rbind(df_combined,temp)}
}


df_combined %>%
  dplyr::filter(!is.na(padj))%>%
  dplyr::group_by(contrast,base,signChange) %>%
  dplyr::filter(padj<0.05 & abs(log2FoldChange)>1)%>%
  dplyr::summarise(P0.05Fold2=length(padj),
                   growthPhase=unique(growthPhase),
                   data_type=unique(pick_data))->df_summary_sta_mrna

df_combined %>%
  dplyr::filter(base == "baseNa",
                padj<0.05, abs(log2FoldChange)>1)->temp
sta_mrna_Na=unique(as.vector(temp$gene_name))
df_combined %>%
  dplyr::filter(base == "baseMg",
                padj<0.05, abs(log2FoldChange)>1)->temp
sta_mrna_Mg=unique(as.vector(temp$gene_name))
df_combined %>%
  dplyr::filter(base == "glucose", 
                padj<0.05, abs(log2FoldChange)>1)->temp
sta_mrna_Carb=unique(as.vector(temp$gene_name))
###*****************************


###*****************************
# c) Exp // protein
for(counter01 in 1: length(resultList_Exp_protein))
{
  temp=read.csv(file = paste0("../c_results/DeSeq2_diffGene_batch_Results//",
                              resultList_Exp_protein[counter01]))
  if(counter01==1){df_combined=temp}
  if(counter01!=1){df_combined=rbind(df_combined,temp)}
}

df_combined %>%
  dplyr::filter(!is.na(padj))%>%
  dplyr::group_by(contrast,base,signChange) %>%
  dplyr::filter(padj<0.05 & abs(log2FoldChange)>1)%>%
  dplyr::summarise(P0.05Fold2=length(padj),
                   growthPhase=unique(growthPhase),
                   data_type=unique(pick_data))->df_summary_exp_protein

df_combined %>%
  dplyr::filter(base == "baseNa", 
                padj<0.05, abs(log2FoldChange)>1)->temp
exp_protein_Na=unique(as.vector(temp$gene_name))
df_combined %>%
  dplyr::filter(base == "baseMg",
                padj<0.05, abs(log2FoldChange)>1)->temp
exp_protein_Mg=unique(as.vector(temp$gene_name))
df_combined %>%
  dplyr::filter(base == "glucose", 
                padj<0.05, abs(log2FoldChange)>1)->temp
exp_protein_Carb=unique(as.vector(temp$gene_name))
###*****************************


###*****************************
# d) Sta // protein
for(counter01 in 1: length(resultList_Sta_protein))
{
  temp=read.csv(file = paste0("../c_results/DeSeq2_diffGene_batch_Results//",
                              resultList_Sta_protein[counter01]))
  if(counter01==1){df_combined=temp}
  if(counter01!=1){df_combined=rbind(df_combined,temp)}
}

df_combined %>%
  dplyr::filter(!is.na(padj))%>%
  dplyr::group_by(contrast,base,signChange) %>%
  dplyr::filter(padj<0.05 & abs(log2FoldChange)>1)%>%
  dplyr::summarise(P0.05Fold2=length(padj),
                   growthPhase=unique(growthPhase),
                   data_type=unique(pick_data))->df_summary_sta_protein

df_combined %>%
  dplyr::filter(base == "baseNa", 
                padj<0.05, abs(log2FoldChange)>1)->temp
sta_protein_Na=unique(as.vector(temp$gene_name))
df_combined %>%
  dplyr::filter(base == "baseMg",
                padj<0.05, abs(log2FoldChange)>1)->temp
sta_protein_Mg=unique(as.vector(temp$gene_name))
df_combined %>%
  dplyr::filter(base == "glucose", 
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
  dplyr::mutate(investigated_conditions=ifelse(contrast == "lactate","Lac",NA),
                investigated_conditions=ifelse(contrast == "highNa","High Na",investigated_conditions),
                investigated_conditions=ifelse(contrast == "highMg","High Mg",investigated_conditions),
                investigated_conditions=ifelse(contrast == "lowMg","Low Mg",investigated_conditions),
                investigated_conditions=ifelse(contrast == "gluconate","Glc",investigated_conditions),
                investigated_conditions=ifelse(contrast == "glycerol" ,"Gly",investigated_conditions))->df_summary


df_summary$signChange <- factor(df_summary$signChange, levels = c("-1", "1"))
df_summary$investigated_conditions <- factor(df_summary$investigated_conditions, levels = c("Low Mg", 
                                                                                            "High Mg",
                                                                                            "High Na",
                                                                                            "Gly",
                                                                                            "Glc",
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
cowplot::save_plot(filename = paste0("../c_figures/difExpressedGenesBatch",".pdf"),
                   plot = figBarGraph01, 
                   ncol = 2, nrow = 2, 
                   base_height = 3.075,base_width = 4.5,
                   units = "in",useDingbats=FALSE)

cowplot::save_plot(filename = paste0("../c_figures/difExpressedGenesBatch_mrna",".pdf"),
                   plot = figBarGraph02a, 
                   ncol = 2, nrow = 1, 
                   base_height = 4,base_width = 4.5,
                   units = "in",useDingbats=FALSE)

cowplot::save_plot(filename = paste0("../c_figures/difExpressedGenesBatch_protein",".pdf"),
                   plot = figBarGraph02b, 
                   ncol = 2, nrow = 1, 
                   base_height = 4,base_width = 4.5,
                   units = "in",useDingbats=FALSE)

###*****************************
###*****************************


###*****************************
# Font sizes for Venn Diagrams
categoryNameFontSize=1.8
mainTextFontSizeVector=1.6
###*****************************


###*****************************
###*****************************
# INDIVIDUAL VENN DIAGRAM

exp_mrna_list = list("Carbon \nsource"=exp_mrna_Carb, "Mg stress"=exp_mrna_Mg, "Na stress"=exp_mrna_Na)
venn.grid = venn.diagram(x=exp_mrna_list, filename="../c_figures/exp_mrna_batch_venn.jpeg", 
                         euler.d=FALSE,scaled=FALSE,
                         print.mode=c("raw","percent"),force.unique=TRUE, 
                         fill=c("purple","cyan","orange"), 
                         fontface = "bold", cex=rep(mainTextFontSizeVector,7), cat.fontface="bold", 
                         fontfamily="Helvetica", main.fontfamily =  "Helvetica", sub.fontfamily = "Helvetica", cat.fontfamily = "Helvetica",
                         cat.cex=c(categoryNameFontSize,categoryNameFontSize,categoryNameFontSize),
                         cat.dist=c(0.1, 0.08, 0.05) ,margin = 0.04)

exp_protein_list = list("Carbon \nsource"=exp_protein_Carb, "Mg stress"=exp_protein_Mg, "Na stress"=exp_protein_Na)
venn.grid = venn.diagram(x=exp_protein_list, filename="../c_figures/exp_protein_batch_venn.jpeg", 
                         euler.d=FALSE,scaled=FALSE,
                         print.mode=c("raw","percent"),force.unique=TRUE, 
                         fill=c("purple","cyan","orange"), 
                         fontface = "bold", cex=rep(mainTextFontSizeVector,7), cat.fontface="bold", 
                         fontfamily="Helvetica", main.fontfamily =  "Helvetica", sub.fontfamily = "Helvetica", cat.fontfamily = "Helvetica",
                         cat.cex=c(categoryNameFontSize,categoryNameFontSize,categoryNameFontSize),
                         cat.dist=c(0.1, 0.08, 0.05) ,margin = 0.04)

sta_mrna_list = list("Carbon \nsource"=sta_mrna_Carb, "Mg stress"=sta_mrna_Mg, "Na stress"=sta_mrna_Na)
venn.grid = venn.diagram(x=sta_mrna_list, filename="../c_figures/sta_mrna_batch_venn.jpeg", 
                         euler.d=FALSE,scaled=FALSE,
                         print.mode=c("raw","percent"),force.unique=TRUE, 
                         fill=c("purple","cyan","orange"), 
                         fontface = "bold", cex=rep(mainTextFontSizeVector,7), cat.fontface="bold", 
                         fontfamily="Helvetica", main.fontfamily =  "Helvetica", sub.fontfamily = "Helvetica", cat.fontfamily = "Helvetica",
                         cat.cex=c(categoryNameFontSize,categoryNameFontSize,categoryNameFontSize),
                         cat.dist=c(0.1, 0.08, 0.05) ,margin = 0.04)

sta_protein_list = list("Carbon \nsource"=sta_protein_Carb, "Mg stress"=sta_protein_Mg, "Na stress"=sta_protein_Na)
venn.grid = venn.diagram(x=sta_protein_list, filename="../c_figures/sta_protein_batch_venn.jpeg", 
                         euler.d=FALSE,scaled=FALSE,
                         print.mode=c("raw","percent"),force.unique=TRUE, 
                         fill=c("purple","cyan","orange"), 
                         fontface = "bold", cex=rep(mainTextFontSizeVector,7), cat.fontface="bold", 
                         fontfamily="Helvetica", main.fontfamily =  "Helvetica", sub.fontfamily = "Helvetica", cat.fontfamily = "Helvetica",
                         cat.cex=c(categoryNameFontSize,categoryNameFontSize,categoryNameFontSize),
                         cat.dist=c(0.1, 0.08, 0.05) ,margin = 0.04)


###*****************************
###*****************************
# MULTIPLE VENN DIAGRAM

exp_mrna_list = list("Carbon \nsource"=exp_mrna_Carb, "Mg stress"=exp_mrna_Mg, "Na stress"=exp_mrna_Na)
exp_mrna_fig = venn.diagram(x=exp_mrna_list, filename=NULL, 
                            euler.d=FALSE,scaled=FALSE,
                            print.mode=c("raw","percent"),force.unique=TRUE, 
                            fill=c("purple","cyan","orange"), 
                            fontface = "bold", cex=rep(mainTextFontSizeVector,7), cat.fontface="bold", 
                            fontfamily="Helvetica", main.fontfamily =  "Helvetica", sub.fontfamily = "Helvetica", cat.fontfamily = "Helvetica",
                            cat.cex=c(categoryNameFontSize,categoryNameFontSize,categoryNameFontSize),
                            cat.dist=c(0.1, 0.08, 0.05) ,margin = 0.07,
                            main = "A", main.pos = c(0.05,1), main.fontface = "bold", main.cex = 2.5)

exp_protein_list = list("Carbon \nsource"=exp_protein_Carb, "Mg stress"=exp_protein_Mg, "Na stress"=exp_protein_Na)
exp_protein_fig = venn.diagram(x=exp_protein_list, filename=NULL, 
                               euler.d=FALSE,scaled=FALSE,
                               print.mode=c("raw","percent"),force.unique=TRUE, 
                               fill=c("purple","cyan","orange"), 
                               fontface = "bold", cex=rep(mainTextFontSizeVector,7), cat.fontface="bold", 
                               fontfamily="Helvetica", main.fontfamily =  "Helvetica", sub.fontfamily = "Helvetica",  cat.fontfamily = "Helvetica",
                               cat.cex=c(categoryNameFontSize,categoryNameFontSize,categoryNameFontSize),
                               cat.dist=c(0.1, 0.08, 0.05) ,margin = 0.07,
                               main = "B", main.pos = c(0.05,1), main.fontface = "bold", main.cex = 2.5)

sta_mrna_list = list("Carbon \nsource"=sta_mrna_Carb, "Mg stress"=sta_mrna_Mg, "Na stress"=sta_mrna_Na)
sta_mrna_fig = venn.diagram(x=sta_mrna_list, filename=NULL, 
                            euler.d=FALSE,scaled=FALSE,
                            print.mode=c("raw","percent"),force.unique=TRUE, 
                            fill=c("purple","cyan","orange"), 
                            fontface = "bold", cex=rep(mainTextFontSizeVector,7), cat.fontface="bold", 
                            fontfamily="Helvetica", main.fontfamily =  "Helvetica", sub.fontfamily = "Helvetica", cat.fontfamily = "Helvetica",
                            cat.cex=c(categoryNameFontSize,categoryNameFontSize,categoryNameFontSize),
                            cat.dist=c(0.1, 0.08, 0.05) ,margin = 0.07,
                            main = "C", main.pos = c(0.05,1), main.fontface = "bold", main.cex = 2.5)

sta_protein_list = list("Carbon \nsource"=sta_protein_Carb, "Mg stress"=sta_protein_Mg, "Na stress"=sta_protein_Na)
sta_protein_fig = venn.diagram(x=sta_protein_list, filename=NULL,
                               euler.d=FALSE,scaled=FALSE,
                               print.mode=c("raw","percent"),force.unique=TRUE, 
                               fill=c("purple","cyan","orange"), 
                               fontface = "bold", cex=rep(mainTextFontSizeVector,7), cat.fontface="bold", 
                               fontfamily="Helvetica", main.fontfamily =  "Helvetica", sub.fontfamily = "Helvetica", cat.fontfamily = "Helvetica",
                               cat.cex=c(categoryNameFontSize,categoryNameFontSize,categoryNameFontSize),
                               cat.dist=c(0.1, 0.08, 0.05) ,margin = 0.07,
                               main = "D", main.pos = c(0.05,1), main.fontface = "bold", main.cex = 2.5)


#**************************************
png(filename = '../c_figures/venn_batch.png',width = 6600, height = 6600, units="px", res =500)
pushViewport(viewport(layout=grid.layout(ncol=3,nrow = 3, 
                                         widths = unit(c(5,5,1)/11, "npc"), 
                                         heights = unit(c(1,5,5)/11, "npc")
)))
pushViewport(viewport(layout.pos.row = 2, layout.pos.col=1))
grid.draw(exp_mrna_fig)
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col=2))
grid.draw(exp_protein_fig)
popViewport()

pushViewport(viewport(layout.pos.row = 3, layout.pos.col=1))
grid.draw(sta_mrna_fig)
popViewport()

pushViewport(viewport(layout.pos.row = 3, layout.pos.col=2))
grid.draw(sta_protein_fig)
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col=1))
grid.draw(textGrob(label = "mRNA",gp=gpar(fontsize=35)))
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col=2))
grid.draw(textGrob(label = "Protein",gp=gpar(fontsize=35)))
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col=3))
grid.draw(textGrob(label = "Exponential",gp=gpar(fontsize=35),rot = -90))
popViewport()

pushViewport(viewport(layout.pos.row = 3, layout.pos.col=3))
grid.draw(textGrob(label = "Stationary",gp=gpar(fontsize=35),rot = -90))
popViewport()

popViewport(0)
dev.off()
#**************************************


#**************************************
pdf(file = '../c_figures/venn_batch.pdf',width = 14, height = 14)
pushViewport(viewport(layout=grid.layout(ncol=3,nrow = 3, 
                                         widths = unit(c(5,5,1)/11, "npc"), 
                                         heights = unit(c(1,5,5)/11, "npc")
)))
pushViewport(viewport(layout.pos.row = 2, layout.pos.col=1))
grid.draw(exp_mrna_fig)
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col=2))
grid.draw(exp_protein_fig)
popViewport()

pushViewport(viewport(layout.pos.row = 3, layout.pos.col=1))
grid.draw(sta_mrna_fig)
popViewport()

pushViewport(viewport(layout.pos.row = 3, layout.pos.col=2))
grid.draw(sta_protein_fig)
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col=1))
grid.draw(textGrob(label = "mRNA",gp=gpar(fontsize=35)))
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col=2))
grid.draw(textGrob(label = "Protein",gp=gpar(fontsize=35)))
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col=3))
grid.draw(textGrob(label = "Exponential",gp=gpar(fontsize=35),rot = -90))
popViewport()

pushViewport(viewport(layout.pos.row = 3, layout.pos.col=3))
grid.draw(textGrob(label = "Stationary",gp=gpar(fontsize=35),rot = -90))
popViewport()

popViewport(0)
dev.off()
#**************************************

unlink(x=grep(pattern = "*.log",
              x = dir(path = "."),
              value = TRUE))

unlink(x=paste0("../c_figures/",grep(pattern = "*.log",
                                     x = dir(path = "../c_figures/"),
                                     value = TRUE)))

