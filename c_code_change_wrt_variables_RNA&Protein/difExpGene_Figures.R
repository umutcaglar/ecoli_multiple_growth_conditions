# Generate diffreential analysis of count data by using DESeq2 package

###*****************************
# INITIAL COMMANDS TO RESET THE SYSTEM

rm(list = ls())
if (is.integer(dev.list())){dev.off()}
cat("\014")

# Set Working Directory
setwd('/Users/umut/GitHub/ecoli_multiple_growth_conditions/c_code_change_wrt_variables_RNA&Protein/') # mac computer
###*****************************

###*****************************
# DOWNLOAD LIBRARIES
library("Biobase")
library("DESeq2")
library("dplyr")
library("ggplot2")
require('gridExtra')
require("cowplot")

###*****************************


###*****************************
# LOAD DATA
load(file = "../c_results/significantChanges.Rdata")
dictionaryRNA=read.csv(paste0("../generateDictionary/","nameDictionary_RNA_barrick.csv"))
dictionaryProtein=read.csv(paste0("../generateDictionary/","nameDictionary_Protein.csv"))
###*****************************


###*****************************
# RNA Data
resHighMgRNAExp %>%
  dplyr::mutate(gene_id=row.names(resHighMgRNAExp),
                signChange=sign(log2FoldChange),
                dataName="resHighMgRNA")%>%
  left_join(.,dictionaryRNA)->resHighMgRNAExp
resLowMgRNAExp%>%
  dplyr::mutate(gene_id=row.names(resLowMgRNAExp),
                signChange=sign(log2FoldChange),
                dataName="resLowMgRNA")%>%
  left_join(.,dictionaryRNA)->resLowMgRNAExp
resHighNaRNAExp%>%
  dplyr::mutate(gene_id=row.names(resHighNaRNAExp),
                signChange=sign(log2FoldChange),
                dataName="resHighNaRNA")%>%
  left_join(.,dictionaryRNA)->resHighNaRNAExp
resGluRNAExp%>%
  dplyr::mutate(gene_id=row.names(resGluRNAExp),
                signChange=sign(log2FoldChange),
                dataName="resGluRNA")%>%
  left_join(.,dictionaryRNA)->resGluRNAExp
resGlyRNAExp%>%
  dplyr::mutate(gene_id=row.names(resGlyRNAExp),
                signChange=sign(log2FoldChange),
                dataName="resGlyRNA")%>%
  left_join(.,dictionaryRNA)->resGlyRNAExp
resLacRNAExp%>%
  dplyr::mutate(gene_id=row.names(resLacRNAExp),
                signChange=sign(log2FoldChange),
                dataName="resLacRNA")%>%
  left_join(.,dictionaryRNA)->resLacRNAExp
###*****************************


###*****************************
# Count Significant Changes
RNAResultsExp=rbind_list(resHighMgRNAExp,resLowMgRNAExp,resHighNaRNAExp,resGluRNAExp,resGlyRNAExp,resLacRNAExp)
RNAResultsExp %>%
  dplyr::filter(padj<0.05) %>%
  dplyr::group_by(signChange,dataName) ->RNAResultsExp
###*****************************


###*****************************
# RNA Data
resHighMgProteinExp %>%
  dplyr::mutate(Protein_id=row.names(resHighMgProteinExp),
                signChange=sign(log2FoldChange),
                dataName="resHighMgProtein")%>%
  left_join(.,dictionaryProtein)->resHighMgProteinExp
resLowMgProteinExp%>%
  dplyr::mutate(Protein_id=row.names(resLowMgProteinExp),
                signChange=sign(log2FoldChange),
                dataName="resLowMgProtein")%>%
  left_join(.,dictionaryProtein)->resLowMgProteinExp
resHighNaProteinExp%>%
  dplyr::mutate(Protein_id=row.names(resHighNaProteinExp),
                signChange=sign(log2FoldChange),
                dataName="resHighNaProtein")%>%
  left_join(.,dictionaryProtein)->resHighNaProteinExp
resGlyProteinExp%>%
  dplyr::mutate(Protein_id=row.names(resGlyProteinExp),
                signChange=sign(log2FoldChange),
                dataName="resGlyProtein")%>%
  left_join(.,dictionaryProtein)->resGlyProteinExp
###*****************************


###*****************************
# Count Significant Changes
ProteinResultsExp=rbind_list(resHighMgProteinExp,resLowMgProteinExp,resHighNaProteinExp,resGlyProteinExp)
ProteinResultsExp %>%
  dplyr::filter(padj<0.05) %>%
  dplyr::group_by(signChange,dataName) ->ProteinResultsExp
###*****************************


###*****************************
# Figure RNA
RNAResultsExp$dataName <- factor(RNAResultsExp$dataName, levels = c("resLowMgRNA",
                                                                    "resHighMgRNA",
                                                                    "resHighNaRNA",
                                                                    "resGlyRNA",
                                                                    "resGluRNA",
                                                                    "resLacRNA"))
figExpRNA=ggplot(RNAResultsExp, aes(x=dataName, fill=as.factor(signChange))) +
  geom_bar(position="dodge")+
  scale_fill_manual(values = c("blue","red"),
                    name="mRNAs",
                    breaks=c("-1", "1"),
                    labels=c("Down Regulated", "Up Regulated"))+
  theme_bw()+
  xlab("Test Condition") + ylab("Count") +
  scale_x_discrete(labels = c("resGluRNA"="Glu","resGlyRNA" = "Gly","resHighMgRNA" = "High Mg",
                              "resHighNaRNA" = "High Na","resLacRNA"="Lac","resLowMgRNA" = "Low Mg"))+
  theme(panel.grid.minor.x = element_blank(),
        legend.position="none")+
  ylim(0, 1500)

print(figExpRNA)
###*****************************


###*****************************
# Figure Protein
ProteinResultsExp$dataName <- factor(ProteinResultsExp$dataName, levels = c("resLowMgProtein",
                                                                            "resHighMgProtein",
                                                                            "resHighNaProtein",
                                                                            "resGlyProtein"))
figExpPro=ggplot(ProteinResultsExp, aes(x=dataName, fill=as.factor(signChange))) +
  geom_bar(position="dodge")+
  scale_fill_manual(values = c("blue","red"),
                    name="Proteins",
                    breaks=c("-1", "1"),
                    labels=c("Down Regulated", "Up Regulated"))+
  theme_bw()+
  xlab("Test Condition") + ylab("Count")+
  scale_x_discrete(labels = c("resGlyProtein" = "Gly","resHighMgProtein" = "High Mg",
                              "resHighNaProtein" = "High Na","resLowMgProtein" = "Low Mg"))+
  theme(panel.grid.minor.x = element_blank(),
        legend.position="none")+
  ylim(0, 1500)
print(figExpPro)
###*****************************


###*****************************
###*****************************


###*****************************
###*****************************
# RNA Data
resHighMgRNASta %>%
  dplyr::mutate(gene_id=row.names(resHighMgRNASta),
                signChange=sign(log2FoldChange),
                dataName="resHighMgRNA")%>%
  left_join(.,dictionaryRNA)->resHighMgRNASta
resLowMgRNASta%>%
  dplyr::mutate(gene_id=row.names(resLowMgRNASta),
                signChange=sign(log2FoldChange),
                dataName="resLowMgRNA")%>%
  left_join(.,dictionaryRNA)->resLowMgRNASta
resHighNaRNASta%>%
  dplyr::mutate(gene_id=row.names(resHighNaRNASta),
                signChange=sign(log2FoldChange),
                dataName="resHighNaRNA")%>%
  left_join(.,dictionaryRNA)->resHighNaRNASta
resGluRNASta%>%
  dplyr::mutate(gene_id=row.names(resGluRNASta),
                signChange=sign(log2FoldChange),
                dataName="resGluRNA")%>%
  left_join(.,dictionaryRNA)->resGluRNASta
resGlyRNASta%>%
  dplyr::mutate(gene_id=row.names(resGlyRNASta),
                signChange=sign(log2FoldChange),
                dataName="resGlyRNA")%>%
  left_join(.,dictionaryRNA)->resGlyRNASta
resLacRNASta%>%
  dplyr::mutate(gene_id=row.names(resLacRNASta),
                signChange=sign(log2FoldChange),
                dataName="resLacRNA")%>%
  left_join(.,dictionaryRNA)->resLacRNASta
###*****************************


###*****************************
# Count Significant Changes
RNAResultsSta=rbind_list(resHighMgRNASta,resLowMgRNASta,resHighNaRNASta,resGluRNASta,resGlyRNASta,resLacRNASta)
RNAResultsSta %>%
  dplyr::filter(padj<0.05) %>%
  dplyr::group_by(signChange,dataName) ->RNAResultsSta
###*****************************


###*****************************
# RNA Data
resHighMgProteinSta %>%
  dplyr::mutate(Protein_id=row.names(resHighMgProteinSta),
                signChange=sign(log2FoldChange),
                dataName="resHighMgProtein")%>%
  left_join(.,dictionaryProtein)->resHighMgProteinSta
resLowMgProteinSta%>%
  dplyr::mutate(Protein_id=row.names(resLowMgProteinSta),
                signChange=sign(log2FoldChange),
                dataName="resLowMgProtein")%>%
  left_join(.,dictionaryProtein)->resLowMgProteinSta
resHighNaProteinSta%>%
  dplyr::mutate(Protein_id=row.names(resHighNaProteinSta),
                signChange=sign(log2FoldChange),
                dataName="resHighNaProtein")%>%
  left_join(.,dictionaryProtein)->resHighNaProteinSta
resGlyProteinSta%>%
  dplyr::mutate(Protein_id=row.names(resGlyProteinSta),
                signChange=sign(log2FoldChange),
                dataName="resGlyProtein")%>%
  left_join(.,dictionaryProtein)->resGlyProteinSta
###*****************************


###*****************************
# Count Significant Changes
ProteinResultsSta=rbind_list(resHighMgProteinSta,resLowMgProteinSta,resHighNaProteinSta,resGlyProteinSta)
ProteinResultsSta %>%
  dplyr::filter(padj<0.05) %>%
  dplyr::group_by(signChange,dataName) ->ProteinResultsSta
###*****************************


###*****************************
# Figure RNA
RNAResultsSta$dataName <- factor(RNAResultsSta$dataName, levels = c("resLowMgRNA",
                                                                    "resHighMgRNA",
                                                                    "resHighNaRNA",
                                                                    "resGlyRNA",
                                                                    "resGluRNA",
                                                                    "resLacRNA"))
figStaRNA=ggplot(RNAResultsSta, aes(x=dataName, fill=as.factor(signChange))) +
  geom_bar(position="dodge")+
  scale_fill_manual(values = c("blue","red"),
                    name="mRNAs",
                    breaks=c("-1", "1"),
                    labels=c("Down Regulated", "Up Regulated"))+
  theme_bw()+
  xlab("Test Condition") + ylab("Count")+
  scale_x_discrete(labels = c("resGluRNA"="Glu","resGlyRNA" = "Gly","resHighMgRNA" = "High Mg",
                              "resHighNaRNA" = "High Na","resLacRNA"="Lac","resLowMgRNA" = "Low Mg"))+
  theme(panel.grid.minor.x = element_blank(),
        legend.position=c(0.7,0.8))+
  ylim(0, 1500)

print(figStaRNA)

###*****************************


###*****************************
# Figure Protein
ProteinResultsSta$dataName <- factor(ProteinResultsSta$dataName, levels = c("resLowMgProtein",
                                             "resHighMgProtein",
                                             "resHighNaProtein",
                                             "resGlyProtein"))
figStaPro=ggplot(ProteinResultsSta, aes(x=dataName, fill=as.factor(signChange))) +
  geom_bar(position="dodge")+
  scale_fill_manual(values = c("blue","red"),
                    name="Proteins",
                    breaks=c("-1", "1"),
                    labels=c("Down Regulated", "Up Regulated"))+
  theme_bw()+
  xlab("Test Condition") + ylab("Count")+
  scale_x_discrete(labels = c("resGlyProtein" = "Gly",
                              "resLowMgProtein" = "Low Mg",
                              "resHighMgProtein" = "High Mg",
                              "resHighNaProtein" = "High Na"))+
  theme(panel.grid.minor.x = element_blank(),
        legend.position="none")+
  ylim(0, 1500)
print(figStaPro)
###*****************************


fig_Counts<-arrangeGrob(figExpRNA, figExpPro, figStaRNA, figStaPro,ncol=2)
print(fig_Counts)
ggsave(fig_Counts,
       filename = paste0("../c_figures/counts",".pdf"),
       width = 8*2,
       height = 6*2,
       units = "in",
       useDingbats=FALSE)

save_plot(filename = paste0("../c_figures/ExpRNACounts",".pdf"),
          plot = figExpRNA,
          base_height = 4.1,
          base_width = 4.1,
          units = "in",
          useDingbats=FALSE)

save_plot(filename = paste0("../c_figures/ExpProCounts",".pdf"),
          plot = figExpPro,
          base_height = 4.1,
          base_width = 4.1,
          units = "in",
          useDingbats=FALSE)

save_plot(filename = paste0("../c_figures/StaRNACounts",".pdf"),
          plot = figStaRNA,
          base_height = 4.1,
          base_width = 4.1,
          units = "in",
          useDingbats=FALSE)

save_plot(filename = paste0("../c_figures/StaProteinCounts",".pdf"),
          plot = figStaPro,
          base_height = 4.1,
          base_width = 4.1,
          units = "in",
          useDingbats=FALSE)


###*****************************
###*****************************
###*****************************
# Generating set data

RNAResultsExp %>%
  dplyr::group_by()%>%
  dplyr::mutate(SetNames=as.character(ifelse(dataName %in% c("resLowMgRNA","resHighMgRNA"),"Mg",NA)),
                SetNames=as.character(ifelse(dataName %in% c("resHighNaRNA"),"Na",SetNames)),
                SetNames=as.character(ifelse(dataName %in% c("resGlyRNA","resGluRNA","resLacRNA"),"Carbon",SetNames)))->RNAResultsExp

RNAResultsExp%>%dplyr::filter(SetNames=="Mg")->temp;
data.frame(unique(temp$gene_name))->temp
temp%>%dplyr::mutate("Mg"=1)->RNAResultsExpMg


RNAResultsExp%>%dplyr::filter(SetNames=="Na")->temp;
data.frame(unique(temp$gene_name))->temp
temp%>%dplyr::mutate("Na"=1)->RNAResultsExpNa

RNAResultsExp%>%dplyr::filter(SetNames=="Carbon")->temp;
data.frame(unique(temp$gene_name))->temp
temp%>%dplyr::mutate("Carbon"=1)->RNAResultsExpCarbon

dplyr::full_join(x = RNAResultsExpMg,y=RNAResultsExpNa)->RNAResultsExpSet
dplyr::full_join(x = RNAResultsExpSet,y=RNAResultsExpCarbon)->RNAResultsExpSet


RNAResultsSta %>%
  dplyr::group_by()%>%
  dplyr::mutate(SetNames=as.character(ifelse(dataName %in% c("resLowMgRNA","resHighMgRNA"),"Mg",NA)),
                SetNames=as.character(ifelse(dataName %in% c("resHighNaRNA"),"Na",SetNames)),
                SetNames=as.character(ifelse(dataName %in% c("resGlyRNA","resGluRNA","resLacRNA"),"Carbon",SetNames)))->RNAResultsSta

RNAResultsSta%>%dplyr::filter(SetNames=="Mg")->temp;
data.frame(unique(temp$gene_name))->temp
temp%>%dplyr::mutate("Mg"=1)->RNAResultsStaMg


RNAResultsSta%>%dplyr::filter(SetNames=="Na")->temp;
data.frame(unique(temp$gene_name))->temp
temp%>%dplyr::mutate("Na"=1)->RNAResultsStaNa

RNAResultsSta%>%dplyr::filter(SetNames=="Carbon")->temp;
data.frame(unique(temp$gene_name))->temp
temp%>%dplyr::mutate("Carbon"=1)->RNAResultsStaCarbon

dplyr::full_join(x = RNAResultsStaMg,y=RNAResultsStaNa)->RNAResultsStaSet
dplyr::full_join(x = RNAResultsStaSet,y=RNAResultsStaCarbon)->RNAResultsStaSet

#****
ProteinResultsExp %>%
  dplyr::group_by()%>%
  dplyr::mutate(SetNames=as.character(ifelse(dataName %in% c("resLowMgProtein","resHighMgProtein"),"Mg",NA)),
                SetNames=as.character(ifelse(dataName %in% c("resHighNaProtein"),"Na",SetNames)),
                SetNames=as.character(ifelse(dataName %in% c("resGlyProtein","resGluProtein","resLacProtein"),"Carbon",SetNames)))->ProteinResultsExp

ProteinResultsExp%>%dplyr::filter(SetNames=="Mg")->temp;
data.frame(unique(temp$gene_name))->temp
temp%>%dplyr::mutate("Mg"=1)->ProteinResultsExpMg


ProteinResultsExp%>%dplyr::filter(SetNames=="Na")->temp;
data.frame(unique(temp$gene_name))->temp
temp%>%dplyr::mutate("Na"=1)->ProteinResultsExpNa

ProteinResultsExp%>%dplyr::filter(SetNames=="Carbon")->temp;
data.frame(unique(temp$gene_name))->temp
temp%>%dplyr::mutate("Carbon"=1)->ProteinResultsExpCarbon

dplyr::full_join(x = ProteinResultsExpMg,y=ProteinResultsExpNa)->ProteinResultsExpSet
dplyr::full_join(x = ProteinResultsExpSet,y=ProteinResultsExpCarbon)->ProteinResultsExpSet


ProteinResultsSta %>%
  dplyr::group_by()%>%
  dplyr::mutate(SetNames=as.character(ifelse(dataName %in% c("resLowMgProtein","resHighMgProtein"),"Mg",NA)),
                SetNames=as.character(ifelse(dataName %in% c("resHighNaProtein"),"Na",SetNames)),
                SetNames=as.character(ifelse(dataName %in% c("resGlyProtein","resGluProtein","resLacProtein"),"Carbon",SetNames)))->ProteinResultsSta

ProteinResultsSta%>%dplyr::filter(SetNames=="Mg")->temp;
data.frame(unique(temp$gene_name))->temp
temp%>%dplyr::mutate("Mg"=1)->ProteinResultsStaMg


ProteinResultsSta%>%dplyr::filter(SetNames=="Na")->temp;
data.frame(unique(temp$gene_name))->temp
temp%>%dplyr::mutate("Na"=1)->ProteinResultsStaNa

ProteinResultsSta%>%dplyr::filter(SetNames=="Carbon")->temp;
data.frame(unique(temp$gene_name))->temp
temp%>%dplyr::mutate("Carbon"=1)->ProteinResultsStaCarbon

dplyr::full_join(x = ProteinResultsStaMg,y=ProteinResultsStaNa)->ProteinResultsStaSet
dplyr::full_join(x = ProteinResultsStaSet,y=ProteinResultsStaCarbon)->ProteinResultsStaSet


#***************************
# Generate tables of genes for individual sections of 3 section venn diagrams

RNAResultsExpSet[is.na(RNAResultsExpSet)]<-0
RNAResultsExpSet%>%filter(Mg==1 & Na==1 & Carbon==1)->temp; nrow(temp)->RNAExp_Mg_Na_Ca
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/RNAExp_Mg_Na_Ca.csv",
            row.names = F,col.names = F,quote = F)
RNAResultsExpSet%>%filter(Mg==0 & Na==1 & Carbon==1)->temp; nrow(temp)->RNAExp_Na_Ca
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/RNAExp_Na_Ca.csv",
            row.names = F,col.names = F,quote = F)
RNAResultsExpSet%>%filter(Mg==1 & Na==0 & Carbon==1)->temp; nrow(temp)->RNAExp_Mg_Ca
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/RNAExp_Mg_Ca.csv",
            row.names = F,col.names = F,quote = F)
RNAResultsExpSet%>%filter(Mg==1 & Na==1 & Carbon==0)->temp; nrow(temp)->RNAExp_Mg_Na
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/RNAExp_Mg_Na.csv",
            row.names = F,col.names = F,quote = F)
RNAResultsExpSet%>%filter(Mg==0 & Na==0 & Carbon==1)->temp; nrow(temp)->RNAExp_Ca
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/RNAExp_Ca.csv",
            row.names = F,col.names = F,quote = F)
RNAResultsExpSet%>%filter(Mg==0 & Na==1 & Carbon==0)->temp; nrow(temp)->RNAExp_Na
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/RNAExp_Na.csv",
            row.names = F,col.names = F,quote = F)
RNAResultsExpSet%>%filter(Mg==1 & Na==0 & Carbon==0)->temp; nrow(temp)->RNAExp_Mg
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/RNAExp_Mg.csv",
            row.names = F,col.names = F,quote = F)

RNAResultsStaSet[is.na(RNAResultsStaSet)]<-0
RNAResultsStaSet%>%filter(Mg==1 & Na==1 & Carbon==1)->temp; nrow(temp)->RNASta_Mg_Na_Ca
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/RNASta_Mg_Na_Ca.csv",
            row.names = F,col.names = F,quote = F)
RNAResultsStaSet%>%filter(Mg==0 & Na==1 & Carbon==1)->temp; nrow(temp)->RNASta_Na_Ca
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/RNASta_Na_Ca.csv",
            row.names = F,col.names = F,quote = F)
RNAResultsStaSet%>%filter(Mg==1 & Na==0 & Carbon==1)->temp; nrow(temp)->RNASta_Mg_Ca
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/RNASta_Mg_Ca.csv",
            row.names = F,col.names = F,quote = F)
RNAResultsStaSet%>%filter(Mg==1 & Na==1 & Carbon==0)->temp; nrow(temp)->RNASta_Mg_Na
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/RNASta_Mg_Na.csv",
            row.names = F,col.names = F,quote = F)
RNAResultsStaSet%>%filter(Mg==0 & Na==0 & Carbon==1)->temp; nrow(temp)->RNASta_Ca
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/RNASta_Ca.csv",
            row.names = F,col.names = F,quote = F)
RNAResultsStaSet%>%filter(Mg==0 & Na==1 & Carbon==0)->temp; nrow(temp)->RNASta_Na
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/RNASta_Na.csv",
            row.names = F,col.names = F,quote = F)
RNAResultsStaSet%>%filter(Mg==1 & Na==0 & Carbon==0)->temp; nrow(temp)->RNASta_Mg
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/RNASta_Mg.csv",
            row.names = F,col.names = F,quote = F)
#***

ProteinResultsExpSet[is.na(ProteinResultsExpSet)]<-0
ProteinResultsExpSet%>%filter(Mg==1 & Na==1 & Carbon==1)->temp; nrow(temp)->ProteinExp_Mg_Na_Ca
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/ProteinExp_Mg_Na_Ca.csv",
            row.names = F,col.names = F,quote = F)
ProteinResultsExpSet%>%filter(Mg==0 & Na==1 & Carbon==1)->temp; nrow(temp)->ProteinExp_Na_Ca
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/ProteinExp_Na_Ca.csv",
            row.names = F,col.names = F,quote = F)
ProteinResultsExpSet%>%filter(Mg==1 & Na==0 & Carbon==1)->temp; nrow(temp)->ProteinExp_Mg_Ca
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/ProteinExp_Mg_Ca.csv",
            row.names = F,col.names = F,quote = F)
ProteinResultsExpSet%>%filter(Mg==1 & Na==1 & Carbon==0)->temp; nrow(temp)->ProteinExp_Mg_Na
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/ProteinExp_Mg_Na.csv",
            row.names = F,col.names = F,quote = F)
ProteinResultsExpSet%>%filter(Mg==0 & Na==0 & Carbon==1)->temp; nrow(temp)->ProteinExp_Ca
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/ProteinExp_Ca.csv",
            row.names = F,col.names = F,quote = F)
ProteinResultsExpSet%>%filter(Mg==0 & Na==1 & Carbon==0)->temp; nrow(temp)->ProteinExp_Na
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/ProteinExp_Na.csv",
            row.names = F,col.names = F,quote = F)
ProteinResultsExpSet%>%filter(Mg==1 & Na==0 & Carbon==0)->temp; nrow(temp)->ProteinExp_Mg
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/ProteinExp_Mg.csv",
            row.names = F,col.names = F,quote = F)

ProteinResultsStaSet[is.na(ProteinResultsStaSet)]<-0
ProteinResultsStaSet%>%filter(Mg==1 & Na==1 & Carbon==1)->temp; nrow(temp)->ProteinSta_Mg_Na_Ca
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/ProteinSta_Mg_Na_Ca.csv",
            row.names = F,col.names = F,quote = F)
ProteinResultsStaSet%>%filter(Mg==0 & Na==1 & Carbon==1)->temp; nrow(temp)->ProteinSta_Na_Ca
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/ProteinSta_Na_Ca.csv",
            row.names = F,col.names = F,quote = F)
ProteinResultsStaSet%>%filter(Mg==1 & Na==0 & Carbon==1)->temp; nrow(temp)->ProteinSta_Mg_Ca
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/ProteinSta_Mg_Ca.csv",
            row.names = F,col.names = F,quote = F)
ProteinResultsStaSet%>%filter(Mg==1 & Na==1 & Carbon==0)->temp; nrow(temp)->ProteinSta_Mg_Na
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/ProteinSta_Mg_Na.csv",
            row.names = F,col.names = F,quote = F)
ProteinResultsStaSet%>%filter(Mg==0 & Na==0 & Carbon==1)->temp; nrow(temp)->ProteinSta_Ca
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/ProteinSta_Ca.csv",
            row.names = F,col.names = F,quote = F)
ProteinResultsStaSet%>%filter(Mg==0 & Na==1 & Carbon==0)->temp; nrow(temp)->ProteinSta_Na
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/ProteinSta_Na.csv",
            row.names = F,col.names = F,quote = F)
ProteinResultsStaSet%>%filter(Mg==1 & Na==0 & Carbon==0)->temp; nrow(temp)->ProteinSta_Mg
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/ProteinSta_Mg.csv",
            row.names = F,col.names = F,quote = F)

#***************************


#***************************
# Generate tables of genes for whole sets
# RNA Exp
RNAResultsExpSet%>%filter(Carbon==1)->temp; nrow(temp)->RNAExp_Ca_Comp
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/RNAExp_Ca_Comp.csv",
            row.names = F,col.names = F,quote = F)
RNAResultsExpSet%>%filter(Na==1)->temp; nrow(temp)->RNAExp_Na_Comp
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/RNAExp_Na_Comp.csv",
            row.names = F,col.names = F,quote = F)
RNAResultsExpSet%>%filter(Mg==1)->temp; nrow(temp)->RNAExp_Mg_Comp
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/RNAExp_Mg_Comp.csv",
            row.names = F,col.names = F,quote = F)


RNAResultsStaSet%>%filter(Carbon==1)->temp; nrow(temp)->RNASta_Ca_Comp
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/RNASta_Ca_Comp.csv",
            row.names = F,col.names = F,quote = F)
RNAResultsStaSet%>%filter(Na==1)->temp; nrow(temp)->RNASta_Na_Comp
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/RNASta_Na_Comp.csv",
            row.names = F,col.names = F,quote = F)
RNAResultsStaSet%>%filter(Mg==1)->temp; nrow(temp)->RNASta_Mg_Comp
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/RNASta_Mg_Comp.csv",
            row.names = F,col.names = F,quote = F)


ProteinResultsExpSet%>%filter(Carbon==1)->temp; nrow(temp)->ProteinExp_Ca_Comp
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/ProteinExp_Ca_Comp.csv",
            row.names = F,col.names = F,quote = F)
ProteinResultsExpSet%>%filter(Na==1)->temp; nrow(temp)->ProteinExp_Na_Comp
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/ProteinExp_Na_Comp.csv",
            row.names = F,col.names = F,quote = F)
ProteinResultsExpSet%>%filter(Mg==1)->temp; nrow(temp)->ProteinExp_Mg_Comp
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/ProteinExp_Mg_Comp.csv",
            row.names = F,col.names = F,quote = F)

ProteinResultsStaSet%>%filter(Carbon==1)->temp; nrow(temp)->ProteinSta_Ca_Comp
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/ProteinSta_Ca_Comp.csv",
            row.names = F,col.names = F,quote = F)
ProteinResultsStaSet%>%filter(Na==1)->temp; nrow(temp)->ProteinSta_Na_Comp
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/ProteinSta_Na_Comp.csv",
            row.names = F,col.names = F,quote = F)
ProteinResultsStaSet%>%filter(Mg==1)->temp; nrow(temp)->ProteinSta_Mg_Comp
write.table(x = as.vector(temp$unique.temp.gene_name.),file = "../c_results/ProteinSta_Mg_Comp.csv",
            row.names = F,col.names = F,quote = F)



