# Analayze Go Annotations Molecular function DESeq + DAVID

# Aim of the code is to find the genes in kegg pathways and send files to generate figures


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
require("Biobase") 
require("DESeq2")
require("dplyr")
require("tidyr")
require("ggplot2")
require("RColorBrewer")
require("scales")
require("cowplot")
require("ggrepel")
###*****************************


###*****************************
# Load Files
source("../c_code_change_wrt_variables_RNA&Protein/data_naming_functions.R")
###*****************************


###****************************
# Download the DAVID input and output
dataName=name_data(initialValue="genes0.05", # can be "genes0.05"
                   dataType = "mrna", # can be "rna", "mrna", "protein", "protein_wo_NA"
                   badDataSet = "set02", # can be "set00",set01","set02", "set03"
                   # referenceParameters can be a vector like
                   # c("growthPhase", "Mg_mM_Levels", "Na_mM_Levels", "carbonSource", "experiment")
                   referenceParameters=c("growthPhase",
                                         "Mg_mM_Levels", 
                                         "Na_mM_Levels", 
                                         "carbonSource", 
                                         "experiment"),
                   # referenceLevels can be a vector like
                   # c("exponential", "baseMg", "baseNa", "glucose", "glucose_time_course")
                   referenceLevels=c("exponential",
                                     "baseMg", 
                                     "baseNa", 
                                     "glucose", 
                                     "glucose_time_course"),
                   experimentVector = c("allEx"), # can be "Stc","Ytc","Nas","Agr","Ngr","Mgl","Mgh" // "allEx"
                   carbonSourceVector = "S", # can be any sub combination of "SYAN"
                   MgLevelVector = c("baseMg","highMg"), # can be "lowMg","baseMg","highMg" // "allMg"
                   NaLevelVector = c("baseNa"), # can be "baseNa","highNa" // "allNa"
                   growthPhaseVector = c("exponential"), # can be "exponential","stationary","late_stationary" // "allPhase"
                   filterGenes = "noFilter", # can be "noFilter", "meanFilter", "maxFilter", "sdFilter" 
                   threshold=NA, # the threshold value for "meanFilter", "maxFilter", "sdFilter"
                   roundData=TRUE,
                   sumTechnicalReplicates=TRUE,
                   deSeqSfChoice="p1Sf", # can be "regSf", "p1Sf"
                   normalizationMethodChoice= "noNorm", # can be "vst", "rlog", "log10", "noNorm"
                   test_for = "Mg_mM_Levels")  # works only if normalizationMethodChoice == noNorm
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")

objectName=paste(dataName$objectName,collapse = "_")

GoMF_input<-read.csv(file = paste0("../c_results/DeSeq2_diffGene_Results/",objectName,".csv"),
                     header = TRUE)
objectName_df=dataName$objectName
objectName_df$initial="resDf"
objectName_df=paste(objectName_df,collapse = "_")
GoMF_input_df<-read.csv(file = paste0("../c_results/DeSeq2_diffGene_Results/",objectName_df,".csv"),
                        header = TRUE)

GoMF_result<-read.table(file=paste0("../c_results/david_results/",objectName,"_mf.txt"),
                        sep = "\t",header = TRUE)

names(GoMF_result)[names(GoMF_result) == 'Term'] <- 'MF_Name'
###*****************************


###*****************************
# load reference libraries
dataType=dataName$objectName$pick_data

if(dataType %in% c("rna","mrna"))
{
  MF_david_ref_2009_tidy<-read.table(file = paste0("../d_Mf_Pathway_Analyze/ReferenceFiles/",
                                                   "MF_david_mrna_ref_2009_tidy.txt"),
                                     sep = "\t",header = TRUE,fill = TRUE, quote = "")
}


if(dataType %in% c("protein_wo_NA", "protein"))
{
  MF_david_ref_2009_tidy<-read.table(paste0("../d_Mf_Pathway_Analyze/ReferenceFiles/",
                                            "MF_david_protein_ref_2009_tidy.txt"),
                                     sep = "\t",header = TRUE,fill = TRUE, quote = "")
}
###*****************************


###*****************************
# Find MF output in ref 
# do the consistency check if everything is in

for(counter01 in 1:nrow(GoMF_result))
{
  theGoMF=as.vector(GoMF_result[["MF_Name"]][counter01])
  MF_david_ref_2009_tidy %>%
    dplyr::filter(MF_Name==theGoMF)->temp
  
  m2<-as.vector(temp$ID)
  m1<-strsplit(x=as.vector(GoMF_result$Genes[counter01]),split = ",")
  m1<-sub(" ","",noquote(m1[[1]]))
  m1<-paste0(tolower(substr(start = 1,stop = 3,x = m1)),substr(start = 4,stop = 4,x = m1))
  print(all(m1 %in% m2))
}
###*****************************


###*****************************
# Find go annotations (MF) // 
# select go annotation (MF) related genes 
# // find the genes in our data related with Go annotation//
GoMF_result ->GoMF_result_short
GoMF_list=as.vector(GoMF_result_short$MF_Name)
inputGenes=unique(as.vector(GoMF_input_df$gene_name))
GoMF_input_df %>%
  dplyr::select(ID=gene_name,padj_gene=padj, log2=log2FoldChange, signChange)%>%
  dplyr::mutate(score_gene=-signChange*log10(padj_gene))->GoMF_input_df_narrow

GoMF_result_short %>%
  dplyr::select(MF_Name, FDR_GoMF=FDR)%>%
  dplyr::mutate(FDR_GoMF=FDR_GoMF/100)->GoMF_result_narrow

MF_david_ref_2009_tidy %>%
  dplyr::select(ID, MF_Name)%>%
  dplyr::filter(MF_Name %in% GoMF_list)%>%
  dplyr::group_by(MF_Name)%>%
  dplyr::summarise(numGenesInCat=length(MF_Name))->numGenesInCat_df


MF_david_ref_2009_tidy %>%
  dplyr::select(ID, MF_Name)%>%
  dplyr::filter(MF_Name %in% GoMF_list,
                ID %in% inputGenes) %>%
  dplyr::left_join(.,GoMF_input_df_narrow) %>%
  dplyr::left_join(.,GoMF_result_narrow) %>%
  dplyr::left_join(.,numGenesInCat_df) %>%
  dplyr::distinct()%>%
  dplyr::filter(!is.na(padj_gene))%>%
  dplyr::group_by(MF_Name)%>%
  dplyr::arrange(desc(score_gene))->selectedDf



# add limitations for figures
selectedDf %>%
  dplyr::filter(padj_gene<0.05,FDR_GoMF<0.05) %>%
  dplyr::mutate(abs_score=abs(score_gene))%>%
  dplyr::group_by(MF_Name,signChange)%>%
  dplyr::arrange(abs_score)%>%
  dplyr::mutate(rank=signChange*seq(1,n()))%>%
  dplyr::group_by(MF_Name)%>%
  dplyr::mutate(numSigP=(max(rank)+abs(max(rank)))/2,
                numSigN=abs(min(rank)))%>%
  dplyr::group_by(ID,MF_Name)%>%
  dplyr::mutate(MF_Name_long=paste0(sub(".*~","",MF_Name),
                                    "\n padj:",
                                    sprintf("%.5f", FDR_GoMF),
                                    " N( -",numSigN,"/ +",numSigP,"/ ",numGenesInCat,")"))%>%
  dplyr::mutate(MF_Name_short=paste0(sub(".*~","",MF_Name)))%>%
  dplyr::group_by(MF_Name)%>%
  dplyr::arrange(desc(score_gene)) -> selectedDf

selectedDf %>%
  dplyr::group_by(MF_Name_long)%>%
  dplyr::summarise(FDR_GoMF=unique(FDR_GoMF))%>%
  dplyr::arrange(FDR_GoMF)->summary_df

as.vector(summary_df$MF_Name_long)


selectedDf$MF_Name_long <- factor(selectedDf$MF_Name_long, 
                                  levels = rev(as.vector(summary_df$MF_Name_long)))
###*****************************


###*****************************
# Generate simple Data Frame
# Additional Parameters
maxPathway=5
maxGene=10

if(length(unique(as.vector(selectedDf$FDR_GoMF)))<maxPathway)
{maxPathway=length(unique(as.vector(selectedDf$FDR_GoMF)))}

FDR_GoMFTopn=sort(unique(as.vector(selectedDf$FDR_GoMF)))[maxPathway]
selectedDf %>%
  dplyr::group_by()%>%
  dplyr::filter(FDR_GoMF<=FDR_GoMFTopn) %>%
  dplyr::group_by(MF_Name) %>%
  dplyr::top_n(n=maxGene, wt = padj_gene)%>%
  dplyr::group_by(MF_Name,signChange)%>%
  dplyr::arrange(abs_score)%>%
  dplyr::mutate(rank=signChange*seq(1,n()))->selectedDf_simp

selectedDf_simp %>%
  dplyr::group_by(MF_Name_short)%>%
  dplyr::summarise(FDR_GoMF=unique(FDR_GoMF))%>%
  dplyr::arrange(FDR_GoMF)->summary_df_simp

as.vector(summary_df_simp$MF_Name_short)


selectedDf_simp$MF_Name_short <- factor(selectedDf_simp$MF_Name_short, 
                                        levels = rev(as.vector(summary_df_simp$MF_Name_short)))
###*****************************


###*****************************
# Generate Figures 
# a) Complex figure
scaleHigh=max(abs(selectedDf$score_gene))
scaleMid=0
scaleLow=-max(abs(selectedDf$score_gene))

fig01=ggplot( selectedDf, aes( x=rank,y=MF_Name_long)) +
  geom_tile(aes(fill=score_gene))+
  scale_fill_gradientn(colours=c("Blue","Grey50","Red"),
                       values=rescale(c(scaleLow,scaleMid,scaleHigh)),
                       limits=c(scaleLow,scaleHigh),
                       guide = guide_colorbar(title = "-sign(cor)*P_log10"))+
  geom_text(aes(label=ID),size=3, colour="White", fontface="bold")+
  theme_bw()+
  scale_x_continuous(breaks=min(selectedDf$rank):max(selectedDf$rank))+
  ggtitle(paste0(objectName,"_mf"))+
  theme(axis.line.y = element_blank(),
        legend.position="bottom",
        axis.title.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())

print(fig01)

# # b) Simple Figure
# scaleHigh_simp=max(abs(selectedDf_simp$score_gene))
# scaleMid_simp=0
# scaleLow_simp=-max(abs(selectedDf_simp$score_gene))
# 
# fig02=ggplot( selectedDf_simp, aes( x=rank,y=MF_Name_short)) +
#   geom_tile(aes(fill=score_gene))+
#   scale_fill_gradientn(colours=c("Blue","Grey50","Red"),
#                        values=rescale(c(scaleLow_simp,scaleMid_simp,scaleHigh_simp)),
#                        limits=c(scaleLow_simp,scaleHigh_simp),
#                        guide = guide_colorbar(title = "-sign(cor)*P_log10",barwidth = 12))+
#   geom_text(aes(label=ID),size=3, colour="White", fontface="bold")+
#   theme_bw()+
#   scale_x_continuous(breaks=min(selectedDf_simp$rank):max(selectedDf_simp$rank))+
#   theme(axis.line.y = element_blank(),
#         legend.position="bottom",
#         axis.title.y = element_blank(),
#         panel.grid.minor=element_blank(),
#         panel.grid.major.x=element_blank())
# 
# print(fig02)

# c) simple figure with geom_point
scaleHigh_simp=max(abs(selectedDf_simp$log2))
scaleMid_simp=0
scaleLow_simp=-max(abs(selectedDf_simp$log2))

minimumFold=min(selectedDf_simp$log2)
if(minimumFold>-1){minimumFold=-1}
maximumFold=max(selectedDf_simp$log2)
if(maximumFold<1){maximumFold=1}

fig03=ggplot(selectedDf_simp, aes( x=log2,y=MF_Name_short)) +
  geom_point(aes(colour = log2),size=2.5)+
  geom_vline(xintercept = c(log2(1/2),log2(2)), colour="orange", linetype = "longdash")+
  geom_vline(xintercept = c(log2(1)), colour="black", linetype = "longdash")+
  geom_text_repel(aes(label=ID),size=5, colour="Black", fontface="bold")+
  scale_colour_gradientn(colours=c("Blue","Grey50","Red"),
                         values=rescale(c(scaleLow_simp,scaleMid_simp,scaleHigh_simp)),
                         limits=c(scaleLow_simp,scaleHigh_simp),
                         guide = guide_colorbar(title = "log2FoldChange",barwidth = 12))+
  theme_bw()+
  scale_x_continuous(breaks=seq(floor(minimumFold),ceiling(maximumFold)))+
  xlab("fold change")+
  theme(axis.line.y = element_blank(),
        legend.position="bottom",
        axis.title.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())

print(fig03)

GoMF_input_df %>%
  dplyr::mutate(score_gene=-signChange*log10(padj))->GoMF_input_df

fig04=ggplot2::ggplot(GoMF_input_df, aes( x=log2FoldChange,y=score_gene))+
  geom_point()

print(fig04)
###*****************************


###*****************************
# Save Files
selectedDf<-cbind(selectedDf, 
                  unique(GoMF_input_df[,c("pick_data","growthPhase","test_for","vs")]),
                  df_category="MF")
write.csv(x = selectedDf, file = paste0("../d_results/",objectName,"_mf.csv"))
###*****************************


###*****************************
# Save Figures

# Detailed Figure
colWidth=ifelse(max(selectedDf$rank)-min(selectedDf$rank)+1<16, 
                16, 
                max(selectedDf$rank)-min(selectedDf$rank))
rowWidth=ifelse(nrow(summary_df)*1<3,3,nrow(summary_df)*1)
cowplot::save_plot(filename = paste0("../d_figures/",objectName,"_mf.pdf"),
                   plot = fig01,
                   base_height = rowWidth,
                   base_width = colWidth,
                   limitsize = FALSE)

# Save simple figure
rowWidth=ifelse(nrow(summary_df_simp)*1<3,3,nrow(summary_df_simp)*1)

cowplot::save_plot(filename = paste0("../d_figures/simple",objectName,"_mf.pdf"),
                   plot = fig03,
                   base_height = rowWidth,
                   ncol=2,
                   limitsize = FALSE)






