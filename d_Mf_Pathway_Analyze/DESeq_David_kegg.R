# Analayze Kegg Pathways DESeq + DAVID

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
source("../a_code_dataPreperation_RNA&Protein/data_naming_functions.R")
###*****************************


###*****************************
# Download the DAVID input and output
dataName=name_data(initialValue="genes_P0.05Fold2", # can be "genes0.05", "genes_P0.05Fold2"
                   dataType = "mrna", # can be "rna", "mrna", "protein", "protein_wo_NA"
                   badDataSet = "set00", # can be "set00",set01","set02", "set03"
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
                   carbonSourceVector = "SYAN", # can be any sub combination of "SYAN"
                   MgLevelVector = c("allMg"), # can be "lowMg","baseMg","highMg" // "allMg"
                   NaLevelVector = c("allNa"), # can be "baseNa","highNa" // "allNa"
                   growthPhaseVector = c("exponential"), # can be "exponential","stationary","late_stationary" // "allPhase"
                   filterGenes = "noFilter", # can be "noFilter", "meanFilter", "maxFilter", "sdFilter" 
                   threshold=NA, # the threshold value for "meanFilter", "maxFilter", "sdFilter"
                   roundData=TRUE,
                   sumTechnicalReplicates=TRUE,
                   deSeqSfChoice="p1Sf", # can be "regSf", "p1Sf"
                   normalizationMethodChoice= "noNorm", # can be "vst", "rlog", "log10", "noNorm"
                   test_for = "Mg_mM_Levels")  # works only if normalizationMethodChoice == noNorm
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")

# DeSeq2 parameters
objectName_df=dataName$objectName

objectName_df$test_for=paste0("_batchNumberPLUS",gsub("^_","",objectName_df$test_for))
test_base="baseMg"
test_contrast="highMg"
objectName_df$contrast=paste0("_",test_contrast,"VS",test_base)

objectName=paste(objectName_df,collapse = "_")

kegg_input<-read.csv(file = paste0("../c_results/",objectName,".csv"),
                     header = TRUE)

objectName_df$initial="resDf"
objectName=paste(objectName_df,collapse = "_")
kegg_input_df<-read.csv(file = paste0("../c_results/",objectName,".csv"),
                        header = TRUE)

objectName_df$initial="ez_P0.05Fold2"
objectName=paste(objectName_df,collapse = "_")
kegg_result<-read.csv(file=paste0("../c_results/",objectName,"_kegg.csv"),header = TRUE)
kegg_result<-kegg_result[grep("eum",as.vector(kegg_result$Term)),]
browser()
###*****************************


###*****************************
# load reference libraries
dataType=dataName$objectName$pick_data

if(dataType %in% c("rna","mrna"))
{
  kegg_pathway_david_ref_2009_tidy<-read.table(file = paste0("../d_Mf_Pathway_Analyze/ReferenceFiles/",
                                                             "kegg_pathway_david_mrna_ref_2009_tidy.txt"),
                                               sep = "\t",header = TRUE,fill = TRUE, quote = "")
}


if(dataType %in% c("protein_wo_NA", "protein"))
{
  kegg_pathway_david_ref_2009_tidy<-read.table(paste0("../d_Mf_Pathway_Analyze/ReferenceFiles/",
                                                      "kegg_pathway_david_protein_ref_2009_tidy.txt"),
                                               sep = "\t",header = TRUE,fill = TRUE, quote = "")
}
###*****************************

browser()
###*****************************
# Find KEGG output pathways in ref 
# do the consistency check if everything is in

for(counter01 in 1:nrow(kegg_result))
{
  thePathway=as.vector(kegg_result[["Term"]][counter01])
  kegg_pathway_david_ref_2009_tidy %>%
    dplyr::filter(KEGG_PATH_Name==thePathway)->temp
  
  m2<-as.vector(temp$ID)
  m1<-strsplit(x=as.vector(kegg_result$Genes[counter01]),split = ",")
  m1<-sub(" ","",noquote(m1[[1]]))
  m1<-paste0(tolower(substr(start = 1,stop = 3,x = m1)),substr(start = 4,stop = 4,x = m1))
  print(all(m1 %in% m2))
}
###*****************************


###*****************************
# Find pathways // 
# select pathway related genes 
# // find the gens in our data related with pathway//

kegg_result ->kegg_result_short
pathway_list=as.vector(kegg_result_short$Term)
inputGenes=unique(as.vector(kegg_input_df$gene_name))
kegg_input_df %>%
  dplyr::select(ID=gene_name, padj_gene=padj, log2=log2FoldChange, signChange)%>%
  dplyr::mutate(score_gene=-signChange*log10(padj_gene))->kegg_input_df_narrow

kegg_result_short %>%
  dplyr::select(KEGG_Path=Term, FDR_KEGG_Path=FDR)%>%
  dplyr::mutate(FDR_KEGG_Path=FDR_KEGG_Path/100)->kegg_result_narrow

kegg_pathway_david_ref_2009_tidy %>%
  dplyr::select(ID, KEGG_Path=KEGG_PATH_Name)%>%
  dplyr::filter(KEGG_Path %in% pathway_list)%>%
  dplyr::group_by(KEGG_Path)%>%
  dplyr::summarise(numGenesInCat=length(KEGG_Path))->numGenesInCat_df

kegg_pathway_david_ref_2009_tidy %>%
  dplyr::select(ID, KEGG_Path=KEGG_PATH_Name)%>%
  dplyr::filter(KEGG_Path %in% pathway_list,
                ID %in% inputGenes) %>%
  dplyr::left_join(.,kegg_input_df_narrow) %>%
  dplyr::left_join(.,kegg_result_narrow) %>%
  dplyr::left_join(.,numGenesInCat_df) %>%
  dplyr::filter(!is.na(padj_gene))%>%
  dplyr::group_by(KEGG_Path)%>%
  dplyr::arrange(desc(score_gene))->selectedDf

# add limitations for figures
selectedDf %>%
  dplyr::filter(padj_gene<0.05,FDR_KEGG_Path<0.05, abs(log2)>1) %>%
  dplyr::mutate(abs_score=abs(score_gene))%>%
  dplyr::group_by(KEGG_Path,signChange)%>%
  dplyr::arrange(abs_score)%>%
  dplyr::mutate(rank=signChange*seq(1,n()))%>%
  dplyr::group_by(KEGG_Path)%>%
  dplyr::mutate(numSigP=(max(rank)+abs(max(rank)))/2,
                numSigN=abs(min(rank)))%>%
  dplyr::group_by(ID,KEGG_Path)%>%
  dplyr::mutate(KEGG_Path_long=paste0(sub(".*:","",KEGG_Path),
                                      "\n padj:",
                                      sprintf("%.5f", FDR_KEGG_Path),
                                      " N( -",numSigN,"/ +",numSigP,"/ ",numGenesInCat,")"))%>%
  dplyr::mutate(KEGG_Path_short=paste0(sub(".*:","",KEGG_Path)))%>%
  dplyr::group_by(KEGG_Path)%>%
  dplyr::arrange(desc(score_gene))-> selectedDf


selectedDf %>%
  dplyr::group_by(KEGG_Path_long)%>%
  dplyr::summarise(FDR_KEGG_Path=unique(FDR_KEGG_Path))%>%
  dplyr::arrange(FDR_KEGG_Path)->summary_df

as.vector(summary_df$KEGG_Path_long)


selectedDf$KEGG_Path_long <- factor(selectedDf$KEGG_Path_long, 
                                    levels = rev(as.vector(summary_df$KEGG_Path_long)))
###*****************************


###*****************************
# Generate simple Data Frame
# Additional Parameters
maxPathway=10
maxGene=15

if(length(unique(as.vector(selectedDf$FDR_KEGG_Path)))<maxPathway)
{maxPathway=length(unique(as.vector(selectedDf$FDR_KEGG_Path)))}

FDR_KEGG_PathTopn=sort(unique(as.vector(selectedDf$FDR_KEGG_Path)))[maxPathway]

selectedDf %>%
  dplyr::group_by()%>%
  dplyr::filter(FDR_KEGG_Path<=FDR_KEGG_PathTopn) %>%
  dplyr::group_by(KEGG_Path) %>%
  dplyr::arrange(desc(abs_score))%>%
  dplyr::top_n(n=maxGene, wt = abs_score)%>%
  dplyr::group_by(KEGG_Path,signChange)%>%
  dplyr::arrange(abs_score)%>%
  dplyr::mutate(rank=signChange*seq(1,n()))->selectedDf_simp

selectedDf_simp %>%
  dplyr::group_by(KEGG_Path_short)%>%
  dplyr::summarise(FDR_KEGG_Path=unique(FDR_KEGG_Path))%>%
  dplyr::arrange(FDR_KEGG_Path)->summary_df_simp

as.vector(summary_df_simp$KEGG_Path_short)


selectedDf_simp$KEGG_Path_short <- factor(selectedDf_simp$KEGG_Path_short, 
                                          levels = rev(as.vector(summary_df_simp$KEGG_Path_short)))
###*****************************


###*****************************
# Generate Figures 
# a) Complex figure
scaleHigh=max(abs(selectedDf$score_gene))
scaleMid=0
scaleLow=-max(abs(selectedDf$score_gene))

fig01=ggplot( selectedDf, aes( x=rank,y=KEGG_Path_long)) +
  geom_tile(aes(fill=score_gene))+
  scale_fill_gradientn(colours=c("Blue","Grey50","Red"),
                       values=rescale(c(scaleLow,scaleMid,scaleHigh)),
                       limits=c(scaleLow,scaleHigh),
                       guide = guide_colorbar(title = "-sign(cor)*P_log10"))+
  geom_text(aes(label=ID),size=3, colour="White", fontface="bold")+
  theme_bw()+
  scale_x_continuous(breaks=min(selectedDf$rank):max(selectedDf$rank))+
  ggtitle(paste0(objectName,"_kegg"))+
  theme(axis.line.y = element_blank(),
        legend.position="bottom",
        axis.title.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())

print(fig01)



# b) simple figure with geom point

minimumFold=min(selectedDf_simp$log2)
if(minimumFold>-1){minimumFold=-1}
maximumFold=max(selectedDf_simp$log2)
if(maximumFold<1){maximumFold=1}

fig03=ggplot(selectedDf_simp, aes( x=log2,y=KEGG_Path_short)) +
  geom_point(colour="blue", size=2.5)+
  geom_vline(xintercept = c(log2(1/2),log2(2)), colour="orange", linetype = "longdash")+
  geom_vline(xintercept = c(log2(1)), colour="black", linetype = "longdash")+
  geom_text_repel(aes(label=ID),size=3, colour="Black", fontface="plain")+
  theme_bw()+
  scale_x_continuous(breaks=seq(floor(minimumFold),ceiling(maximumFold)))+
  xlab("Log2 Fold Change")+
  theme(axis.line.y = element_blank(),
        legend.position="bottom",
        axis.title.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14))


###*****************************


###*****************************
# Save Files
selectedDf<-cbind(selectedDf, 
                  unique(kegg_input_df[,c("pick_data","growthPhase","test_for","vs")]),
                  df_category="kegg")
write.csv(x = selectedDf, file = paste0("../d_results/",objectName,"_kegg.csv"))
###*****************************

###*****************************
# Save Figures

# Detailed Figure
colWidth=ifelse(max(selectedDf$rank)-min(selectedDf$rank)+1<16, 
                16, 
                max(selectedDf$rank)-min(selectedDf$rank))
rowWidth=ifelse(nrow(summary_df)*1<3,3,nrow(summary_df)*1)
cowplot::save_plot(filename = paste0("../d_figures/",objectName,"_kegg.pdf"),
                   plot = fig01,
                   nrow=1.5,
                   ncol=2,
                   base_height = rowWidth,
                   base_width = colWidth,
                   limitsize = FALSE)

# Save simple figure
rowWidth=ifelse(nrow(summary_df_simp)*1<3,3,nrow(summary_df_simp)*1)

cowplot::save_plot(filename = paste0("../d_figures/simple",objectName,"_kegg.pdf"),
                   plot = fig03,
                   base_height = rowWidth,
                   ncol=3,
                   nrow=1.2,
                   limitsize = FALSE)









