# This file is for plotting supplementary figures for KEGG pathways 
# that seems to be significantly changing under different experiments.

# This file is for single pathway. Next file will do this atomatically in loop


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
require("stringr")

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
# Generate filename and load files
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

# The file is the list of significantly altered genes (without the information how much they are altered)
objectName_df$initial="genes_P0.05Fold2"
objectName=paste(objectName_df,collapse = "_")
mf_input<-read.csv(file = paste0("../c_results/",objectName,".csv"),header = TRUE) 

# the file is the list of whole genes and includes information aout how much they are altered
objectName_df$initial="resDf"
objectName=paste(objectName_df,collapse = "_")
mf_input_df<-read.csv(file = paste0("../c_results/",objectName,".csv"),header = TRUE)

# The file is the output of RDAVIDWebService for KEGG pathways with given input "kegg_input"
# The only problem is it gives the ez ids for genes
objectName_df$initial="ez_P0.05Fold2"
objectName=paste(objectName_df,collapse = "_")
mf_result<-read.csv(file=paste0("../c_results/david_results/",objectName,"_mf_o.csv"),header = TRUE)

# the file makes the transition between official gene names and ez name
if(objectName_df$pick_data=="mrna")
{
  officialGeneSymbol_entrezGeneId<-read.csv(file="../generateDictionary/rna_tidy_eColi_ez.csv")
  officialGeneSymbol_entrezGeneId %>%
    dplyr::filter(Species=="Escherichia coli str. K-12 substr. MG1655")->officialGeneSymbol_entrezGeneId
}

if(objectName_df$pick_data=="protein")
{
  officialGeneSymbol_entrezGeneId<-read.csv(file="../generateDictionary/protein_tidy_eColi_ez.csv")
  officialGeneSymbol_entrezGeneId %>%
    dplyr::filter(Species=="Escherichia coli str. K-12 substr. MG1655")->officialGeneSymbol_entrezGeneId
}
###*****************************
browser()

###*****************************
# Put the MF output into a more useful format

# 1. one needs to divide the column of genes into multiple columns and turn the data into tidy format
additionalColumnNames<-sprintf("gene_%02d", 1:(max(str_count(mf_result$Genes, ","))+1))


mf_result %>%
  dplyr::rename(MF=Term, FDR_MF=FDR)%>%
  dplyr::mutate(FDR_MF=FDR_MF/100)%>%
  dplyr::mutate(MF_Short=gsub("*.*~","",MF))%>%
  tidyr::separate(col = Genes, into = additionalColumnNames, sep = ",")->mf_result_divided

mf_result_divided %>%
  tidyr::gather(key = gene_number, value = gene_name_ez, 
                dplyr::starts_with("gene_"), na.rm = TRUE)->mf_result_tidy
mf_result_tidy$gene_name_ez=gsub(" ", "", mf_result_tidy$gene_name_ez)

# 2. combine the entrez gene ides in kegg result with official gene symbols
officialGeneSymbol_entrezGeneId$To=as.character(officialGeneSymbol_entrezGeneId$To)
officialGeneSymbol_entrezGeneId %>%
  dplyr::select(gene_name_ez=To, gene_name=From)%>%
  dplyr::left_join(mf_result_tidy,.)->mf_result_tidy


# 3. Now add the information about individual genes (p.adj and log2 change) into tidy data
mf_input_df %>%
  dplyr::select(gene_name, padj_gene=padj, log2=log2FoldChange, signChange)%>%
  dplyr::mutate(score_gene=-signChange*log10(padj_gene))->mf_input_df_narrow 
# includes info abut how much each gene altered in terms of p.adj and log2 change

dplyr::left_join(mf_result_tidy,mf_input_df_narrow)->mf_result_tidy

# 4. filter the tidy data padj_gene<0.05,FDR_KEGG_Path<0.05, abs(log2)>1 and add rank
mf_result_tidy%>%
  dplyr::filter(padj_gene<0.05,FDR_MF<0.05, abs(log2)>1) %>%
  dplyr::mutate(abs_score=abs(score_gene))%>%
  dplyr::group_by(MF,signChange)%>%
  dplyr::arrange(abs_score)%>%
  dplyr::mutate(rank=signChange*seq(1,n()))%>%
  dplyr::group_by(MF)%>%
  dplyr::mutate(numSigP=(max(rank)+abs(max(rank)))/2,
                numSigN=abs(min(rank)))%>%
  dplyr::group_by(gene_name,MF)%>%
  dplyr::mutate(MF_long=paste0(sub(".*~","",MF),
                                      "\n padj:",
                                      sprintf("%.5f", FDR_MF),
                                      " N( -",numSigN,"/ +",numSigP,"/ ",Pop.Hits,")"))%>%
  dplyr::mutate(MF_short=paste0(sub(".*~","",MF)))%>%
  dplyr::group_by(MF)%>%
  dplyr::arrange(desc(score_gene))->mf_tidy_organized

# 5. order the factors for KEGG_Path_long
mf_tidy_organized %>% 
  dplyr::group_by(MF_long)%>%
  dplyr::summarise(FDR_MF=unique(FDR_MF))%>%
  dplyr::arrange(FDR_MF)->mf_summary


mf_tidy_organized$MF_long <- factor(mf_tidy_organized$MF_long, 
                                             levels = rev(as.vector(mf_summary$MF_long)))
###*****************************
browser()

###*****************************
# Generate simple Data Frame
# Additional Parameters
maxPathway=10
maxGene=15

if(length(unique(as.vector(mf_tidy_organized$FDR_MF)))<maxPathway)
{maxPathway=length(unique(as.vector(mf_tidy_organized$FDR_MF)))}

FDR_MFTopn=sort(unique(as.vector(mf_tidy_organized$FDR_MF)))[maxPathway]

mf_tidy_organized %>%
  dplyr::group_by()%>%
  dplyr::filter(FDR_MF<=FDR_MFTopn) %>%
  dplyr::group_by(MF) %>%
  dplyr::arrange(desc(abs_score))%>%
  dplyr::top_n(n=maxGene, wt = abs_score)%>%
  dplyr::group_by(MF,signChange)%>%
  dplyr::arrange(abs_score)%>%
  dplyr::mutate(rank=signChange*seq(1,n()))->mf_tidy_organized_simp

mf_tidy_organized_simp %>%
  dplyr::group_by(MF_short)%>%
  dplyr::summarise(FDR_MF=unique(FDR_MF))%>%
  dplyr::arrange(FDR_MF)->mf_organized_summary

mf_tidy_organized_simp$MF_short <- factor(mf_tidy_organized_simp$MF_short, 
                                                   levels = rev(as.vector(mf_organized_summary$MF_short)))
###*****************************
browser()

###*****************************
# simple figure with geom point

minimumFold=min(mf_tidy_organized_simp$log2)
if(minimumFold>-1){minimumFold=-1}
maximumFold=max(mf_tidy_organized_simp$log2)
if(maximumFold<1){maximumFold=1}

fig02=ggplot(mf_tidy_organized_simp, aes( x=log2,y=MF_short)) +
  geom_point(colour="blue", size=2.5)+
  geom_vline(xintercept = c(log2(1/2),log2(2)), colour="orange", linetype = "longdash")+
  geom_vline(xintercept = c(log2(1)), colour="black", linetype = "longdash")+
  geom_text_repel(aes(label=gene_name),size=3, colour="Black", fontface="plain")+
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

print(fig02)
###*****************************
browser()

###*****************************
# Save Files
mf_tidy_organized %>%
  dplyr::mutate(pick_data=unique(mf_input_df[,c("pick_data")]),
                growthPhase=unique(mf_input_df[,c("growthPhase")]),
                test_for=unique(mf_input_df[,c("test_for")]),
                base=unique(mf_input_df[,c("base")]),
                contrast=unique(mf_input_df[,c("contrast")]),
                df_category="kegg")->mf_tidy_organized


write.csv(x = mf_tidy_organized, file = paste0("../d_results/",objectName,"_mf_o.csv"))
###*****************************


###*****************************
# Save figure
rowWidth=ifelse(nrow(mf_organized_summary)*1<3,3,nrow(mf_organized_summary)*1)

cowplot::save_plot(filename = paste0("../d_figures/simple",objectName,"_mf_o.pdf"),
                   plot = fig02,
                   base_height = rowWidth,
                   ncol=3,
                   nrow=1.2,
                   limitsize = FALSE)
###*****************************
