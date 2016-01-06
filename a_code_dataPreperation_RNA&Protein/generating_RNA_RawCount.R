# GENERATION OF RAW RNA COUNT MATRICIES AND META_RNA FILE 
# this is the main start file for the pipe line for machine learning analysis.
# Outputs of the file will be several rna data frames in form of csv files.
# The data frames for
# a) rna (includes all types of RNA's i.e mRNA, tRNA, RF do not include depleted rRNA and meta information)
# b) mRNA
# c) tRNA
# d) tRNA
# e) meta file for representing information about each column
# f) a whole rna matix includes rRNA, tRNA, RF, mRNA, metaRNA
# The code also generates a file for representing overal properties of samples
#***************************************


#***************************************
# Initial Command to Reset the System & set working directory
rm(list = ls())
if (is.integer(dev.list())){dev.off()}
cat("\014")

# Set Working Directory
setwd('/Users/umut/GitHub/AG3C_Analyze/Code') # mac computer
#***************************************

#***************************************
#Load Libraries and functions
require(dplyr)
require("grid")
require("gridExtra")
library("ggplot2")

#order function
source("fun_order.R")
#***************************************



#***************************************
# load and combine AG3C_RNA_data to produce RNA matrix
#locationOfRNADataSets="../RNA_reads_02_02_2015/"
locationOfRNADataSets="../RNA_reads_03_04_2015/"
RNA_FileList=dir(path=locationOfRNADataSets)

###---REMOVE 83---###
RNA_FileList<-RNA_FileList[grepl("MURI_83",RNA_FileList)==FALSE]
#***

RawRNA_FileList=grep("raw_rna_count", RNA_FileList, value = TRUE)
RawRNA_FileList<-RawRNA_FileList[-grep("ND", RawRNA_FileList)] #get rid of not depleted 
RawRNA_colNames=gsub("_raw_rna_count.txt","",RawRNA_FileList)
RawRNA_colNames=sprintf("MURI_%03d",
                             as.numeric(gsub("MURI_","",RawRNA_colNames)))

for(counter01 in 1: length(RawRNA_FileList)){
  print(counter01)
  
  temp<-read.table(file=paste0(locationOfRNADataSets,RawRNA_FileList[counter01]),
                   header=TRUE)
  
  colnames(temp)<-c("gene_ID",RawRNA_colNames[counter01])
  if(counter01==1){
    rnaMatrix=temp
  }
  else{
    rnaMatrix=left_join(rnaMatrix,temp)
  }
  
}

remove(temp)
rnaMatrix<-rnaMatrix[sort(colnames(rnaMatrix))]

# Generate geneType Column
rnaMatrix %>% 
  mutate(gene_Type = ifelse(grepl("_r",rnaMatrix[["gene_ID"]]) , "rRNA", NA),
         gene_Type = ifelse(grepl("_t",rnaMatrix[["gene_ID"]]) , "tRNA", gene_Type),
         gene_Type = ifelse(grepl("RF",rnaMatrix[["gene_ID"]]) , "RF", gene_Type),
         gene_Type = ifelse(grepl("_0",rnaMatrix[["gene_ID"]]) , "mRNA", gene_Type),
         gene_Type = ifelse(is.na(gene_Type) , "meta", gene_Type)
         )->rnaMatrix

newOrder<-moveme(names(rnaMatrix), "gene_Type first")
rnaMatrix<-rnaMatrix[newOrder]
rnaMatrix<-dplyr::group_by(rnaMatrix,gene_Type)
# ******************************


# ******************************
# Generating Sub Matrices
rnaMatrix_meta<-dplyr::filter(rnaMatrix,gene_Type=="meta")
rnaMatrix_rRNA<-dplyr::filter(rnaMatrix,gene_Type=="rRNA")
rnaMatrix_mRNA<-dplyr::filter(rnaMatrix,gene_Type=="mRNA")
rnaMatrix_tRNA<-dplyr::filter(rnaMatrix,gene_Type=="tRNA")
rnaMatrix_RNA<-dplyr::filter(rnaMatrix,gene_Type!="rRNA" & gene_Type!="meta")
# ******************************



# *******************************

savedFilename=paste0("../Processed_RNA/","rnaRawData.Rda")
savingList=c("rnaMatrix",
             "rnaMatrix_meta",
             "rnaMatrix_mRNA",
             "rnaMatrix_rRNA",
             "rnaMatrix_tRNA",
             "rnaMatrix_RNA")
save(list=savingList,file=savedFilename)



# Save as csv files
savedFilename=paste0("../Processed_RNA/","rnaMatrix.csv")
write.csv(rnaMatrix, file = savedFilename, row.names = FALSE)

savedFilename=paste0("../Processed_RNA/","rnaMatrix_meta.csv")
write.csv(rnaMatrix_meta, file = savedFilename, row.names = FALSE)

savedFilename=paste0("../Processed_RNA/","rnaMatrix_mRNA.csv")
write.csv(rnaMatrix_mRNA, file = savedFilename, row.names = FALSE)

savedFilename=paste0("../Processed_RNA/","rnaMatrix_rRNA.csv")
write.csv(rnaMatrix_rRNA, file = savedFilename, row.names = FALSE)

savedFilename=paste0("../Processed_RNA/","rnaMatrix_tRNA.csv")
write.csv(rnaMatrix_tRNA, file = savedFilename, row.names = FALSE)

savedFilename=paste0("../Processed_RNA/","rnaMatrix_RNA.csv")
write.csv(rnaMatrix_RNA, file = savedFilename, row.names = FALSE)

## Drawing Distribution Graph
rnaMatrix_RNA %>%
  group_by() %>%
  select(-gene_ID, -gene_Type) ->p
as.data.frame(colSums(p))->p2
as.data.frame(log2(colSums(p)))->p
colnames(p)<-"log2amount"

fig01<-ggplot(p,aes_string(x = "log2amount")) +
  geom_line(aes(y=..density..), stat="density") +
  xlab("log2 amounts")+
  ggtitle("log2 amount density")+
  theme_classic()+
  xlim(0,30) +
  theme( axis.text.x=element_text(size=16),
         axis.text.y=element_text(size=16),
         axis.title.x=element_text(size=16),
         axis.title.y=element_text(size=16),
         legend.title=element_text(size=14),
         legend.text=element_text(size=14) )
plot(fig01)

ggsave(fig01, 
       filename = paste0("../figures/log2AmountDensity.pdf"), 
       width = 8, 
       height = 6, 
       units = "in", 
       useDingbats=FALSE)
