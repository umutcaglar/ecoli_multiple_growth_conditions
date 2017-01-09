# Comparing the names 


###*****************************
# INITIAL COMMANDS TO RESET THE SYSTEM
rm(list = ls())
if (is.integer(dev.list())){dev.off()}
cat("\014")
set.seed(14159)

###*****************************


###*****************************
# Set Working Directory
# One needs to arrange the correct pathway if this is not umut's computer ;)
if(as.vector(Sys.info()["effective_user"]=="umut"))
{setwd(paste0("/Users/umut/GitHub/ecoli_multiple_growth_conditions/generateDictionary/"))} # mac computer
###*****************************


###*****************************
# REQUIRED LIBRARIES
# Data tracking
require("magrittr")
require("DESeq2")
require("dplyr")
require("tidyr")

# Graphing
require("ggplot2")
require("cowplot")
###*****************************


###*****************************
# Install data
mrna_protein<-read.csv(file = "measured_rna&Proteins.csv")
dictionary=read.csv(file = "../generateDictionary/nameDictionary_RNA&Protein.csv")
###*****************************


###*****************************
# Divide rna - protein
mrna_protein %>% dplyr::filter(type=="RNA")->rna_names
colnames(rna_names)[2]<-"mRNA_ID"
mrna_protein %>% dplyr::filter(type=="Protein")->protein_names
colnames(protein_names)[2]<-"Protein_id"
###*****************************


###*****************************
#Combine each with dictionary
dplyr::left_join(rna_names, dictionary)->rna_names
dplyr::left_join(protein_names, dictionary)->protein_names
###*****************************


###*****************************
# find non intersecting ones
protein_names %>%
  dplyr::group_by(mRNA_ID) %>%
  dplyr::mutate(length=n()) %>%
  dplyr::filter(length!=1)->non_intersecting_protein_names

rna_names %>%
  dplyr::group_by(Protein_id) %>%
  dplyr::mutate(length=n()) %>%
  dplyr::filter(length!=1)->non_intersecting_rna_names
###*****************************


###*****************************
# find non intersecting ones
protein_names %>%
  dplyr::group_by(mRNA_ID) %>%
  dplyr::mutate(length=n()) %>%
  dplyr::filter(length==1)->intersecting_protein_names

rna_names %>%
  dplyr::group_by(Protein_id) %>%
  dplyr::mutate(length=n()) %>%
  dplyr::filter(length==1)->intersecting_rna_names

length(as.vector(intersecting_protein_names$mRNA_ID))
length(as.vector(intersecting_protein_names$mRNA_ID))
setdiff(as.vector(intersecting_protein_names$mRNA_ID),as.vector(intersecting_protein_names$mRNA_ID))
###*****************************


###*****************************
# find unique ones
setdiff(intersecting_rna_names$mRNA_ID, intersecting_protein_names$mRNA_ID)->unique_in_rna
setdiff(intersecting_protein_names$mRNA_ID, intersecting_rna_names$mRNA_ID)->unique_in_protein

intersecting_rna_names %>% 
  dplyr::filter(mRNA_ID %in% unique_in_rna) %>% 
  dplyr::select(mRNA_ID, type, Protein_id, gene_name)->unique_in_rnaDF

intersecting_protein_names %>% 
  dplyr::filter(mRNA_ID %in% unique_in_protein) %>% 
  dplyr::select(mRNA_ID, type, Protein_id, gene_name)->unique_in_proteinDF
###*****************************


###*****************************
# List of problematic genes
oldNewDifference<-c(as.vector(non_intersecting_rna_names$mRNA_ID),
                     as.vector(non_intersecting_protein_names$Protein_id),
                     as.vector(unique_in_rnaDF$mRNA_ID),
                     as.vector(unique_in_rnaDF$Protein_id),
                     as.vector(unique_in_rnaDF$gene_name),
                     as.vector(unique_in_proteinDF$mRNA_ID),
                     as.vector(unique_in_proteinDF$Protein_id),
                     as.vector(unique_in_proteinDF$gene_name))
###*****************************


###*****************************
write.csv(x = oldNewDifference, file = "oldNewDifference.csv")
###*****************************