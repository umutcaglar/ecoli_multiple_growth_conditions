# Organize reference files for kegg pathway and go analyze

# Aim of the code is to organize the files downladed from internet for reference.

# the first 2 files are downloaded from DAVID they are the kegg pathway annotations that david use
# they were generated in 2009 so a bit old
# the other file associates go annotations with genes.


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
              "d_Mf_Pathway_Analyze/ReferenceFiles/"))} # mac computer
###*****************************


###*****************************
# DOWNLOAD LIBRARIES
require("dplyr")
require("tidyr")
###*****************************


###*****************************
# Organize "kegg_pathway_david_mrna_ref_2009"

kegg_pathway_david_mrna_ref_2009<-as.data.frame(read.table(file = "kegg_pathway_david_mrna_ref_2009.txt",
                                                           sep = "\t",header = TRUE, quote = ""))
kegg_pathway_david_mrna_ref_2009[["KEGG_PATHWAY"]]=
  as.character(kegg_pathway_david_mrna_ref_2009[["KEGG_PATHWAY"]])

kegg_pathway_david_mrna_ref_2009 %>%
  dplyr::group_by(ID)%>%
  dplyr::mutate(KEGG_PATHWAY=gsub(", ","__",KEGG_PATHWAY))%>%
  dplyr::mutate(numCommas=length(gregexpr(",",KEGG_PATHWAY)[[1]]))->kegg_pathway_david_mrna_ref_2009

maxCommas_rna<-max(kegg_pathway_david_mrna_ref_2009[["numCommas"]])
colNameVector_rna<-sprintf("KEGG_Path_%03d",seq(1,maxCommas_rna))

kegg_pathway_david_mrna_ref_2009 %>%
  tidyr::separate(col = KEGG_PATHWAY, 
                  into = colNameVector_rna,
                  sep = ",")->kegg_pathway_david_mrna_ref_2009

kegg_pathway_david_mrna_ref_2009<-as.data.frame(apply(kegg_pathway_david_mrna_ref_2009, 
                                                      2, 
                                                      function(x) gsub("^$|^ $|\"", NA, x)))


kegg_pathway_david_mrna_ref_2009 %>%
  tidyr::gather(key="KEGG_Path",
                value = "KEGG_PATH_Name",
                ...=get(colNameVector_rna[1]):get(colNameVector_rna[maxCommas_rna])) %>%
  dplyr::arrange(Gene.Name) %>%
  dplyr::filter(!is.na(as.character(KEGG_PATH_Name)))->kegg_pathway_david_mrna_ref_2009


# Only keep "eco"
kegg_pathway_david_mrna_ref_2009 %>%
  dplyr::mutate(KEGG_PATH_Name=gsub("__",", ",KEGG_PATH_Name))%>%
  dplyr::filter(grepl("^eco",KEGG_PATH_Name)) %>%
  dplyr::arrange(ID,KEGG_Path)->kegg_pathway_david_mrna_ref_2009

kegg_pathway_david_mrna_ref_2009 %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(count = seq(1,n())) %>%
  dplyr::select(-Gene.Name)->kegg_pathway_david_mrna_ref_2009

write.table(x = kegg_pathway_david_mrna_ref_2009, 
            file = "kegg_pathway_david_mrna_ref_2009_tidy.txt",
            sep="\t",row.names = FALSE, quote = FALSE)
###*****************************



###*****************************
# Organize "kegg_pathway_david_protein_ref_2009"

kegg_pathway_david_protein_ref_2009<-
  as.data.frame(read.table(file = "kegg_pathway_david_protein_ref_2009.txt",
                                                              sep = "\t",header = TRUE,quote = ""))
kegg_pathway_david_protein_ref_2009[["KEGG_PATHWAY"]]=
  as.character(kegg_pathway_david_protein_ref_2009[["KEGG_PATHWAY"]])


kegg_pathway_david_protein_ref_2009 %>%
  dplyr::group_by(ID)%>%
  dplyr::mutate(KEGG_PATHWAY=gsub(", ","__",KEGG_PATHWAY))%>%
  dplyr::mutate(numCommas=length(gregexpr(",",KEGG_PATHWAY)[[1]]))->kegg_pathway_david_protein_ref_2009

maxCommas_protein<-max(kegg_pathway_david_protein_ref_2009[["numCommas"]])
colNameVector_protein<-sprintf("KEGG_Path_%03d",seq(1,maxCommas_protein))

kegg_pathway_david_protein_ref_2009 %>%
  tidyr::separate(col = KEGG_PATHWAY, 
                  into = colNameVector_protein,
                  sep = ",")->kegg_pathway_david_protein_ref_2009

kegg_pathway_david_protein_ref_2009<-as.data.frame(apply(kegg_pathway_david_protein_ref_2009, 
                                                      2, 
                                                      function(x) gsub("^$|^ $|\"", NA, x)))


kegg_pathway_david_protein_ref_2009 %>%
  tidyr::gather(key="KEGG_Path",
                value = "KEGG_PATH_Name",
                ...=get(colNameVector_protein[1]):get(colNameVector_protein[maxCommas_protein])) %>%
  dplyr::arrange(Gene.Name) %>%
  dplyr::filter(!is.na(as.character(KEGG_PATH_Name)))->kegg_pathway_david_protein_ref_2009


# Only keep "eco"
kegg_pathway_david_protein_ref_2009 %>%
  dplyr::mutate(KEGG_PATH_Name=gsub("__",", ",KEGG_PATH_Name))%>%
  dplyr::filter(grepl("^eco",KEGG_PATH_Name)) %>%
  dplyr::arrange(ID,KEGG_Path)->kegg_pathway_david_protein_ref_2009

kegg_pathway_david_protein_ref_2009 %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(count = seq(1,n())) %>%
  dplyr::select(-Gene.Name)->kegg_pathway_david_protein_ref_2009

write.table(x = kegg_pathway_david_protein_ref_2009, 
            file = "kegg_pathway_david_protein_ref_2009_tidy.txt",
            sep="\t",row.names = FALSE, quote = FALSE)
###*****************************

m1=unique(as.vector(kegg_pathway_david_mrna_ref_2009$ID))
m2=unique(as.vector(kegg_pathway_david_protein_ref_2009$ID))
setdiff(m1,m2)
setdiff(m2,m1)