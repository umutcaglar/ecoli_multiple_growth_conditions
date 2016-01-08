# Comparing dictionaries

# Comparing genes between similar conditions.

#****************************
# Initial Command to Reset the System
rm(list = ls())
if (is.integer(dev.list())){dev.off()}
cat("\014")

# Set Working Directory
setwd('/Users/umut/GitHub/AG3C_Analyze/generateDictionary/') # mac computer
#****************************


#****************************
# Loading Libraries
require("dplyr")
require("tidyr")
#****************************


#****************************
# load 2 dictionaries
dictionary_barrick=read.csv(paste0("./nameDictionary_RNA_barrick.csv")) # gene name dictionary barrick
dictionary_internet=read.csv(paste0("./nameDictionary_RNA.internet.csv")) # gene name dictionary barrick
#****************************


#****************************
# make dictionaries unique
dictionary_barrick_unique=unique(dictionary_barrick)
dictionary_internet_unique=unique(dictionary_internet)
#****************************


#****************************
# Fix barrick unique (barick includes one empty gene id)
row.names(dictionary_barrick_unique) <- NULL 
dictionary_barrick_unique[dictionary_barrick_unique$gene_id=="",] # it is line 927
dictionary_barrick_unique<-dictionary_barrick_unique[-which(dictionary_barrick_unique$gene_id==""),] # remove it
row.names(dictionary_barrick_unique) <- NULL 
#****************************


#****************************
# Fix internet unique 
#(internet includes a row gene_id=NIL //gene_name=INFERRED_GENE GENE)
row.names(dictionary_internet_unique) <- NULL 
dictionary_internet_unique[dictionary_internet_unique$gene_id=="NIL",] # it is line 4360
dictionary_internet_unique<-dictionary_internet_unique[-which(dictionary_internet_unique$gene_id=="NIL"),] # remove it
row.names(dictionary_internet_unique) <- NULL 

# internet gene_names column contains silly "|" character
dictionary_internet_unique$gene_name<-gsub("\\|","",as.vector(dictionary_internet_unique$gene_name))
#****************************



#****************************
# n rows in unique dictionaries
nrow(dictionary_barrick_unique) #4485
nrow(dictionary_internet_unique) #4365
#****************************


#****************************
# what are unique gene id's
barrick_temp=as.vector(dictionary_barrick_unique$gene_id)
internet_temp=as.vector(dictionary_internet_unique$gene_id)
unique_gene_id_barrick=setdiff(barrick_temp, internet_temp) # many
unique_gene_id_internet=setdiff(internet_temp, barrick_temp) # "ECB_00830" "ECB_04280"
#****************************


#****************************
# Problematic cases in barrick
list_dup_barrick=as.vector(dictionary_barrick_unique$gene_name[duplicated(dictionary_barrick_unique $gene_name)])
dictionary_barrick_unique %>%
  dplyr::filter(gene_name %in% list_dup_barrick) %>%
  dplyr::arrange(gene_id)-> problematicBarrick
#****************************


#****************************
# Problematic cases in internet
list_dup_internet=as.vector(dictionary_internet_unique$gene_name[duplicated(dictionary_internet_unique $gene_name)])
dictionary_internet_unique %>%
  dplyr::filter(gene_name %in% list_dup_internet) %>%
  dplyr::arrange(gene_id)-> problematicInternet
#****************************


#****************************
# change colnames
colnames(dictionary_barrick_unique) <- c("gene_id","gene_name_barrick")
colnames(dictionary_internet_unique) <- c("gene_id","gene_name_internet")
#****************************


#****************************
# Join Data Frames
dictionary_join=merge(x = dictionary_barrick_unique, 
                      y = dictionary_internet_unique, 
                      by = "gene_id", all = TRUE)
#****************************


#****************************
# Find Nonmaching names
dictionary_join %>%
  dplyr::filter(as.vector(gene_name_barrick) != as.vector(gene_name_internet))->different_names
#****************************


#****************************
# correct barricks dictionary with internet one
dictionary_join$gene_name_barrick <- as.vector(dictionary_join$gene_name_barrick )
dictionary_join %>%
  dplyr::mutate(gene_name_barrickC=ifelse(gene_name_barrick %in% list_dup_barrick & !is.na(gene_name_internet ), 
                                          gene_name_internet, 
                                          gene_name_barrick)) %>%
  dplyr::mutate(gene_name_barrickC=ifelse(is.na(gene_name_barrick), gene_name_internet, gene_name_barrickC))-> dictionary_join
#****************************


#****************************
# Problematic cases in correction
list_dup_barrickC=as.vector(dictionary_join$gene_name_barrickC[duplicated(dictionary_join$gene_name_barrickC)])
dictionary_join %>%
  dplyr::filter(gene_name_barrickC %in% list_dup_barrickC) %>%
  dplyr::arrange(gene_id)-> problematicJoin
#****************************


#****************************
# Hand correct
dictionary_join[which(dictionary_join$gene_id=="ECB_03537"),]["gene_name_barrickC"]="yeeV"
dictionary_join[which(dictionary_join$gene_id=="ECB_04146"),]["gene_name_barrickC"]="yis1b"
#****************************


#****************************
# What is changed
dictionary_join %>% 
  dplyr::filter(gene_name_barrick %in% list_dup_barrick) %>%
  dplyr::arrange(gene_name_barrick)->changedLines
#****************************


#****************************
dictionary_join %>% 
  dplyr::select(gene_id,gene_name=gene_name_barrickC)->outPut
#****************************


#****************************
write.csv(x = outPut, 
          file = "nameDictionaryCombined_RNA.csv", 
          quote = F, 
          row.names = F, 
          col.names = T)

