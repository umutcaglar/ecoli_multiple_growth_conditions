# Gnerate ENTREZ_GENE_ID

# Comparing dictionaries

# Comparing genes between similar conditions.

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


#****************************
# Loading Libraries
require("dplyr")
require("tidyr")
#****************************


#****************************
# load 2 EZ dictionaries i 2 pieces
rna1_ez=read.table(file = "RNA_list1_ez.txt",header = TRUE,sep = "\t",quote="\"",fill = TRUE)
rna2_ez=read.table(file = "RNA_list2_ez.txt",header = TRUE,sep = "\t",quote="\"",fill = TRUE)
protein1_ez=read.table(file = "Protein_list1_ez.txt",header = TRUE,sep = "\t",quote="\"",fill = TRUE)
protein2_ez=read.table(file = "Protein_list2_ez.txt",header = TRUE,sep = "\t",quote="\"",fill = TRUE)
#****************************


#****************************
# Combine
dplyr::rbind_list(rna1_ez,rna2_ez)->rna_tidy_ez
dplyr::rbind_list(protein1_ez,protein2_ez)->protein_tidy_ez
#****************************


#****************************
# filter Ecoli
rna_tidy_ez%>%
  dplyr::filter(grepl("Escherichia coli str. K-12 substr. MG1655",Species))->rna_tidy_eColi_ez

protein_tidy_ez%>%
  dplyr::filter(grepl("Escherichia coli str. K-12 substr. MG1655",Species))->protein_tidy_eColi_ez
#****************************


#****************************
# make lists Untidy
rna_tidy_eColi_ez[["Species"]]=factor(rna_tidy_eColi_ez[["Species"]])
rna_tidy_eColi_ez%>%
  group_by(Species,Gene.Name)%>%
  mutate(row = 1:n()) %>%
  tidyr::spread(key = Species, value = To)->rna_eColi_ez


protein_tidy_eColi_ez[["Species"]]=factor(protein_tidy_eColi_ez[["Species"]])
protein_tidy_eColi_ez%>%
  group_by(Species,Gene.Name)%>%
  mutate(row = 1:n()) %>%
  tidyr::spread(key = Species, value = To)->protein_eColi_ez
#****************************


#****************************
# Save files
write.csv(x = protein_tidy_ez,file = "protein_tidy_ez.csv")
write.csv(x = rna_tidy_ez,file = "rna_tidy_ez.csv")

write.csv(x = protein_eColi_ez,file = "protein_eColi_ez.csv")
write.csv(x = rna_eColi_ez,file = "rna_eColi_ez.csv")

write.csv(x = protein_tidy_eColi_ez,file = "protein_tidy_eColi_ez.csv")
write.csv(x = rna_tidy_eColi_ez,file = "rna_tidy_eColi_ez.csv")
#****************************
#****************************