# Play with results of DAVID Gene Ontology Results

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
require("tidyr")
###*****************************

###*****************************
# LOAD DATA
dataTypeChoice="mrna" # "mrna" or "protein"
dataTimeChoice="stationary" # "stationary" or "exponential"
setNameChoice="Ca" # "Mg", "Na", "Ca"

if(dataTypeChoice=="mrna") {dataType="RNA"}
if(dataTypeChoice=="protein") {dataType="Protein"}

if(dataTimeChoice=="exponential") {dataTime="Exp"}
if(dataTimeChoice=="stationary") {dataTime="Sta"}

if(setNameChoice=="Mg") {setName="Mg"}
if(setNameChoice=="Na") {setName="Na"}
if(setNameChoice=="Ca") {setName="Ca"}

fileName=paste0(dataType,dataTime,"_",setName ,"_comp")
davidOutput<-read.table(file = paste0("../c_results/david_results/",fileName,".txt"),
                        fill = TRUE,
                        sep = "\t",
                        header = TRUE)

# #****************************
# # Test
# davidOutput<-read.table(file = paste0("../c_results/david_results/","test",".txt"),
#                         fill = TRUE,
#                         sep = "\t",
#                         header = TRUE)
# #****************************

davidOutput=as.data.frame(davidOutput)
colnames(davidOutput)[which(names(davidOutput) == "X.")] <- "percentage"

davidOutput %>%  
  tidyr::separate(col = Term,into = c("pathway.Number","pathway.Name"),sep = "\\:") ->davidOutput

davidOutput %>%
  dplyr::filter(grepl("eco",pathway.Number)) %>%
  dplyr::mutate(FDR_umut=p.adjust(p = PValue, method = "fdr"))->davidOutputFiltered

davidOutputFiltered %>%
  dplyr::filter(FDR_umut<0.05)->davidOutputFiltered

assign(x = fileName,value = davidOutputFiltered)


save(file =paste0("../c_results/david_results/",fileName,"_filtered.rData"),
     list = fileName)
write.csv(file =paste0("../c_results/david_results/",fileName,"_filtered.csv"),
          x = davidOutputFiltered, quote = TRUE, row.names = FALSE)
###*****************************