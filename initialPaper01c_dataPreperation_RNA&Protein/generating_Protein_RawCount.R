#***************************
# Initial Command to Reset the System
rm(list = ls())
if (is.integer(dev.list())){dev.off()}
cat("\014")
#***************************


#***************************
# Set Working Directory
setwd('/Users/umut/GitHub/AG3C_Analyze/Code') # mac computer
#***************************


#***************************
#Load Libraries
require(dplyr)
#***************************

#***************************
#order function
source("fun_order.R")
#***************************

#***************************
# load AG3C_Protein_data
# Data came with seperate files 
# from these fileswe generate individual files where each one represents a replicate

#***************************
# NaCl stress
locationOfProteinDataSets="../Protein_reads_10_28_2015/"
data<-read.csv(file=paste0(locationOfProteinDataSets,"Protein_Counts_NaCl.csv"),header=TRUE,  fill = TRUE)
#data[is.na(data)]<-0
data %>% group_by(gene_id) -> data


# MURI stands for Multidisciplinary University Research Initiative
for(counter01 in 69:84){
  dataNumber=counter01;
  colNumber=grep(paste0(dataNumber),names(data))
  data %>% dplyr::select(colNumber) -> temp
  colnames(temp) <-c( "gene_id" , "count" )
  temp$count=as.numeric(temp$count)
  fileName=paste0("MURI_",sprintf("%03d", dataNumber),"_0_raw_protein_count.txt")
  write.table(temp,file=paste0(locationOfProteinDataSets,fileName),sep="\t",row.names=FALSE,quote = F)
  remove(temp)
}
#***************************


#*************************
# Glucose time course

data<-read.table(file=paste0(locationOfProteinDataSets,"Protein_Counts_GlucoseTimeCourse.csv"),header=TRUE,  fill = TRUE,sep =",",quote = "\"'")
colnames(data)[1]<-"gene_id"
data[[1]]=paste0("Y",sub("*.Y","",as.vector(data[[1]])))
data[[ncol(data)]]=as.numeric(data[[ncol(data)]])
#data[is.na(data)]<-0
data %>% group_by(gene_id) -> data

# generate names
nameList=as.numeric(gsub("sample","",gsub("_.*","",names(data))))
#nameList <- nameList[!is.na(nameList)]


for(counter01 in 2:length(data)){
  colNumber=grep(paste0(nameList[counter01],"_"),names(data))
  data %>% dplyr::select(colNumber) -> temp
  colnames(temp) <-c( "gene_id" , "count" )
  temp$count=as.numeric(temp$count)
  fileName=paste0("MURI_",sprintf("%03d", nameList[counter01]),"_0_raw_protein_count.txt")
  write.table(temp,file=paste0(locationOfProteinDataSets,fileName),sep="\t",row.names=FALSE,quote = F,)
  remove(temp)
}
#***************************


#*************************
# Gluconate time course
data<-read.table(file=paste0(locationOfProteinDataSets,"Protein_Counts_GluconateTimeCourse.csv"),header=F,  fill = TRUE, sep = ",",skip=1)
nameList<-c("gene_id","034_0","035_0","036_0","037_0",'038_0',"039_0","040_0","041_0",
                  "042_0","042_1","043_0","043_1","044_0","045_0","046_0","053_0","053_1","054_0","055_0","055_1","056_0")
colnames(data)<-nameList
#data[is.na(data)]<-0
data %>% group_by(gene_id) -> data

for(counter01 in 2:length(data)){
  data %>% dplyr::select(counter01) -> temp
  colnames(temp) <-c( "gene_id" , "count" )
  temp$count=as.numeric(temp$count)
  fileName=paste0("MURI_",nameList[counter01],"_raw_protein_count.txt")
  write.table(temp,file=paste0(locationOfProteinDataSets,fileName),sep="\t",row.names=FALSE,quote = F)
  remove(temp)
}
#***************************


#***************************
# Mg Time Course
data<-read.table(file=paste0(locationOfProteinDataSets,"Protein_Counts_Mg.csv"),header=T,  fill = TRUE, sep = ",",skip=0)
nameList<-colnames(data)
nameList=nameList[2:length(nameList)]
nameList=sub(pattern = "sample_",replacement = "",x = nameList)
nameList=c("gene_id",nameList)
colnames(data)<-nameList
#data[is.na(data)]<-0
data %>% group_by(gene_id) -> data

for(counter01 in 2:length(data)){
  print(nameList[counter01])
  data %>% dplyr::select(counter01) -> temp
  colnames(temp) <-c( "gene_id" , "count" )
  temp$count=as.numeric(temp$count)
  fileName=paste0("MURI_",nameList[counter01],"_0_raw_protein_count.txt")
  write.table(temp,file=paste0(locationOfProteinDataSets,fileName),sep="\t",row.names=FALSE,quote = F)
  remove(temp)
}
#***************************


#*********************Q******************
# load and combine AG3C_Protein_data to produce RNA matrix
Protein_FileList=dir(path=locationOfProteinDataSets)

RawProtein_FileList=grep("raw_protein_count", Protein_FileList, value = TRUE)
RawProtein_colNames=gsub("_raw_protein_count.txt","",RawProtein_FileList)


for(counter01 in 1: length(RawProtein_FileList)){
  print(counter01)
  
  temp<-read.table(file=paste0(locationOfProteinDataSets,RawProtein_FileList[counter01]),header=TRUE)
  row.names(temp)<-NULL
  colnames(temp)<-c("gene_id",RawProtein_colNames[counter01])
  temp[[2]]<-as.numeric(temp[[2]])
  
  # replace zeros with NA in the imported data
  temp2=temp[[2]]; temp2[is.na(temp2)]<-0;temp[[2]]=temp2;
  if(counter01==1){
    proteinMatrix=temp
  }
  else{
    proteinMatrix=full_join(proteinMatrix,temp)
  }
  
}

remove(temp)
proteinMatrix<-proteinMatrix[sort(colnames(proteinMatrix))]

# remove the MURI_082_0 from protein data
proteinMatrix %>% dplyr::select(-MURI_082_0)->proteinMatrix

# put the proteins in to order and remove contaminants
proteinMatrix %>%
  dplyr::arrange(gene_id)%>%
  dplyr::filter(!grepl('CON', gene_id))->proteinMatrix

# generate a protein matrix without NA
proteinMatrix_wo_NA <- proteinMatrix[rowSums(is.na(proteinMatrix)) == 0,]

# Save as csv files
savedFilename=paste0("../Processed_Protein/","proteinMatrix.csv")
write.csv(proteinMatrix, file = savedFilename, row.names = FALSE)

savedFilename=paste0("../Processed_Protein/","proteinMatrix_wo_NA.csv")
write.csv(proteinMatrix_wo_NA, file = savedFilename, row.names = FALSE)