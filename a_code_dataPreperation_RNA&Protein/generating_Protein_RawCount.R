#***************************
# Initial Command to Reset the System
rm(list = ls())
if (is.integer(dev.list())){dev.off()}
cat("\014")
#***************************


###*****************************
# Set Working Directory
# One needs to arrange the correct pathway if this is not umut's computer ;)
if(as.vector(Sys.info()["effective_user"]=="umut"))
{setwd(paste0("/Users/umut/GitHub/ecoli_multiple_growth_conditions/",
              "a_code_dataPreperation_RNA&Protein//"))} # mac computer
###*****************************


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
locationOfProteinDataSets="../Protein_reads_06_16_2016/"
data<-read.csv(file=paste0(locationOfProteinDataSets,"Protein_Counts_NaCl.csv"),
               header=TRUE,  fill = TRUE)
#data[is.na(data)]<-0
data %>% group_by(gene_id) -> data
NaCl_allzero_rows=which(rowSums((!is.na(data))+0)==1)

# MURI stands for Multidisciplinary University Research Initiative
for(counter01 in 69:84){
  dataNumber=counter01;
  colNumber=grep(paste0(dataNumber),names(data))
  data %>% dplyr::select(colNumber) -> temp
  colnames(temp) <-c( "gene_id" , "count" )
  temp$count=as.numeric(temp$count)
  fileName=paste0("MURI_",sprintf("%03d", dataNumber),"_0_raw_protein_count.txt")
  write.table(temp,file=paste0(locationOfProteinDataSets,fileName),
              sep="\t",row.names=FALSE,quote = F)
  remove(temp)
}
#***************************


#*************************
# Glucose time course

data<-read.table(file=paste0(locationOfProteinDataSets,
                             "Protein_Counts_GlucoseTimeCourse.csv"),
                 header=TRUE,  fill = TRUE,sep =",",quote = "\"'")
colnames(data)[1]<-"gene_id"
data[[1]]=paste0("Y",sub("*.Y","",as.vector(data[[1]])))
data[[ncol(data)]]=as.numeric(data[[ncol(data)]])
#data[is.na(data)]<-0
data %>% group_by(gene_id) -> data
GlucoseTC_allzero_rows=which(rowSums((!is.na(data))+0)==1)


# generate names
nameList=as.numeric(gsub("sample","",gsub("_.*","",names(data))))
#nameList <- nameList[!is.na(nameList)]


for(counter01 in 2:length(data)){
  colNumber=grep(paste0(nameList[counter01],"_"),names(data))
  data %>% dplyr::select(colNumber) -> temp
  colnames(temp) <-c( "gene_id" , "count" )
  temp$count=as.numeric(temp$count)
  fileName=paste0("MURI_",sprintf("%03d", nameList[counter01]),
                  "_0_raw_protein_count.txt")
  write.table(temp,file=paste0(locationOfProteinDataSets,fileName),
              sep="\t",row.names=FALSE,quote = F)
  remove(temp)
}
#***************************


#*************************
# Glycerol time course
data<-read.table(file=paste0(locationOfProteinDataSets,
                             "Protein_Counts_GlycerolTimeCourse.csv"),
                 header=T,  fill = TRUE, sep = ",")
nameList=sprintf("%03d",as.numeric(gsub("Sample","",
                                        gsub("_.*","",colnames(data[2:ncol(data)])))))
nameList<-paste0(nameList,"_0")
nameListLong<-paste0("sample_",nameList)
nameList<-c("gene_id",nameList)
nameListLong<-c("gene_id",nameListLong)
colnames(data)<-nameListLong
#data[is.na(data)]<-0
data %>% group_by(gene_id) -> data
GlycerolTC_allzero_rows=which(rowSums((!is.na(data))+0)==1)


for(counter01 in 2:length(data)){
  data %>% dplyr::select(counter01) -> temp
  colnames(temp) <-c( "gene_id" , "count" )
  temp$count=as.numeric(temp$count)
  fileName=paste0("MURI_",nameList[counter01],"_raw_protein_count.txt")
  write.table(temp,file=paste0(locationOfProteinDataSets,fileName),sep="\t",
              row.names=FALSE,quote = F)
  remove(temp)
}
#***************************


#***************************
# MgSO4 Data
data<-read.table(file=paste0(locationOfProteinDataSets,"Protein_Counts_Mg.csv"),
                 header=T,  fill = TRUE, sep = ",",skip=0)
nameList<-colnames(data)
nameList=nameList[2:length(nameList)]
nameList=sub(pattern = "sample_",replacement = "",x = nameList)
nameList=c("gene_id",nameList)
colnames(data)<-nameList
#data[is.na(data)]<-0
data %>% group_by(gene_id) -> data
MgSO4_allzero_rows=which(rowSums((!is.na(data))+0)==1)

for(counter01 in 2:length(data)){
  print(nameList[counter01])
  data %>% dplyr::select(counter01) -> temp
  colnames(temp) <-c( "gene_id" , "count" )
  temp$count=as.numeric(temp$count)
  fileName=paste0("MURI_",nameList[counter01],"_0_raw_protein_count.txt")
  write.table(temp,file=paste0(locationOfProteinDataSets,fileName),sep="\t",
              row.names=FALSE,quote = F)
  remove(temp)
}
#***************************


#***************************
# Gluconate data
data<-read.table(file=paste0(locationOfProteinDataSets,"Protein_Counts_Gluconate.csv"),
                 header=T,  fill = TRUE, sep = ",",skip=0)

temp01<-gsub("Sample\\.","",colnames(data)[2:ncol(data)])
temp01<-gsub("\\.stationary|\\.EXP","",temp01)
temp01<-gsub("_R.*","",temp01)
temp01<-paste0(sprintf("%03d",as.numeric(gsub("\\..*","",temp01))),
               ".",
               gsub("*.*\\.","",temp01))

temp01<-gsub("\\.a","_0",temp01)
temp01<-gsub("\\.b","_1",temp01)

nameList=c("gene_id",temp01)
Gluconate_allzero_rows=which(rowSums((!is.na(data))+0)==1)


for(counter01 in 2:length(data)){
  print(nameList[counter01])
  data %>% dplyr::select(1,counter01) -> temp
  colnames(temp) <-c( "gene_id" , "count" )
  temp$count=as.numeric(temp$count)
  fileName=paste0("MURI_",nameList[counter01],"_raw_protein_count.txt")
  write.table(temp,file=paste0(locationOfProteinDataSets,fileName),
              sep="\t",row.names=FALSE,quote = F)
  remove(temp)
}
#***************************


#***************************
# Lactate data
data<-read.table(file=paste0(locationOfProteinDataSets,"Protein_Counts_Lactate.csv"),
                 header=T,  fill = TRUE, sep = ",",skip=0)
temp01<-gsub("Sample","",colnames(data)[2:ncol(data)])
temp01<-paste0(sprintf("%03d",as.numeric(gsub("_.*","",temp01))),
               "_",
               gsub("*.*_","",temp01))
nameList=c("gene_id",temp01)
Lactate_allzero_rows=which(rowSums((!is.na(data))+0)==1)


for(counter01 in 2:length(data)){
  print(nameList[counter01])
  data %>% dplyr::select(1,counter01) -> temp
  colnames(temp) <-c( "gene_id" , "count" )
  temp$count=as.numeric(temp$count)
  fileName=paste0("MURI_",nameList[counter01],"_raw_protein_count.txt")
  write.table(temp,file=paste0(locationOfProteinDataSets,fileName),
              sep="\t",row.names=FALSE,quote = F)
  remove(temp)
}
#***************************


#***************************
# Missing data
missingData1<-read.table(file=paste0(locationOfProteinDataSets,"MURI_52_1a.csv"),
                         header=T,  fill = TRUE, sep = "\t",skip=0)
missingData1 <- missingData1[-nrow(missingData1),c(1,4)]
colnames(missingData1) <-c( "gene_id" , "count" )
missingData1$count=as.numeric(missingData1$count)
fileName=paste0("MURI_","052_0","_raw_protein_count.txt")
write.table(missingData1,file=paste0(locationOfProteinDataSets,fileName),
            sep="\t",row.names=FALSE,quote = F)


missingData2<-read.table(file=paste0(locationOfProteinDataSets,"MURI_57_1a.csv"),
                         header=T,  fill = TRUE, sep = "\t",skip=0)
missingData2 <- missingData2[-nrow(missingData2),c(1,4)]
colnames(missingData2) <-c( "gene_id" , "count" )
missingData2$count=as.numeric(missingData2$count)
fileName=paste0("MURI_","057_0","_raw_protein_count.txt")
write.table(missingData2,file=paste0(locationOfProteinDataSets,fileName),
            sep="\t",row.names=FALSE,quote = F)


missingData3<-read.table(file=paste0(locationOfProteinDataSets,"MURI_106_1a.csv"),
                         header=T,  fill = TRUE, sep = "\t",skip=0)
missingData3 <- missingData3[-nrow(missingData3),c(1,4)]
colnames(missingData3) <-c( "gene_id" , "count" )
missingData3$count=as.numeric(missingData3$count)
fileName=paste0("MURI_","106_0","_raw_protein_count.txt")
write.table(missingData3,file=paste0(locationOfProteinDataSets,fileName),
            sep="\t",row.names=FALSE,quote = F)


missingData4<-read.table(file=paste0(locationOfProteinDataSets,"MURI_110_1b.csv"),
                         header=T,  fill = TRUE, sep = "\t",skip=0)
missingData4 <- missingData4[-nrow(missingData4),c(1,4)]
colnames(missingData4) <-c( "gene_id" , "count" )
missingData4$count=as.numeric(missingData4$count)
fileName=paste0("MURI_","110_0","_raw_protein_count.txt")
write.table(missingData4,file=paste0(locationOfProteinDataSets,fileName),
            sep="\t",row.names=FALSE,quote = F)
#***************************



#*********************Q******************
# load and combine AG3C_Protein_data to produce RNA matrix
Protein_FileList=dir(path=locationOfProteinDataSets)

RawProtein_FileList=grep("raw_protein_count", Protein_FileList, value = TRUE)
RawProtein_colNames=gsub("_raw_protein_count.txt","",RawProtein_FileList)


for(counter01 in 1: length(RawProtein_FileList)){
  print(counter01)
  
  temp<-read.table(file=paste0(locationOfProteinDataSets,RawProtein_FileList[counter01]),
                   header=TRUE)
  row.names(temp)<-NULL
  colnames(temp)<-c("gene_id",RawProtein_colNames[counter01])
  temp[[2]]<-as.numeric(temp[[2]])
  
  # replace NA with zero in the imported data
  temp2=temp[[2]]; temp2[is.na(temp2)]<-0;temp[[2]]=temp2;
  if(counter01==1){
    proteinMatrix=temp
  }
  else{
    proteinMatrix=full_join(proteinMatrix,temp)
  }
  print(c(nrow(temp),sum(!is.na(proteinMatrix[,counter01+1]))))
}
print(colSums(is.na(proteinMatrix)))


remove(temp)
proteinMatrix<-proteinMatrix[sort(colnames(proteinMatrix))]

# remove the MURI_082_0 from protein data
# in protein data MURI_82 is all NA so need to be fixed
proteinMatrix %>% dplyr::select(-MURI_082_0)->proteinMatrix

# put the proteins in to order and remove contaminants
proteinMatrix %>%
  dplyr::arrange(gene_id)%>%
  dplyr::filter(!grepl('CON', gene_id))->proteinMatrix


# generate a protein matrix without NA
proteinMatrix_w_NA=proteinMatrix
proteinMatrix_wo_NA <- proteinMatrix[rowSums(is.na(proteinMatrix)) == 0,]
proteinMatrix[is.na(proteinMatrix)]<-0

# Save as csv files
savedFilename=paste0("../a_results/","proteinMatrix.csv")
write.csv(proteinMatrix, file = savedFilename, row.names = FALSE)

savedFilename=paste0("../a_results/","proteinMatrix_w_NA.csv")
write.csv(proteinMatrix_w_NA, file = savedFilename, row.names = FALSE)

savedFilename=paste0("../a_results/","proteinMatrix_wo_NA.csv")
write.csv(proteinMatrix_wo_NA, file = savedFilename, row.names = FALSE)