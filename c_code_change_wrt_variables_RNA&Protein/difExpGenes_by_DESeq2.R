# Generate diffreential analysis of count data by using DESeq2 package

###*****************************
# INITIAL COMMANDS TO RESET THE SYSTEM

rm(list = ls())
if (is.integer(dev.list())){dev.off()}
cat("\014")

# Set Working Directory
setwd('/Users/umut/GitHub/ecoli_multiple_growth_conditions/c_code_change_wrt_variables_RNA&Protein') # mac computer
###*****************************

###*****************************
# DOWNLOAD LIBRARIES
library("Biobase") 
library("DESeq2")
library("dplyr")
###*****************************


###*****************************
# LOAD DATA
unnormalized_mrna_data=read.csv(file=paste0("../Processed_RNA/","rnaMatrix_mRNA.csv")) #The main unnormalized mRNA matrix
unnormalized_protein_data_wo_NA=read.csv(file=paste0("../Processed_Protein/","proteinMatrix_wo_NA.csv")) # The unnormalized protein matrix contains proteins that are measured for all datasets
meta_rna=read.csv(paste0("../a_results/","metaRNA.csv")) #Main Meta data Frame
meta_protein=read.csv(paste0("../a_results/","metaProtein.csv")) #Main Meta data Frame
###*****************************


###*****************************
# PREPERATION FOR mRNA
rownames(unnormalized_mrna_data)<-unnormalized_mrna_data[["gene_ID"]]
unnormalized_mrna_data %>% 
  dplyr::select(-c(gene_Type, gene_ID))->unnormalized_mrna_data

# Get rid of bad data sets
badDataFilterOutDataSet=c("MURI_029","MURI_067","MURI_075","MURI_084",
                          "MURI_136","MURI_138")
badDataRemainingDataSets=base::setdiff(as.vector(meta_rna$dataSet),badDataFilterOutDataSet)

unnormalized_mrna_data=unnormalized_mrna_data[,badDataRemainingDataSets]
mainDataFrameRNA=as.matrix(unnormalized_mrna_data)
remove(unnormalized_mrna_data)

meta_rna %>%
  dplyr::filter(dataSet %in% badDataRemainingDataSets)->meta_rna
conditionRNA=meta_rna

remove(meta_rna)
###*****************************


###*****************************
# PREPERATION FOR Protein
rownames(unnormalized_protein_data_wo_NA)<-unnormalized_protein_data_wo_NA[["gene_id"]]
unnormalized_protein_data_wo_NA %>% 
  dplyr::select(-c(gene_id))->unnormalized_protein_data_wo_NA
mainDataFrameProtein=unnormalized_protein_data_wo_NA

# Calculate the average of repeated conditions
mainDataFrameProtein<-as.data.frame(mainDataFrameProtein)
mainDataFrameProtein %>% mutate(gene_id=row.names(mainDataFrameProtein))->mainDataFrameProteinTidy
endVal=ncol(mainDataFrameProtein)

mainDataFrameProteinTidy %>%
  tidyr::gather(conditionProtein,read,1:endVal)%>%
  dplyr::mutate(conditionShort=gsub("_[[:digit:]]$","",conditionProtein))->mainDataFrameProteinTidy

mainDataFrameProteinTidy %>% 
  group_by(gene_id,conditionShort)%>%
  dplyr::summarise(read=mean(read))->mainDataFrameProteinSummary

mainDataFrameProteinSummary %>%
  dplyr::group_by()%>%  
  tidyr::spread(key = conditionShort, value = read)->mainDataFrameProteinSummary

row.names(mainDataFrameProteinSummary)<-mainDataFrameProteinSummary$gene_id
mainDataFrameProteinSummary %>%
  dplyr::select(-gene_id)->mainDataFrameProteinSummary
mainDataFrameProtein=as.matrix(round(mainDataFrameProteinSummary))

remove(mainDataFrameProteinSummary,
       mainDataFrameProteinTidy,
       unnormalized_protein_data_wo_NA)
###*****************************


###*****************************
#Reduce number of rows in condition and get rid of repeats
conditionProtein=meta_protein
conditionProtein=conditionProtein[which(substr(as.vector(conditionProtein$dataSet),10,10) %in% c(0)),]
conditionProtein$dataSet<-gsub("_[[:digit:]]$","",as.vector(conditionProtein$dataSet))
conditionProtein$sampleNum<-gsub("_[[:digit:]]$","",as.vector(conditionProtein$sampleNum))
remove(meta_protein)
###*****************************


###*****************************
# P1 Function
generate_p1_object=function(unnormalized_sampleData, meta_condition){
  
  # Generate "p1 Data Object"
  unnormalized_sampleData_p1=unnormalized_sampleData+1;
  obj_p1=DESeqDataSetFromMatrix(countData = unnormalized_sampleData_p1,
                                      colData = meta_condition,design = ~ condition)
  levelList<-as.character(unique(meta_condition[["condition"]]));
  colData(obj_p1)$condition=factor(colData(obj_p1)$condition, levels=levelList)
  
  # Calculate Size Factors
  obj_p1=estimateSizeFactors(obj_p1)
  sizeFactors_p1=sizeFactors(obj_p1)
  
  # Generate "regular Data Object"
  obj=DESeqDataSetFromMatrix(countData = unnormalized_sampleData, colData = meta_condition,design = ~ condition)
  levelList<-as.character(unique(meta_condition[["condition"]]));
  colData(obj)$condition=factor(colData(obj)$condition, levels=levelList)
  
  # Import Size Factors From "p1 Data Object"
  sizeFactors(obj) <- sizeFactors_p1
  
  return(obj)
}
###*****************************


###*****************************
# test for High Mg Exp vs Glucose Exp RNA
conditionRNA %>%
  filter((experiment=="glucose_time_course" & growthPhase=="exponential") |
           (Mg_mM_Levels=="highMg" & growthPhase=="exponential"))%>%
  mutate(condition=Mg_mM_Levels)->conditionSpesificRNA
dataSetList=as.vector(conditionSpesificRNA$dataSet)
rownames(conditionSpesificRNA)<-dataSetList

mainDataFrameSpesificRNA<-mainDataFrameRNA[,dataSetList]
obj<-generate_p1_object(mainDataFrameSpesificRNA,conditionSpesificRNA)
obj$condition <- relevel(obj$condition, ref="baseMg")
obj<-DESeq2::DESeq(obj)
resHighMgRNAExp<-as.data.frame(DESeq2::results(obj))
head(resHighMgRNAExp)
###*****************************


###*****************************
# test for Low Mg Exp vs Glucose Exp RNA
conditionRNA %>%
  filter((experiment=="glucose_time_course" & growthPhase=="exponential") |
           (Mg_mM_Levels=="lowMg" & growthPhase=="exponential"))%>%
  mutate(condition=Mg_mM_Levels)->conditionSpesificRNA
dataSetList=as.vector(conditionSpesificRNA$dataSet)
rownames(conditionSpesificRNA)<-dataSetList

mainDataFrameSpesificRNA<-mainDataFrameRNA[,dataSetList]
obj<-generate_p1_object(mainDataFrameSpesificRNA,conditionSpesificRNA)
obj$condition <- relevel(obj$condition, ref="baseMg")
obj<-DESeq2::DESeq(obj)
resLowMgRNAExp<-as.data.frame(DESeq2::results(obj))
head(resLowMgRNAExp)
###*****************************


###*****************************
# test for High Na Exp vs Glucose Exp RNA
conditionRNA %>%
  filter((experiment=="glucose_time_course" & growthPhase=="exponential") |
           (Na_mM_Levels=="highNa" & growthPhase=="exponential"))%>%
  mutate(condition=Na_mM_Levels)->conditionSpesificRNA
dataSetList=as.vector(conditionSpesificRNA$dataSet)
rownames(conditionSpesificRNA)<-dataSetList

mainDataFrameSpesificRNA<-mainDataFrameRNA[,dataSetList]
obj<-generate_p1_object(mainDataFrameSpesificRNA,conditionSpesificRNA)
obj$condition <- relevel(obj$condition, ref="baseNa")
obj<-DESeq2::DESeq(obj)
resHighNaRNAExp<-as.data.frame(DESeq2::results(obj))
head(resHighNaRNAExp)
###*****************************


###*****************************
# test for glycerol Exp vs Glucose Exp RNA
conditionRNA %>%
  filter((experiment=="glucose_time_course" & growthPhase=="exponential") |
           (carbonSource=="glycerol" & growthPhase=="exponential"))%>%
  mutate(condition=carbonSource)->conditionSpesificRNA
dataSetList=as.vector(conditionSpesificRNA$dataSet)
rownames(conditionSpesificRNA)<-dataSetList

mainDataFrameSpesificRNA<-mainDataFrameRNA[,dataSetList]
obj<-generate_p1_object(mainDataFrameSpesificRNA,conditionSpesificRNA)
obj$condition <- relevel(obj$condition, ref="glucose")
obj<-DESeq2::DESeq(obj)
resGlyRNAExp<-as.data.frame(DESeq2::results(obj))
head(resGlyRNAExp)
###*****************************


###*****************************
# test for gluconate Exp vs Glucose Exp RNA
conditionRNA %>%
  filter((experiment=="glucose_time_course" & growthPhase=="exponential") |
           (carbonSource=="gluconate"& growthPhase=="exponential"))%>%
  mutate(condition=carbonSource)->conditionSpesificRNA
dataSetList=as.vector(conditionSpesificRNA$dataSet)
rownames(conditionSpesificRNA)<-dataSetList

mainDataFrameSpesificRNA<-mainDataFrameRNA[,dataSetList]
obj<-generate_p1_object(mainDataFrameSpesificRNA,conditionSpesificRNA)
obj$condition <- relevel(obj$condition, ref="glucose")
obj<-DESeq2::DESeq(obj)
resGluRNAExp<-as.data.frame(DESeq2::results(obj))
head(resGluRNAExp)
###*****************************


###*****************************
# test for lactate Exp vs Glucose Exp RNA
conditionRNA %>%
  filter((experiment=="glucose_time_course" & growthPhase=="exponential") |
           (carbonSource=="lactate" & growthPhase=="exponential"))%>%
  mutate(condition=carbonSource)->conditionSpesificRNA
dataSetList=as.vector(conditionSpesificRNA$dataSet)
rownames(conditionSpesificRNA)<-dataSetList

mainDataFrameSpesificRNA<-mainDataFrameRNA[,dataSetList]
obj<-generate_p1_object(mainDataFrameSpesificRNA,conditionSpesificRNA)
obj$condition <- relevel(obj$condition, ref="glucose")
obj<-DESeq2::DESeq(obj)
resLacRNAExp<-as.data.frame(DESeq2::results(obj))
head(resLacRNAExp)
###*****************************
###*****************************


###*****************************
###*****************************
# test for High Mg Exp vs Glucose Exp Protein
conditionProtein %>%
  filter((experiment=="glucose_time_course" & growthPhase=="exponential") |
           (Mg_mM_Levels=="highMg" & growthPhase=="exponential"))%>%
  mutate(condition=Mg_mM_Levels)->conditionSpesificProtein
dataSetList=as.vector(conditionSpesificProtein$dataSet)
rownames(conditionSpesificProtein)<-dataSetList

mainDataFrameSpesificProtein<-mainDataFrameProtein[,dataSetList]
obj<-generate_p1_object(mainDataFrameSpesificProtein,conditionSpesificProtein)
obj$condition <- relevel(obj$condition, ref="baseMg")
obj<-DESeq2::DESeq(obj)
resHighMgProteinExp<-as.data.frame(DESeq2::results(obj))
head(resHighMgProteinExp)
###*****************************


###*****************************
# test for Low Mg Exp vs Glucose Exp Protein
conditionProtein %>%
  filter((experiment=="glucose_time_course" & growthPhase=="exponential") |
           (Mg_mM_Levels=="lowMg" & growthPhase=="exponential"))%>%
  mutate(condition=Mg_mM_Levels)->conditionSpesificProtein
dataSetList=as.vector(conditionSpesificProtein$dataSet)
rownames(conditionSpesificProtein)<-dataSetList

mainDataFrameSpesificProtein<-mainDataFrameProtein[,dataSetList]
obj<-generate_p1_object(mainDataFrameSpesificProtein,conditionSpesificProtein)
obj$condition <- relevel(obj$condition, ref="baseMg")
obj<-DESeq2::DESeq(obj)
resLowMgProteinExp<-as.data.frame(DESeq2::results(obj))
head(resLowMgProteinExp)
###*****************************

###*****************************
# test for High Na Exp vs Glucose Exp Protein
conditionProtein %>%
  filter((experiment=="glucose_time_course" & growthPhase=="exponential") |
           (Na_mM_Levels=="highNa" & growthPhase=="exponential"))%>%
  mutate(condition=Na_mM_Levels)->conditionSpesificProtein
dataSetList=as.vector(conditionSpesificProtein$dataSet)
rownames(conditionSpesificProtein)<-dataSetList

mainDataFrameSpesificProtein<-mainDataFrameProtein[,dataSetList]
obj<-generate_p1_object(mainDataFrameSpesificProtein,conditionSpesificProtein)
obj$condition <- relevel(obj$condition, ref="baseNa")
obj<-DESeq2::DESeq(obj)
resHighNaProteinExp<-as.data.frame(DESeq2::results(obj))
head(resHighNaProteinExp)
###*****************************


###*****************************
# test for Glycerol Exp vs Glucose Exp Protein
conditionProtein %>%
  filter((experiment=="glucose_time_course" & growthPhase=="exponential") |
           (carbonSource=="glycerol" & growthPhase=="exponential"))%>%
  mutate(condition=carbonSource)->conditionSpesificProtein
dataSetList=as.vector(conditionSpesificProtein$dataSet)
rownames(conditionSpesificProtein)<-dataSetList

mainDataFrameSpesificProtein<-mainDataFrameProtein[,dataSetList]
obj<-generate_p1_object(mainDataFrameSpesificProtein,conditionSpesificProtein)
obj$condition <- relevel(obj$condition, ref="glucose")
obj<-DESeq2::DESeq(obj)
resGlyProteinExp<-as.data.frame(DESeq2::results(obj))
head(resGlyProteinExp)
###*****************************


##*****************************
###*****************************


###*****************************
###*****************************
###*****************************
# test for High Mg Sta vs Glucose Sta RNA
conditionRNA %>%
  filter((experiment=="glucose_time_course" & growthPhase=="stationary") |
           (Mg_mM_Levels=="highMg" & growthPhase=="stationary"))%>%
  mutate(condition=Mg_mM_Levels)->conditionSpesificRNA
dataSetList=as.vector(conditionSpesificRNA$dataSet)
rownames(conditionSpesificRNA)<-dataSetList

mainDataFrameSpesificRNA<-mainDataFrameRNA[,dataSetList]
obj<-generate_p1_object(mainDataFrameSpesificRNA,conditionSpesificRNA)
obj$condition <- relevel(obj$condition, ref="baseMg")
obj<-DESeq2::DESeq(obj)
resHighMgRNASta<-as.data.frame(DESeq2::results(obj))
head(resHighMgRNASta)
###*****************************


###*****************************
# test for Low Mg Sta vs Glucose Sta RNA
conditionRNA %>%
  filter((experiment=="glucose_time_course" & growthPhase=="stationary") |
           (Mg_mM_Levels=="lowMg" & growthPhase=="stationary"))%>%
  mutate(condition=Mg_mM_Levels)->conditionSpesificRNA
dataSetList=as.vector(conditionSpesificRNA$dataSet)
rownames(conditionSpesificRNA)<-dataSetList

mainDataFrameSpesificRNA<-mainDataFrameRNA[,dataSetList]
obj<-generate_p1_object(mainDataFrameSpesificRNA,conditionSpesificRNA)
obj$condition <- relevel(obj$condition, ref="baseMg")
obj<-DESeq2::DESeq(obj)
resLowMgRNASta<-as.data.frame(DESeq2::results(obj))
head(resLowMgRNASta)
###*****************************


###*****************************
# test for High Na Sta vs Glucose Sta RNA
conditionRNA %>%
  filter((experiment=="glucose_time_course" & growthPhase=="stationary") |
           (Na_mM_Levels=="highNa" & growthPhase=="stationary"))%>%
  mutate(condition=Na_mM_Levels)->conditionSpesificRNA
dataSetList=as.vector(conditionSpesificRNA$dataSet)
rownames(conditionSpesificRNA)<-dataSetList

mainDataFrameSpesificRNA<-mainDataFrameRNA[,dataSetList]
obj<-generate_p1_object(mainDataFrameSpesificRNA,conditionSpesificRNA)
obj$condition <- relevel(obj$condition, ref="baseNa")
obj<-DESeq2::DESeq(obj)
resHighNaRNASta<-as.data.frame(DESeq2::results(obj))
head(resHighNaRNASta)
###*****************************


###*****************************
# test for glycerol Sta vs Glucose Sta
conditionRNA %>%
  filter((experiment=="glucose_time_course" & growthPhase=="stationary") |
           (carbonSource=="glycerol" & growthPhase=="stationary"))%>%
  mutate(condition=carbonSource)->conditionSpesificRNA
dataSetList=as.vector(conditionSpesificRNA$dataSet)
rownames(conditionSpesificRNA)<-dataSetList

mainDataFrameSpesificRNA<-mainDataFrameRNA[,dataSetList]
obj<-generate_p1_object(mainDataFrameSpesificRNA,conditionSpesificRNA)
obj$condition <- relevel(obj$condition, ref="glucose")
obj<-DESeq2::DESeq(obj)
resGlyRNASta<-as.data.frame(DESeq2::results(obj))
head(resGlyRNASta)
###*****************************


###*****************************
# test for gluconate Sta vs Glucose Sta RNA
conditionRNA %>%
  filter((experiment=="glucose_time_course" & growthPhase=="stationary") |
           (carbonSource=="gluconate"& growthPhase=="stationary"))%>%
  mutate(condition=carbonSource)->conditionSpesificRNA
dataSetList=as.vector(conditionSpesificRNA$dataSet)
rownames(conditionSpesificRNA)<-dataSetList

mainDataFrameSpesificRNA<-mainDataFrameRNA[,dataSetList]
obj<-generate_p1_object(mainDataFrameSpesificRNA,conditionSpesificRNA)
obj$condition <- relevel(obj$condition, ref="glucose")
obj<-DESeq2::DESeq(obj)
resGluRNASta<-as.data.frame(DESeq2::results(obj))
head(resGluRNASta)
###*****************************


###*****************************
# test for lactate Sta vs Glucose Sta RNA
conditionRNA %>%
  filter((experiment=="glucose_time_course" & growthPhase=="stationary") |
           (carbonSource=="lactate" & growthPhase=="stationary"))%>%
  mutate(condition=carbonSource)->conditionSpesificRNA
dataSetList=as.vector(conditionSpesificRNA$dataSet)
rownames(conditionSpesificRNA)<-dataSetList

mainDataFrameSpesificRNA<-mainDataFrameRNA[,dataSetList]
obj<-generate_p1_object(mainDataFrameSpesificRNA,conditionSpesificRNA)
obj$condition <- relevel(obj$condition, ref="glucose")
obj<-DESeq2::DESeq(obj)
resLacRNASta<-as.data.frame(DESeq2::results(obj))
head(resLacRNASta)
###*****************************
###*****************************


###*****************************
# test for High Mg Sta vs Glucose Sta Protein
conditionProtein %>%
  filter((experiment=="glucose_time_course" & growthPhase=="stationary") |
           (Mg_mM_Levels=="highMg" & growthPhase=="stationary"))%>%
  mutate(condition=Mg_mM_Levels)->conditionSpesificProtein
dataSetList=as.vector(conditionSpesificProtein$dataSet)
rownames(conditionSpesificProtein)<-dataSetList

mainDataFrameSpesificProtein<-mainDataFrameProtein[,dataSetList]
obj<-generate_p1_object(mainDataFrameSpesificProtein,conditionSpesificProtein)
obj$condition <- relevel(obj$condition, ref="baseMg")
obj<-DESeq2::DESeq(obj)
resHighMgProteinSta<-as.data.frame(DESeq2::results(obj))
head(resHighMgProteinSta)
###*****************************


###*****************************
# test for Low Mg Sta vs Glucose Sta Protein
conditionProtein %>%
  filter((experiment=="glucose_time_course" & growthPhase=="stationary") |
           (Mg_mM_Levels=="lowMg" & growthPhase=="stationary"))%>%
  mutate(condition=Mg_mM_Levels)->conditionSpesificProtein
dataSetList=as.vector(conditionSpesificProtein$dataSet)
rownames(conditionSpesificProtein)<-dataSetList

mainDataFrameSpesificProtein<-mainDataFrameProtein[,dataSetList]
obj<-generate_p1_object(mainDataFrameSpesificProtein,conditionSpesificProtein)
obj$condition <- relevel(obj$condition, ref="baseMg")
obj<-DESeq2::DESeq(obj)
resLowMgProteinSta<-as.data.frame(DESeq2::results(obj))
head(resLowMgProteinSta)
###*****************************

###*****************************
# test for High Na Sta vs Glucose Sta Protein
conditionProtein %>%
  filter((experiment=="glucose_time_course" & growthPhase=="stationary") |
           (Na_mM_Levels=="highNa" & growthPhase=="stationary"))%>%
  mutate(condition=Na_mM_Levels)->conditionSpesificProtein
dataSetList=as.vector(conditionSpesificProtein$dataSet)
rownames(conditionSpesificProtein)<-dataSetList

mainDataFrameSpesificProtein<-mainDataFrameProtein[,dataSetList]
obj<-generate_p1_object(mainDataFrameSpesificProtein,conditionSpesificProtein)
obj$condition <- relevel(obj$condition, ref="baseNa")
obj<-DESeq2::DESeq(obj)
resHighNaProteinSta<-as.data.frame(DESeq2::results(obj))
head(resHighNaProteinSta)
###*****************************


###*****************************
# test for Glycerol Sta vs Glucose Sta Protein
conditionProtein %>%
  filter((experiment=="glucose_time_course" & growthPhase=="stationary") |
           (carbonSource=="glycerol" & growthPhase=="stationary"))%>%
  mutate(condition=carbonSource)->conditionSpesificProtein
dataSetList=as.vector(conditionSpesificProtein$dataSet)
rownames(conditionSpesificProtein)<-dataSetList

mainDataFrameSpesificProtein<-mainDataFrameProtein[,dataSetList]
obj<-generate_p1_object(mainDataFrameSpesificProtein,conditionSpesificProtein)
obj$condition <- relevel(obj$condition, ref="glucose")
obj<-DESeq2::DESeq(obj)
resGlyProteinSta<-as.data.frame(DESeq2::results(obj))
head(resGlyProteinSta)
###*****************************

###*****************************
###*****************************
save(resHighMgRNAExp,resLowMgRNAExp,resHighNaRNAExp,resGluRNAExp,resGlyRNAExp,resLacRNAExp,
     resHighMgProteinExp,resLowMgProteinExp,resHighNaProteinExp,resGlyProteinExp,
     resHighMgRNASta,resLowMgRNASta,resHighNaRNASta,resGluRNASta,resGlyRNASta,resLacRNASta,
     resHighMgProteinSta,resLowMgProteinSta,resHighNaProteinSta,resGlyProteinSta,
     file="../c_results/significantChanges.Rdata")
#
