# Test DeSeq on batch


###*****************************
# Do the DeSeq2 test
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
test_for="Mg_mM_Levels"
DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ",test_for))
differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
(res_Mg <- DESeq2::results(object = differentialGeneAnalResults, pAdjustMethod ="fdr"))
mcols(res_Mg, use.names=TRUE)  
res_df_Mg<-as.data.frame(res_Mg)

DESeq2::summary.DESeqResults(object = res_Mg,alpha = 0.05)
###*****************************




###*****************************
# Do the DeSeq2 test
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
test_for="Mg_mM_Levels"
DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ",test_for," + batchNumber"))
differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
(res_MgB <- DESeq2::results(object = differentialGeneAnalResults, pAdjustMethod ="fdr"))
mcols(res_MgB, use.names=TRUE)  
res_df_MgB<-as.data.frame(res_MgB)

DESeq2::summary.DESeqResults(object = res_MgB,alpha = 0.05)
###*****************************



###*****************************
# Do the DeSeq2 test
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
test_for="Mg_mM_Levels"
DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ",test_for," + batchNumber"))
differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
(res_MgBconhb <- DESeq2::results(object = differentialGeneAnalResults, pAdjustMethod ="fdr", contrast = c(test_for,"highMg","baseMg")))
mcols(res_MgBconhb, use.names=TRUE)  
res_df_MgBconhb<-as.data.frame(res_MgBconhb)

DESeq2::summary.DESeqResults(object = res_MgBconhb,alpha = 0.05)
###*****************************


###*****************************
# Do the DeSeq2 test
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
test_for="Mg_mM_Levels"
DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ",test_for," + batchNumber"))
differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
(res_MgBconlb <- DESeq2::results(object = differentialGeneAnalResults, pAdjustMethod ="fdr", contrast = c(test_for,"lowMg","baseMg")))
mcols(res_MgBconlb, use.names=TRUE)  
res_df_MgBconlb<-as.data.frame(res_MgBconlb)

DESeq2::summary.DESeqResults(object = res_MgBconlb,alpha = 0.05)
###*****************************







###*****************************
# Do the DeSeq2 test
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
test_for="Na_mM_Levels"
DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ",test_for))
differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
(res_Na <- DESeq2::results(object = differentialGeneAnalResults, pAdjustMethod ="fdr"))
mcols(res_Na, use.names=TRUE)  
res_df_Na<-as.data.frame(res_Na)

DESeq2::summary.DESeqResults(object = res_Na,alpha = 0.05)
###*****************************


###*****************************
# Do the DeSeq2 test
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
test_for="Na_mM_Levels"
DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ",test_for))
differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
(res_NaCon <- DESeq2::results(object = differentialGeneAnalResults, pAdjustMethod ="fdr", contrast = c(test_for,"highNa","baseNa")))
mcols(res_NaCon, use.names=TRUE)  
res_df_NaCon<-as.data.frame(res_NaCon)

DESeq2::summary.DESeqResults(object = res_NaCon,alpha = 0.05)
###*****************************




###*****************************
# Do the DeSeq2 test
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
test_for="Na_mM_Levels"
DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ",test_for," + batchNumber"))
differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
(res_NaB <- DESeq2::results(object = differentialGeneAnalResults, pAdjustMethod ="fdr"))
mcols(res_NaB, use.names=TRUE)  
res_df_NaB<-as.data.frame(res_NaB)

DESeq2::summary.DESeqResults(object = res_NaB,alpha = 0.05)
###*****************************


###*****************************
# Do the DeSeq2 test
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
test_for="Na_mM_Levels"
DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ",test_for," + batchNumber"))
differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
(res_NaBCon <- DESeq2::results(object = differentialGeneAnalResults, pAdjustMethod ="fdr", contrast = c(test_for,"highNa","baseNa")))
mcols(res_NaBCon, use.names=TRUE)  
res_df_NaBCon<-as.data.frame(res_NaBCon)

DESeq2::summary.DESeqResults(object = res_NaBCon,alpha = 0.05)
###*****************************


###*****************************
# Do the DeSeq2 test
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
test_for="Na_mM_Levels"
DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ","batchNumber + ",test_for))
differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
(res_NaBrevCon <- DESeq2::results(object = differentialGeneAnalResults, pAdjustMethod ="fdr", contrast = c(test_for,"highNa","baseNa")))
mcols(res_NaBrevCon, use.names=TRUE)  
res_df_NaBrevCon<-as.data.frame(res_NaBrevCon)

DESeq2::summary.DESeqResults(object = res_NaBrevCon,alpha = 0.05)
###*****************************



###*****************************
# Do the DeSeq2 test
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
test_for="Na_mM_Levels"
DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ",test_for," + batchNumber"))
differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
(res_NaBConrev <- DESeq2::results(object = differentialGeneAnalResults, pAdjustMethod ="fdr", contrast = c(test_for,"baseNa","highNa")))
mcols(res_NaBConrev, use.names=TRUE)  
res_df_NaBConrev<-as.data.frame(res_NaBConrev)

DESeq2::summary.DESeqResults(object = res_NaBConrev,alpha = 0.05)
###*****************************








###*****************************
# Do the DeSeq2 test
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
test_for="carbonSource"
DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ",test_for))
differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
(res_CarbonCon <- DESeq2::results(object = differentialGeneAnalResults, pAdjustMethod ="fdr", contrast = c(test_for,"glycerol","glucose")))
mcols(res_CarbonCon, use.names=TRUE)  
res_df_CarbonCon<-as.data.frame(res_CarbonCon)

DESeq2::summary.DESeqResults(object = res_CarbonCon,alpha = 0.05)
###*****************************


###*****************************
# Do the DeSeq2 test
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
test_for="carbonSource"
DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ",test_for," + batchNumber"))
differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
(res_CarbonBCon <- DESeq2::results(object = differentialGeneAnalResults, pAdjustMethod ="fdr", contrast = c(test_for,"glycerol","glucose")))
mcols(res_CarbonBCon, use.names=TRUE)  
res_df_CarbonBCon<-as.data.frame(res_CarbonBCon)

DESeq2::summary.DESeqResults(object = res_CarbonBCon,alpha = 0.05)
###*****************************













###*****************************
# Test for totally seperate fake batch (which should not effect anything)
# Do the DeSeq2 test
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
test_for="carbonSource"
DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ",test_for))
differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
(res_CarbonCon <- DESeq2::results(object = differentialGeneAnalResults, pAdjustMethod ="fdr", contrast = c(test_for,"glycerol","glucose")))
mcols(res_CarbonCon, use.names=TRUE)  
res_df_CarbonCon<-as.data.frame(res_CarbonCon)

DESeq2::summary.DESeqResults(object = res_CarbonCon,alpha = 0.05)
###*****************************

###*****************************
# Do the DeSeq2 test
# c("Mg_mM_Levels", "Na_mM_Levels", "growthPhase", "carbonSource")
test_for="carbonSource"
DESeq2::design(deseq_DataObj)<- as.formula(paste0("~ ",test_for," + fake_batch"))
differentialGeneAnalResults<-DESeq2::DESeq(deseq_DataObj)
(res_CarbonfBCon <- DESeq2::results(object = differentialGeneAnalResults, pAdjustMethod ="fdr", contrast = c(test_for,"glycerol","glucose")))
mcols(res_CarbonfBCon, use.names=TRUE)  
res_df_CarbonfBCon<-as.data.frame(res_CarbonfBCon)

DESeq2::summary.DESeqResults(object = res_CarbonfBCon,alpha = 0.05)
###*****************************
