# Aim of the code is to generate figure for Metabolites

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
              "d_Mf_Pathway_Analyze/"))} # mac computer
###*****************************


###*****************************
# DOWNLOAD LIBRARIES
require("dplyr")
require("tidyr")
require("ggplot2")
require("RColorBrewer")
require("scales")
require("cowplot")
###*****************************


###*****************************
#Load Functions
source("../b_code_histogram_RNA&Protein/replace_fun.R")
###*****************************


###*****************************
# Load files
fileList=dir("../c_results/",pattern = "resDf*.*Mg_mM*")
for(counter01 in 1 : length(fileList))
{
  if (counter01 ==1)
  {
    mainLoadedFile=read.csv(file = paste0("../c_results/",fileList[counter01]))
  }
  if (counter01 !=1)
  {
    temp=read.csv(file = paste0("../c_results/",fileList[counter01]),header = TRUE)
    if(any(colnames(temp)=="mRNA_ID"))
    {temp %>% dplyr::select(-mRNA_ID,-b.)->temp}
    mainLoadedFile=rbind(mainLoadedFile,temp)
  }
}
###***************************** 


###*****************************
# Filter Data
mainLoadedFile %>% dplyr::filter(gene_name %in% c("phoQ","phoP"))->mainLoadedFile
###*****************************


###*****************************
# Rename entries
mainLoadedFile$vs=as.vector(mainLoadedFile$vs)
mainLoadedFile %>%
  dplyr::mutate(pick_data=ifelse(pick_data=="mrna","mRNA","Protein"))%>%
  dplyr::mutate(vs=ifelse(vs=="baseMghighMg","High Mg",vs),
                vs=ifelse(vs=="baseMglowMg","Low Mg",vs)) %>%
  dplyr::group_by(pick_data,growthPhase,vs)%>%
  dplyr::mutate(rank=seq(1:n()))->mainLoadedFile
###*****************************


###*****************************
# Generate Tables
fig01<-ggplot(mainLoadedFile, aes( y=rank,x="values"))+
  ylim(2.5,0.5)+
  scale_x_continuous(limits=c(0,1), expand=c(0,0)) +
  geom_text(aes(label=paste(gene_name, 
                            sprintf("P'=%.3e",padj),
                            sprintf("Log2FoldChange=%.3f",log2FoldChange),
                            sep="    "), x=0.05),size=3, hjust=0)+
  facet_grid(growthPhase+vs~pick_data) +
  panel_border() +
  theme(  axis.line = element_blank(), 
          axis.text.x = element_blank(), 
          axis.text.y = element_blank(),
          axis.ticks = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank())

print(fig01)
###*****************************


###*****************************
# Save Figures
cowplot::save_plot(filename = paste0("../d_figures/PhoQ_PhoP.pdf"),
                   plot = fig01,ncol = 1.6,nrow=.9,limitsize = FALSE)