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
              "d_code_Mf_Pathway_Analyze/"))} # mac computer
###*****************************


###*****************************
# DOWNLOAD LIBRARIES
require("dplyr")
require("tidyr")
require("ggplot2")
require("RColorBrewer")
require("scales")
require("cowplot")
require("ggrepel")
###*****************************


###*****************************
#Load Functions
source("../b_code_histogram_RNA&Protein/replace_fun.R")
###*****************************


###*****************************
# Load files
f1_data<-read.csv("../d_results/ez_P0.05Fold2_mrna_trT_set00_StcAllEx_SYAN_baseMgAllMg_baseNaAllNa_Exp_noMatchFilter_p1Sf_noNorm__batchNumberPLUSMg_mM_Levels__highMgVSbaseMg_kegg.csv")
f2_data<-read.csv("../d_results/ez_P0.05Fold2_protein_trT_set00_StcYtcNasAgrNgrMgh_SYAN_baseMgAllMg_baseNaAllNa_Exp_noMatchFilter_p1Sf_noNorm__batchNumberPLUScarbonSource__lactateVSglucose_kegg.csv")
###*****************************


###*****************************
minimumFold=min(f1_data$log2)
if(minimumFold>-1){minimumFold=-1}
maximumFold=max(f1_data$log2)
if(maximumFold<1){maximumFold=1}
###*****************************


###*****************************
fig01=ggplot(f1_data, aes( x=log2,y=KEGG_Path_Short)) +
  geom_point(colour="blue", size=2.5)+
  geom_vline(xintercept = c(log2(1/2),log2(2)), colour="orange", linetype = "longdash")+
  geom_vline(xintercept = c(log2(1)), colour="black", linetype = "longdash")+
  geom_text_repel(aes(label=gene_name),size=3, colour="Black", fontface="plain")+
  theme_bw()+
  scale_x_continuous(breaks=seq(floor(minimumFold),ceiling(maximumFold)))+
  xlab("Log2 Fold Change")+
  theme(axis.line.y = element_blank(),
        legend.position="bottom",
        axis.title.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=16),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14))

print(fig01)
###*****************************


###*****************************
minimumFold=min(f2_data$log2)
if(minimumFold>-1){minimumFold=-1}
maximumFold=max(f2_data$log2)
if(maximumFold<1){maximumFold=1}
###*****************************


###*****************************
fig02=ggplot(f2_data, aes( x=log2,y=KEGG_Path_Short)) +
  geom_point(colour="blue", size=2.5)+
  geom_vline(xintercept = c(log2(1/2),log2(2)), colour="orange", linetype = "longdash")+
  geom_vline(xintercept = c(log2(1)), colour="black", linetype = "longdash")+
  geom_text_repel(aes(label=gene_name),size=3, colour="Black", fontface="plain")+
  theme_bw()+
  scale_x_continuous(breaks=seq(floor(minimumFold),ceiling(maximumFold)))+
  xlab("Log2 Fold Change")+
  theme(axis.line.y = element_blank(),
        legend.position="bottom",
        axis.title.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=16),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14))

print(fig02)
###*****************************


###*****************************
# Combine figures
fig_comb<-cowplot::plot_grid(fig01,fig02,nrow = 2,ncol = 1,align = "v",labels = c("A","B"),scale = .9)

rowWidth=3
cowplot::save_plot(filename = paste0("../d_figures/combined.pdf"),
                   plot = fig_comb,
                   base_height = rowWidth,
                   ncol=2.6,
                   nrow=2.4,
                   limitsize = FALSE)

cowplot::save_plot(filename = paste0("../text/figures/fig8_combined.pdf"),
                   plot = fig_comb,
                   base_height = rowWidth,
                   ncol=2.6,
                   nrow=2.4,
                   limitsize = FALSE)
###*****************************