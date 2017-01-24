# the difference between controling on growth and not controlling on growth

###*****************************
# Set Up java
dyn.load('/Library/Java/JavaVirtualMachines/jdk1.8.0_25.jdk/Contents/Home/jre/lib/server/libjvm.dylib')
###*****************************

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
{setwd('/Users/umut/GitHub/ecoli_multiple_growth_conditions/c_code_change_wrt_variables_RNA&Protein/')} # mac computer
###*****************************


###*****************************
# DOWNLOAD LIBRARIES
require("dplyr")
require("tidyr")

require("ggplot2")
require("cowplot")
require("VennDiagram")
require("gtable")
require("gridExtra")

require("RDAVIDWebService")
require(org.Hs.eg.db)
###*****************************


###*****************************
# Required functions
source("../b_code_histogram_RNA&Protein/replace_fun.R")
###*****************************


###*****************************
# Read significantly altered gene list
data<-read.csv(file = "../c_results/combinedDifferentiallyExpressedGenes_DeSeq.csv")
data %>% 
  dplyr::mutate(test_for=ifelse( test = grepl(pattern = "gluconate", x = testVSbase) , yes =  "Glc", no = NA ),
                test_for=ifelse( test = grepl(pattern = "glycerol", x = testVSbase) , yes = "Gly", no = test_for), 
                test_for=ifelse( test = grepl(pattern = "lactate", x = testVSbase) , yes = "Lac", no = test_for),
                test_for=ifelse( test = grepl(pattern = "highMg", x = testVSbase) , yes = "High Mg", no = test_for),
                test_for=ifelse( test = grepl(pattern = "lowMg", x = testVSbase) , yes = "Low Mg", no = test_for),
                test_for=ifelse( test = grepl(pattern = "highNa", x = testVSbase) , yes = "High Na", no = test_for)) %>%
  dplyr::mutate(controlledForGrowth = ifelse (test = grepl(pattern = "doublingTimeMinutes", x = investigatedEffect), yes = 1 , no = 0)) -> data 
###*****************************

###*****************************
# Divide data into 2 pieces
data %>%
  dplyr::filter(controlledForGrowth==0) %>%
  dplyr::select(genes, dataType, growthPhase, test_for) %>%
  dplyr::mutate(uncontrolled = 1)->uncontrolledGrowth

data %>%
  dplyr::filter(controlledForGrowth==1) %>%
  dplyr::select(genes, dataType, growthPhase, test_for) %>%
  dplyr::mutate(controlled = 1)->controlledGrowth
###*****************************


###*****************************
# join them
dplyr::full_join(uncontrolledGrowth, controlledGrowth) -> joined_data

joined_data$uncontrolled[is.na(joined_data$uncontrolled)]<-0
joined_data$controlled[is.na(joined_data$controlled)]<-0

joined_data %>%
  dplyr::mutate(exist_in = ifelse(test = (controlled==1 & uncontrolled==1) ,yes = "both", NA),
                exist_in = ifelse(test = (controlled==1 & uncontrolled==0) ,yes = "controlled for Growth", exist_in),
                exist_in = ifelse(test = (controlled==0 & uncontrolled==1) ,yes = "not controlled for Growth", exist_in))->joined_data 

joined_data$test_for <- factor(joined_data$test_for, levels=c("Low Mg", "High Mg", "High Na", "Gly", "Glc", "Lac"))
joined_data$exist_in <- factor(joined_data$exist_in, levels=c("controlled for Growth", "both", "not controlled for Growth"))
###*****************************


###*****************************
# change mrna to mRNA
joined_data$dataType<-replace_fun(input_vector = as.vector(joined_data$dataType), 
                                  initialVal = "mrna",
                                  finalVal = "mRNA")
###*****************************


###*****************************
# Figure
fig01<-ggplot(joined_data, aes(x = test_for)) +
  geom_bar(aes(fill = exist_in), position = "fill") +
  facet_grid(growthPhase ~ dataType) +
  scale_fill_discrete(name="Analysis",
                      labels=c("w/ doubling time", "either w/ or w/out \ndoubling time", "w/out doubling time"))+
  theme_bw()+
  xlab("Test Condition") + 
  ylab("Fraction of differentially expressed genes") +
  scale_y_continuous(expand = c(0, 0))+
  theme(panel.margin.y = unit(2, "lines"),
        panel.grid.minor.x = element_blank(),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.key.height=unit(2,"line"),
        #theme(legend.key.width=unit(5,"line")) +
        legend.background = element_rect(fill=alpha(colour = "white",
                                                    alpha = .4)))

print(fig01)
###*****************************


###*****************************
# save figure
cowplot::save_plot(filename = paste0("../c_figures/difference_rtw_GrowthControl",".pdf"),
                   plot = fig01, 
                   ncol = 2, nrow = 2, 
                   base_height = 3.075,base_width = 5.5,
                   units = "in",useDingbats=FALSE)
###*****************************


###*****************************
# Generate 2 data frames related with change
joined_data %>%
  dplyr::filter(dataType=="protein") %>%
  dplyr::filter(exist_in=="controlled for Growth") %>%
  dplyr::filter(test_for %in% c("Gly", "Glc", "Lac")) %>%
  dplyr::filter(growthPhase == "Exp") %>%
  .$genes %>% as.vector(.) %>% sort(.) %>%unique(.)-> changedExp

# write.csv(x = changedExp, file = "../c_results/changedExp.csv")

joined_data %>%
  dplyr::filter(dataType=="protein") %>%
  dplyr::filter(exist_in=="controlled for Growth") %>%
  dplyr::filter(test_for %in% c("Gly", "Glc", "Lac"))%>%
  dplyr::filter(growthPhase == "Sta") %>%
  .$genes %>% as.vector(.) %>% sort(.) %>%unique(.)-> changedSta

# write.csv(x = changedSta, file = "../c_results/changedSta.csv")
###*****************************

###*****************************
joined_data %>%
  dplyr::filter(dataType=="protein") %>%
  dplyr::filter(exist_in=="controlled for Growth") %>%
  dplyr::select(-uncontrolled, -controlled) %>%
  dplyr::filter(test_for %in% c("Gly", "Glc", "Lac")) %>%
  dplyr::filter(growthPhase %in% c("Exp", "Sta")) -> changedGenes_ExpSta

write.csv(x = changedGenes_ExpSta, file = "../c_results/changed_protein_carbonSource_ExpSta.csv")
###*****************************


###*****************************
# Prepeare dictionaries
dictionaryEz=read.csv(file="../generateDictionary/protein_tidy_eColi_ez.csv",row.names = 1)
dictionaryEz%>%dplyr::rename("gene_name"=From, "ez_gene_id"=To)->dictionaryEz
###*****************************

###*****************************
dictionaryEz%>%
  dplyr::filter(gene_name %in% changedExp) %>%
  .$ez_gene_id %>% as.vector(.) %>% unique(.) %>% sort(.)->changedExp_ez

dictionaryEz%>%
  dplyr::filter(gene_name %in% changedSta) %>%
  .$ez_gene_id %>% as.vector(.) %>% unique(.) %>% sort(.)->changedSta_ez
###*****************************


###*****************************
# Connect to david for analyse #EXP
david_d<-DAVIDWebService(email="umut.caglar@utexas.edu",
                         url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
#result<-addList(david_d, davidInputData, idType="OFFICIAL_GENE_SYMBOL", listName="testList", listType="Gene")
result<-addList(david_d, changedExp_ez, idType="ENTREZ_GENE_ID", listName="testList", listType="Gene")
###*****************************


###*****************************
# Set species and backround
selectedSpecie="Escherichia coli str. K-12 substr. MG1655"
backgroundLocation=grep(selectedSpecie,RDAVIDWebService::getBackgroundListNames(david_d))
specieLocation=grep(selectedSpecie,RDAVIDWebService::getSpecieNames(david_d))
setCurrentSpecies(object=david_d, species=specieLocation);
setCurrentBackgroundPosition(object=david_d,position=backgroundLocation)
###*****************************


###*****************************
# KEGG TEST
setAnnotationCategories(david_d, c("KEGG_PATHWAY"))

DavidChanged_exp_kegg<- as.data.frame(getFunctionalAnnotationChart(object=david_d,  threshold=1, count=0L))
DavidChanged_exp_kegg %>% dplyr::mutate(phase="Exp", name="changed", test="Kegg") -> DavidChanged_exp_kegg
# write.csv(x =DavidChanged_exp_kegg ,file = paste0("../c_results/","changedExp_DAVID_kegg",".csv"))
###*****************************


###*****************************
# MF NEW TEST
setAnnotationCategories(david_d, c("GOTERM_MF_ALL"))

DavidChanged_exp_mf<- as.data.frame(getFunctionalAnnotationChart(object=david_d,  threshold=1, count=0L))
DavidChanged_exp_mf %>% dplyr::mutate(phase="Exp", name="changed", test="Mf") -> DavidChanged_exp_mf
# write.csv(x =DavidChanged_exp_mf ,file = paste0("../c_results/","changedExp_DAVID_MF",".csv"))
###*****************************



###*****************************
# Connect to david for analyse #STA
david_d<-DAVIDWebService(email="umut.caglar@utexas.edu",
                         url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
#result<-addList(david_d, davidInputData, idType="OFFICIAL_GENE_SYMBOL", listName="testList", listType="Gene")
result<-addList(david_d, changedSta_ez, idType="ENTREZ_GENE_ID", listName="testList", listType="Gene")
###*****************************


###*****************************
# Set species and backround
selectedSpecie="Escherichia coli str. K-12 substr. MG1655"
backgroundLocation=grep(selectedSpecie,RDAVIDWebService::getBackgroundListNames(david_d))
specieLocation=grep(selectedSpecie,RDAVIDWebService::getSpecieNames(david_d))
setCurrentSpecies(object=david_d, species=specieLocation);
setCurrentBackgroundPosition(object=david_d,position=backgroundLocation)
###*****************************


###*****************************
# KEGG TEST
setAnnotationCategories(david_d, c("KEGG_PATHWAY"))

DavidChanged_sta_kegg<- as.data.frame(getFunctionalAnnotationChart(object=david_d,  threshold=1, count=0L))
DavidChanged_sta_kegg %>% dplyr::mutate(phase="Sta", name="changed", test="Kegg") -> DavidChanged_sta_kegg
# write.csv(x =DavidChanged_sta_kegg ,file = paste0("../c_results/","changedSta_DAVID_KEGG",".csv"))
###*****************************


###*****************************
# MF NEW TEST
setAnnotationCategories(david_d, c("GOTERM_MF_ALL"))

DavidChanged_sta_mf<- as.data.frame(getFunctionalAnnotationChart(object=david_d,  threshold=1, count=0L))
DavidChanged_sta_mf %>% dplyr::mutate(phase="Sta", 
                                     name="changed",
                                     test="Mf") -> DavidChanged_sta_mf
# write.csv(x =DavidChanged_sta_mf ,file = paste0("../c_results/","changedSta_DAVID_MF",".csv"))
###*****************************


###*****************************
# Combine all outputs of kegg
dplyr::bind_rows(DavidChanged_exp_kegg,DavidChanged_exp_mf,DavidChanged_sta_kegg,DavidChanged_sta_mf)->DavidChanged
DavidChanged %>% dplyr::filter(FDR < 0.05) -> DavidChanged
write.csv(x =DavidChanged ,file = paste0("../c_results/","changed_DAVID_P05",".csv"))
###*****************************