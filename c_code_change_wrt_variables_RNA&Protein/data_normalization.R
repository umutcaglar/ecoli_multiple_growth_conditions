## DeSeq Normalization Cleaned 03

# The aim of the code is to generate normalized data matrix.
# The work flow compses of four parts
# Pick up the samples
# pick up the rows
# calculate size factors
# do the normalization


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
require("Biobase") 
require("DESeq2")
require("dplyr")
require("tidyr")
###*****************************


###*****************************
#Load Functions
source("data_normalization_functions.R")
###*****************************


###*****************************
# @Title pick_data
# @Param data name could be 
#  "rna" (that installs all rna data including mRNA tRNA and RF)
#  "mrna" (that installs all mrna data)
#  "protein" (that installs all protein data)
#  "protein_wo_NA" (that installs protein data with removed NA rows)
# @Description installs necessary rawData and metaData
mainData=pick_data(data_name="mrna")
###*****************************


for(counter01 in 1:2)
{
  ###*****************************
  # Remove problematic data sets
  # problematic data sets are categorized in 3 different sets
  #@param problematic_set can be one of the three options.
  # Set 0. No filtering
  # set 1 all problematic rna related files
  # set 2 filtered problematic rna related files
  mainData=filter_bad_quality_data(dataInput = mainData, problematic_set = "set02")
  ###*****************************
  
  
  ###*****************************
  # @Title pick_experiment
  # @ParamexperimentVector Experiment vector can be composed of abbrevations of experiment names
  #   Stc="glucose_time_course", 
  #   Ytc="glycerol_time_course", 
  #   Nas="NaCl_stress", 
  #   Agr="lactate_growth", 
  #   Ngr="gluconate_growth", 
  #   Mgl="MgSO4_stress_low", 
  #   Mgh="MgSO4_stress_high"
  # @Description it filters the experiments from raw data and meta data 
  #      also adds remainng experiments to fileName
  mainData=pick_experiments(dataInput = mainData, experimentVector = c("allEx"))
  ###*****************************
  
  
  ###*****************************
  # @Title pick_carbonSource
  # @ParamexperimentVector Experiment vector can be composed of abbrevations of carbon sources
  #   "S"="glucose_time_course", 
  #   "Y"="glycerol_time_course", 
  #   "A="lactate_growth", 
  #   "N"="gluconate_growth", 
  # @Description it filters the carbon sources from raw data and meta data 
  #      also adds remainng carbon sources to fileName
  mainData=pick_carbonSource(dataInput = mainData, carbonSourceVector="Y")
  ###*****************************
  
  
  ###*****************************
  # @Title pick_MgLevel
  # @Param MgLevelVector MgLevel vector can be composed of abbrevations of MgLevel names
  #   "lowMg"="low Mg levels", 
  #   "midMg"="midium Mg levels", 
  #   "highMg"="high Mg levels"
  # @Description it filters the Mg levels from raw data and meta data 
  #      also adds remainng Mg levels to fileName
  mainData=pick_MgLevel(dataInput = mainData, MgLevelVector=c("allMg"))
  ###*****************************
  
  
  ###*****************************
  # @Title pick_NaLevel
  # @Param NaLevelVector NaLevel vector can be composed of abbrevations of NaLevel names
  #   "baseNa"="midium Na levels", 
  #   "highNa"="high Na levels"
  # @Description it filters the Na levels from raw data and meta data 
  #      also adds remainng Na levels to fileName
  mainData=pick_NaLevel(dataInput = mainData, NaLevelVector=c("baseNa"))
  ###*****************************
  
  
  ###*****************************
  # @Title pick_growthPhase
  # @Param growthPhaseVector growthPhase vector can be composed of abbrevations of growthPhase names
  #   "exponential" = "exponential phase"
  #   "stationary" = "stationary phase", 
  #   "late_stationary" = "late stationary phase"
  # @Description it filters the growth phase from raw data and meta data 
  #      also adds remainng chosen growth phase to fileName
  mainData=pick_growthPhase(dataInput = mainData, growthPhaseVector = c("allPhase"))
  ###*****************************
  
  
  ###*****************************
  # @Title filter_rows
  # @Param fitering method is the way tha we filter out the mon-significant rows
  #   methods are
  #   "noFilter" (basically keeps all the rows)
  #   "meanFilter" (filter out the rows whose mean is below threshold)
  #   "maxFilter" (filter out the rows whose max is below threshold)
  #   "sdFilter" (filter out the rows whose standard deviation is below threshold)
  # @param threshold the threshold related with relevant model. Default is 0.
  # @Description The function filter out the meaningless rows from data
  mainData=filter_rows(dataInput=mainData, filtering_method="noFilter")
  ###*****************************
}





