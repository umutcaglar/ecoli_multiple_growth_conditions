# generating raw file list

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
{setwd(paste0("/Users/umut/Desktop/responses/"))} # mac computer
###*****************************


###*****************************
# Data tracking
require("dplyr")
require("tidyr")
require("broom")
require("stringr")
###*****************************


###*****************************
# Load data
dataFile = read.csv(file = "rawDataList.csv")
colnames(dataFile)<-"V1"
###*****************************


###*****************************
# Seperate to multiple columns
dataFile %>%
  tibble::rownames_to_column()%>%
  dplyr::group_by(rowname)%>%
  dplyr::mutate(short=gsub(pattern = "*.* ", replacement = "", x = V1))->dataFile
###*****************************


###*****************************
dataFile %>%
  dplyr::mutate(MURI=gsub(pattern="_SA*.*", replacement = " ", x = short)) %>%
  dplyr::mutate(L=stringr::str_extract(string = short, pattern = "L...")) %>%
  dplyr::mutate(R=gsub(pattern = "_", replacement = "", 
                       stringr::str_extract(string = short, pattern = "_R._"))) %>%
  dplyr::select(-V1)->dataFile
###*****************************


###*****************************
# Group data
dataFile %>%
  dplyr::group_by(MURI, R) %>%
  dplyr::mutate(L_count=paste0(R,"_L",1:n())) %>%
  dplyr::select(-rowname)->dataFile
###*****************************


###*****************************
# Organize DF
dataFile %>%
  dplyr::group_by() %>%
  dplyr::select(-L,-R) %>%
  tidyr::spread(data = ., key = L_count, value = short)->dataFile 
###*****************************


###*****************************
# generate numbers to reoder data
dataFile %>%
  dplyr::mutate(MURI_num=as.numeric(gsub(pattern = "*.*_",replacement = "",x = MURI))) %>%
  dplyr::arrange(MURI_num)->dataFile 
###*****************************


###*****************************
# save it
write.csv(x = dataFile, file = "organized_raw_files.csv", na = "")
###*****************************









###*****************************
# Load data
dataFile2 = read.csv(file = "checklist.csv")
colnames(dataFile2)<-"V1"
###*****************************


###*****************************
# Seperate to multiple columns
dataFile2 %>%
  tibble::rownames_to_column()%>%
  dplyr::group_by(rowname)%>%
  dplyr::mutate(short = gsub(pattern = "*.* ", replacement = "", x = V1)) %>%
  dplyr::mutate(md5sum = stringr::str_extract(string = V1, pattern = " *.* ") )->dataFile2
###*****************************


###*****************************
dataFile2 %>%
  dplyr::mutate(MURI=gsub(pattern="_SA*.*", replacement = " ", x = short)) %>%
  dplyr::mutate(L=stringr::str_extract(string = short, pattern = "L...")) %>%
  dplyr::mutate(R=gsub(pattern = "_", replacement = "", 
                       stringr::str_extract(string = short, pattern = "_R._"))) %>%
  dplyr::select(-V1)->dataFile2
###*****************************



###*****************************
# generate numbers to reoder data
dataFile2 %>%
  dplyr::mutate(MURI_num=as.numeric(gsub(pattern = "*.*_",replacement = "",x = MURI))) %>%
  dplyr::arrange(MURI_num, R, L)->dataFile2
###*****************************


###*****************************
# save it
write.csv(x = dataFile2, file = "organized_md5.csv", na = "")
###*****************************
