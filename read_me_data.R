library(Quandl)
library(dplyr)
library(purrr)
library(readabs)
library(lubridate)

# metaData(ts)
# Quandl.api_key("e2YZKez7xNnMZs3CYT8_")
# data <- Quandl.datatable("AUSBS/D")
# Quandl.datatable('AUSBS/D',date='1994-11-15', series_id='A85002072C')
# 
# China_GDPPC <- Quandl("WWDI/CHN_NY_GDP_PCAP_KN", type = 'xts')
# 
# Quandl.datatable('AUSBS/D',date='1994-11-15', series_id='A85002072C')
# 
# d <- read_api_dataflows()
# 
# wpi <- read_abs("6345.0")
# read_api("ABS_C16_T10_SA")
# 
# read_api_datastructure("ABORIGINAL_ID_POP_PROJ")

#read macro data from abs using series id
ue <- as.xts(read_abs(series_id = "A84423050A")[c("date","value")])
#collapsing monthly data to quarterly data
ue <- ue[index(ue) == floor_date(as.Date(as.yearqtr(index(ue)), frac = 1), unit = "month")]

gdp <- as.xts(read_abs(series_id = "A2304402X")[c("date","value")])
cpi <- as.xts(read_abs(series_id = "A2325846C")[c("date","value")])

au_me_data <- merge.xts(ue,gdp, cpi)
colnames(au_me_data) <- c("ue", "gdp", "cpi")

#take 1 year change
au_me_data_4qtr_chg <- na.omit(au_me_data/lag(au_me_data, 4) - 1)
colnames(au_me_data_4qtr_chg) <- paste0(colnames(au_me_data_4qtr_chg),"_4qtr_chg")

#take 2 year change
au_me_data_8qtr_chg <- na.omit(au_me_data/lag(au_me_data, 8) - 1)
colnames(au_me_data_8qtr_chg) <- paste0(colnames(au_me_data_8qtr_chg),"_8qtr_chg")

#combine log of unemployment rate with gdp, cpi 8 quarter changes
au_me_data_log_8qtr_chg <- na.omit(merge(log(au_me_data[,1]),au_me_data_8qtr_chg[,c(2:3)]))
colnames(au_me_data_log_8qtr_chg)[1] <- "log_ue"

