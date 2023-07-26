

rm(list=ls())

library(tidyverse)
library(zoo)


# processRawData.R
#
# This script loads data from /data/raw/, processes it, and saves it to the /data/ folder.
#
# Before running, ensure that wastewater data "ww_data_all.csv" and "sites.csv", ....
#  ... as well as reported cases "covid-case-counts-moh.csv" are saved in the /data/raw/ folder.
# 
# Also ensure the working directory is set to the home directory of the repo.
#
# The three required files in can be downloaded from:
#   - https://github.com/ESR-NZ/covid_in_wastewater
#   - https://github.com/minhealthnz/nz-covid-data


# ------------- Process wastewater data ------------- 

df_public = read.csv("data/raw/ww_data_all.csv") %>%
  mutate(date=as.Date(Collected)) %>%
  left_join(select(read.csv("data/raw/sites.csv"), SampleLocation, Population), by = "SampleLocation") %>%
  rename(sample_location=SampleLocation, population=Population)

df_public_daily = df_public %>%
  select(date, sample_location, sars_gcl, population, copies_per_day_per_person) %>%
  mutate(copies_per_day_per_person=as.numeric(copies_per_day_per_person)) %>%
  mutate(implied_total_copies_per_day = copies_per_day_per_person * population) %>%
  group_by(date) %>%
  summarise(total_daily_copies = sum(implied_total_copies_per_day),
            popn_in_catchment = sum(population)) %>%
  mutate(copies_per_person = total_daily_copies/popn_in_catchment,
         wwDataIsValid = !is.na(copies_per_person) & !is.na(popn_in_catchment)) %>%
  rename(nGenomeCopies=copies_per_person, nPopInCatchment=popn_in_catchment)

write.csv(df_public_daily, "data/wastewater.csv", row.names=FALSE)


# ------------- Reported cases ------------- 

df_cases = read.csv("data/raw/covid-case-counts-moh.csv") %>%
  mutate(date = as.Date(Report.Date)) %>%
  group_by(date) %>%
  summarise(nCases=sum(Number.of.cases.reported),
            nBorderCases=sum(Number.of.cases.reported[District=="At the border"]),
            nLocalCases=sum(Number.of.cases.reported[District!="At the border"])) %>%
  arrange(date) %>%
  mutate(nCasesMovingAverage=rollmean(nCases, 7, fill=NA, align="center"))


write.csv(df_cases, "data/cases.csv", row.names=FALSE)
