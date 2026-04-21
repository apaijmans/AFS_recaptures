# #~~ Postprocessing: use custom written function to process recaptures result into excel format

library(tidyverse)
library(here)

# Location recaptures result
# Change accordingly:
file_name <- "results_pup_pup_r159-168_10_02_26"
file_path <- here("Recapture_analysis_server", "Data", "processed")

# Location msat file
msat_name <- paste0("all_msat_genotypes_raw_r159-168.xlsx")
msat_path <- here("Recapture_analysis_server", "Data", "raw")

# Load recaptures result
load(paste0(file_path, "/", file_name, ".Rdata"))

# Source function to converte recapture result to excel
source("Rcode/exportResults.R")

# Wrangle msat data for function
msat_alleles <- readxl::read_excel(paste0(msat_path, "/", msat_name), guess_max = 16000) %>%
  filter(Pup==1 | Female==1 | Male ==1) %>%
  #select(1, 16:93)
  select(SampleID,Pup,Female,Male,PlateNumber,PlateLocation, Matches, Pv9.a:Mang36.b) %>%
  mutate(SampleID = sub("\\?", "NOT", SampleID)) %>%
  mutate(SampleID = sub(" \\(REPEAT\\)", "_REPEAT", SampleID)) %>%
  mutate(SampleID = sub(" \\(POSITIVE CONTROL\\)", "_POSITIVECONTROL", SampleID)) %>%
  mutate(SampleID = sub(" \\(", "_", SampleID)) %>%
  mutate(SampleID = sub("\\)", "", SampleID)) %>%
  mutate(SampleID = gsub(" ", "", SampleID))

# Fix 1 sample ID that causes trouble
msat_alleles[which(msat_alleles$SampleID=="M1"),1] <- "dummy_M1"

# Run function and save recapture result as excel file
exportResults(results, msat_alleles, paste0(here("Recapture_analysis_local", "Data", "processed"), "/", file_name))
