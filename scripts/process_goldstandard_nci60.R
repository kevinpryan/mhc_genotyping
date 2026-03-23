library(tidyverse)
source("scripts/functions/ggroup_mapper.R")
test <- readRDS("downloads/pub/adams_2005/hla_calls.rds")
test$`Cell Line`
colnames(test)[1] <- "sample_id"
# readRDS("testreadRDS("downloads/pub/adams_2005/hla_calls.rds") %>%
#   rename(`Cell line` = "sample_id") #%>%
  #select(-ID, -Tissue) %>%
  #allele_map %>%
  #saveRDS(str_glue("data/gold_standard_nci60.rds"))

test %>% 
  dplyr::select(-ID, -Tissue) %>% 
  allele_map %>% 
  saveRDS(str_glue("data/gold_standard_nci60.rds"))
  
nci60.gold.standard <- readRDS("data/gold_standard_nci60.rds")
