### This script processes drug data from DrugBank ###

## To run this script, you need two files from DrugBank that can be downloaded only after you obtain a license:
## After you get the appropriate license from DrugBank, you can download the necessary file from here:
## - Drug targets for all drugs in DrugBank: https://go.drugbank.com/releases/latest#protein-identifiers --> "Drug Target Identifiers" --> Drug Group "All"
## - Link between DrugBank IDs and drug names: https://go.drugbank.com/releases/latest#external-links --> "External Drug Links" --> Drug Group "All"
## NOTE: download and place the files in the "raw_data" folder

library(data.table)
library(dplyr)
library(stringr)
library(reshape2)

## load data
db_drug_targets = fread("raw_data/all.csv")
db_links = fread("raw_data/drug links.csv")

## each row --> a drug-gene target pair
db_drug_targets = db_drug_targets %>%
  filter(Species == "Humans") %>%
  dplyr::select(drug_target = "Gene Name", "Drug IDs", V14:V149) %>% 
  distinct()
db_drug_targets = melt(db_drug_targets, "drug_target", colnames(db_drug_targets)[2:ncol(db_drug_targets)])
db_drug_targets = db_drug_targets %>% 
  dplyr::select(drug_target, db_id = value) %>%
  distinct()
db_drug_targets = db_drug_targets[-which(db_drug_targets$db_id == ""), ] ; rownames(db_drug_targets) = NULL
db_drug_targets = db_drug_targets[-which(db_drug_targets$drug_target == ""), ] ; rownames(db_drug_targets) = NULL
db_drug_targets = db_drug_targets %>% 
  arrange(db_id)
length(unique(db_drug_targets$db_id)) # 5,800 unique DB_IDs

## add drug name
# keep needed information
db_links = db_links %>% dplyr::select(db_id = "DrugBank ID", drug_name = "Name") %>% distinct()

# annotate with drug names
db_drug_targets = left_join(db_drug_targets, db_links, by = "db_id")

fwrite(db_drug_targets, "processed_data/drugbank_all_drugs_known_targets.txt", sep = "\t", row.names = FALSE)
