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

## -- visualizations -- ##

## number of drugs per Mendelian disease
md_genes = fread("processed_data/md_genes.txt")
md_drugs = left_join(md_genes, db_drug_targets[, 1:2], by = c("causal_gene" = "drug_target"))
md_drugs = md_drugs %>%
  dplyr::select(mendelian_disease, db_id) %>%
  distinct()

# calculate number of drugs per Mendelian disease
md_nr_drugs = md_drugs %>% 
  group_by(mendelian_disease) %>%
  mutate(nr_drugs = if_else(!is.na(db_id), length(db_id), 0)) %>%
  ungroup() %>%
  dplyr::select(mendelian_disease, nr_drugs) %>%
  group_by(mendelian_disease) %>%
  mutate(nr_drugs = max(nr_drugs)) %>%
  ungroup() %>%
  distinct() %>%
  arrange(nr_drugs)
md_nr_drugs$mendelian_disease = factor(md_nr_drugs$mendelian_disease, levels = md_nr_drugs$mendelian_disease, labels = md_nr_drugs$mendelian_disease)

ggplot(md_nr_drugs, aes(x = nr_drugs, y = mendelian_disease)) +
  geom_col(fill = "gray", color = "black") +
  xlab("Number of drugs") +
  ylab("") +
  scale_x_continuous(breaks = seq(0, 120, 5)) +
  theme(axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.3), 
        axis.title = element_text(angle = 0, hjust = 0.5,
                                  margin = margin(t = 3, unit = "cm"),
                                  size = 40, family = "Arial", colour = "black"),
        axis.text = element_text(angle = 0, hjust = 0.5, vjust = 0.5,
                                 margin = margin(l = 0.5, r = 0.2, unit = "cm"),
                                 size = 40, family = "Arial", colour = "black"))