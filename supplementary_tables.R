## This script creates the all the supplementary excel files
## You should run it after you have run all the other scripts in the repository as it depends on files created from them

library(dplyr)
library(data.table)
library(reshape2)
library(openxlsx)

## -- Candidate drugs per complex disease -- ##

md_cd_comorbidity = fread("processed_data/md_cd_comorbidities.txt") %>%
  dplyr::select(-comorbidity)
md_genes = fread("processed_data/md_genes.txt")

cd_categories = fread("raw_data/complex_disease_category.txt")

drug_targets = fread("processed_data/drugbank_all_drugs_known_targets.txt")
md_genes_drugs = left_join(md_genes, drug_targets[, 1:2], by = c("causal_gene" = "drug_target"))
md_genes_drugs = na.omit(md_genes_drugs) ; rownames(md_genes_drugs) = NULL
md_genes_drugs = md_genes_drugs %>% distinct()

cd_candidate_drugs = md_cd_comorbidity %>%
  left_join(cd_categories, by = "complex_disease") %>%
  left_join(md_genes, by = "mendelian_disease") %>%
  arrange(complex_disease, mendelian_disease, causal_gene) %>%
  distinct() %>%
  group_by(complex_disease, mendelian_disease) %>%
  mutate(md_causal_genes = paste0(unique(causal_gene), collapse = ",")) %>%
  ungroup() %>%
  dplyr::select(-causal_gene) %>%
  distinct() %>%
  left_join(md_genes_drugs[, c(1,3)], by = "mendelian_disease") %>%
  arrange(complex_disease, mendelian_disease, db_id) %>%
  group_by(complex_disease, mendelian_disease) %>%
  mutate(candidate_drugs = paste0(unique(db_id), collapse = ",")) %>%
  ungroup() %>%
  dplyr::select(complex_disease, complex_disease_category = disease_category, mendelian_disease, md_causal_genes, candidate_drugs) %>%
  distinct()
cd_candidate_drugs$candidate_drugs = ifelse(cd_candidate_drugs$candidate_drugs == "NA", NA, cd_candidate_drugs$candidate_drugs)

write.xlsx(cd_candidate_drugs, "supplementary_tables/S1_complex_diseases_candidate_drugs.xlsx", overwrite = TRUE)

## -- Ranked list of candidate drugs per complex disease -- ##
# run logistic regression to gain predicted probabilities for each drug-disease pair
md_cd_comorbidities = fread("processed_data/md_cd_comorbidities.txt") %>% 
  dplyr::select(mendelian_disease, complex_disease) %>%
  arrange(mendelian_disease)
# comorbidity matrix
md_cd_comorbidities_matrix = as.data.frame.matrix(table(md_cd_comorbidities[, 2:1]))
md_cd_comorbidities_matrix = md_cd_comorbidities_matrix %>% mutate(complex_disease = rownames(md_cd_comorbidities_matrix), .before = Achromatopsia)
rownames(md_cd_comorbidities_matrix) = NULL
# md genes
md_genes = fread("processed_data/md_genes.txt")
# drug - targets
db_drug_targets = fread("processed_data/drugbank_all_drugs_known_targets.txt") %>%
  dplyr::select(db_id, drug_target) %>% 
  distinct()
# number of targets per drug
drugs_nr_targets = fread("processed_data/drugbank_all_drugs_known_targets.txt") %>%
  dplyr::select(db_id, drug_target) %>%
  group_by(db_id) %>% 
  mutate(total_targets = length(drug_target)) %>%
  dplyr::select(db_id, total_targets) %>% 
  distinct() %>%
  ungroup()
# complex disease categories
complex_disease_categories = fread("raw_data/complex_disease_category.txt")
# indicated/investigated drugs
investigated_indicated_drugs = fread("processed_data/drugs_inv_ind_per_disease.txt")
investigated_indicated_drugs = reshape2::melt(investigated_indicated_drugs, "drugbank_id", colnames(investigated_indicated_drugs)[2:ncol(investigated_indicated_drugs)])
investigated_indicated_drugs = na.omit(investigated_indicated_drugs)
rownames(investigated_indicated_drugs) = NULL
investigated_indicated_drugs = investigated_indicated_drugs %>% dplyr::select(drugbank_id, complex_disease = variable) %>% distinct()
investigated_indicated_drugs$indicated_investigated = 1
# logistic regression to get probabilities
unique_cd = unique(md_cd_comorbidities$complex_disease)
unique_md = unique(md_cd_comorbidities$mendelian_disease)
# drugs targeting md genes
drugs_targeting_md_genes = left_join(md_genes, db_drug_targets, by = c("causal_gene" = "drug_target")) %>%
  na.omit() %>%
  dplyr::select(db_id) %>%
  distinct()
# logistic regression input - each row is a drug-complex disease pair 
log_input_all = data.frame(db_id = rep(drugs_targeting_md_genes$db_id, each = length(unique_cd)))
log_input_all$complex_disease = unique_cd
# information about disease category
log_input_all = left_join(log_input_all, complex_disease_categories, by = "complex_disease")
# information about number of targets for each drug
log_input_all = left_join(log_input_all, drugs_nr_targets, by = "db_id")  
# information about investigated/indicated drugs within each drug-cancer pair
log_input_all = left_join(log_input_all, investigated_indicated_drugs, by = c("complex_disease" = "complex_disease", "db_id" = "drugbank_id"))
log_input_all$indicated_investigated = ifelse(is.na(log_input_all$indicated_investigated), 0, 1)
# information about comorbidity
x = left_join(md_cd_comorbidities, md_genes, by = "mendelian_disease")
x = na.omit(x) ; rownames(x) = NULL
x = x %>% dplyr::select(complex_disease, causal_gene) %>% distinct
x = left_join(x, db_drug_targets, by = c("causal_gene" = "drug_target"))
x = na.omit(x) ; rownames(x) = NULL
x = x %>% dplyr::select(-causal_gene) %>% distinct()
x$comorbidity = 1
log_input_all = left_join(log_input_all, x, by = c("complex_disease", "db_id"))
log_input_all$comorbidity = ifelse(is.na(log_input_all$comorbidity), 0, 1)
# logistic regression
log_reg_all = glm(indicated_investigated ~ total_targets + disease_category + comorbidity,
                  data = log_input_all, 
                  family = binomial())
# predicted probabilities for each drug-disease pair
drug_disease_ranked = data.frame(complex_disease = log_input_all$complex_disease, 
                                 drugbank_id = log_input_all$db_id,
                                 predicted_prob = predict(log_reg_all, type = "response"))
drug_disease_ranked = drug_disease_ranked %>% arrange(complex_disease, desc(predicted_prob))
rownames(drug_disease_ranked) = NULL
drug_disease_ranked = drug_disease_ranked %>%
  group_by(complex_disease) %>%
  mutate(drug_rank2 = match(predicted_prob * -1, sort(unique(predicted_prob) * -1))) %>% # multiply by -1 to have, for each disease, as 1st the drug with the highest predicted probability
  ungroup() %>%
  dplyr::select(complex_disease, drugbank_id, rank = drug_rank2) %>%
  distinct()
drug_disease_ranked = left_join(drug_disease_ranked, investigated_indicated_drugs, by = c("complex_disease" = "complex_disease",
                                                                                          "drugbank_id" = "drugbank_id"))
drug_disease_ranked$indicated_investigated = ifelse(is.na(drug_disease_ranked$indicated_investigated), 0, 1)

write.xlsx(drug_disease_ranked, "supplementary_tables/S2_complex_diseases_candidate_drugs_ranked.xlsx", overwrite = TRUE)

## -- Genetic similarity metrics per Mendelian disease - cancer pair -- ##

genetic_overlap = fread("processed_data/md_cancers_genetic_overlap.txt")
coexpression = fread("processed_data/md_cancers_coexpression.txt")

combined = full_join(genetic_overlap, coexpression, by = c("cancer", "mendelian_disease", "comorbidity"))
colnames(combined) = c("mendelian_disease", "cancer", "pvalue_onesided_geneoverlap", "comorbidity", "genetic_similarity_geneoverlap", "adj_pvalue_coexpression", "genetic_similarity_coexpression")

combined = combined %>%
  dplyr::select(cancer, mendelian_disease, gene_overlap_pvalue = pvalue_onesided_geneoverlap, coexpression_pvalue = adj_pvalue_coexpression) %>%
  distinct()

write.xlsx(combined, "supplementary_tables/S3_cancers_mendelian_diseases_genetic_similarity.xlsx", overwrite = TRUE)

## -- Mendelian diseases and cancers - candidate drugs across different levels of support -- ##

combined$genetic_similarity = ifelse(combined$gene_overlap_pvalue < 0.05 | combined$coexpression_pvalue < 0.05, 1, 0)

combined = combined %>%
  dplyr::select(cancer, mendelian_disease, genetic_similarity) %>%
  distinct()

# add information about comorbidity
md_cd_comorbidity$comorbidity = 1
combined = left_join(combined, md_cd_comorbidity, by = c("cancer" = "complex_disease", "mendelian_disease"))
combined$comorbidity = ifelse(is.na(combined$comorbidity), 0, 1)

# add mendelian disease genes and candidate drugs
combined = combined %>%
  left_join(md_genes, by = "mendelian_disease") %>%
  arrange(cancer, mendelian_disease, causal_gene) %>%
  distinct() %>%
  group_by(cancer, mendelian_disease) %>%
  mutate(md_causal_genes = paste0(unique(causal_gene), collapse = ",")) %>%
  ungroup() %>%
  dplyr::select(-causal_gene) %>%
  distinct() %>%
  left_join(md_genes_drugs[, c(1,3)], by = "mendelian_disease") %>%
  arrange(cancer, mendelian_disease, db_id) %>%
  group_by(cancer, mendelian_disease) %>%
  mutate(candidate_drugs = paste0(unique(db_id), collapse = ",")) %>%
  ungroup() %>%
  dplyr::select(cancer, mendelian_disease, md_causal_genes, comorbidity, genetic_similarity, candidate_drugs) %>%
  distinct()
combined$candidate_drugs = ifelse(combined$candidate_drugs == "NA", NA, combined$candidate_drugs)
# Check genetic similarity between Familial Dysautonomia and cancers (genetic overlap measurement only)
genetic_overlap %>% filter(mendelian_disease == "Familial Dysautonomia") # no genetic similarity --> add this information to the "combined" dataframe
combined$genetic_similarity = ifelse(combined$mendelian_disease == "Familial Dysautonomia", 0, combined$genetic_similarity)

write.xlsx(combined, "supplementary_tables/S4_cancers_candidate_drugs.xlsx", overwrite = TRUE)
