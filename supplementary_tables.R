## This script creates the all the supplementary excel files
## You should run it after you have run all the other scripts in the repository as it depends on files created from them

library(dplyr)
library(data.table)
library(openxlsx)

## -- Candidate drugs per complex disease -- ##

md_cd_comorbidity = fread("processed_data/md_cd_comorbidities.txt") %>%
  dplyr::select(-comorbidity)
md_genes = fread("processed_data/md_genes.txt")

cd_categories = fread("processed_data/complex_disease_category.txt")

drug_targets = fread("processed_data/drugbank_all_drugs_known_targets.txt")
md_genes_drugs = left_join(md_genes, drug_targets[, 1:2], by = c("causal_gene" = "drug_target"))
md_genes_drugs = na.omit(md_genes_drugs) ; rownames(md_genes_drugs) = NULL
md_genes_drugs = md_genes_drugs %>% distinct()

cd_candidate_drugs = md_cd_comorbidity %>%
  left_join(cd_categories, by = "complex_disease") %>%
  left_join(md_genes_drugs, by = "mendelian_disease") %>%
  distinct() %>%
  na.omit() %>%
  arrange(complex_disease, mendelian_disease, db_id) %>%
  group_by(complex_disease, mendelian_disease) %>%
  mutate(md_causal_genes = paste0(unique(causal_gene), collapse = ","),
         candidate_drugs = paste0(db_id, collapse = ",")) %>%
  ungroup() %>%
  dplyr::select(complex_disease, mendelian_disease, md_causal_genes, candidate_drugs, complex_disease_category = disease_category) %>%
  distinct()

write.xlsx(cd_candidate_drugs, "supplementary_tables/S1_complex_diseases_candidate_drugs.xlsx", overwrite = TRUE)

## -- Genetic similarity metrics per Mendelian disease - cancer pair -- ##

genetic_overlap = fread("processed_data/md_cancers_genetic_overlap.txt")
coexpression = fread("processed_data/md_cancers_coexpression.txt")

combined = full_join(genetic_overlap, coexpression, by = c("cancer", "mendelian_disease", "comorbidity"))
colnames(combined) = c("mendelian_disease", "cancer", "pvalue_onesided_geneoverlap", "comorbidity", "genetic_similarity_geneoverlap", "adj_pvalue_coexpression", "genetic_similarity_coexpression")

combined = combined %>%
  dplyr::select(cancer, mendelian_disease, gene_overlap_pvalue = pvalue_onesided_geneoverlap, coexpression_pvalue = adj_pvalue_coexpression) %>%
  distinct()

write.xlsx(combined, "supplementary_tables/S2_cancers_mendelian_diseases_genetic_similarity.xlsx", overwrite = TRUE)

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

write.xlsx(combined, "supplementary_tables/S3_cancers_candidate_drugs.xlsx", overwrite = TRUE)
