## This script creates the all the supplementary excel files
## You should run it after you have run all the other scripts in the repository as it depends on files created from them

library(dplyr)
library(data.table)
library(openxlsx)

## -- Complex diseases and Mendelian genes connected to them through comorbidity -- ##

md_cd_comorbidity = fread("processed_data/md_cd_comorbidities.txt")
md_genes = fread("processed_data/md_genes.txt")
md_cd_comorbidity_genes = left_join(md_cd_comorbidity, md_genes, by = "mendelian_disease")

cd_md_genes = md_cd_comorbidity_genes %>%
  dplyr::select(-comorbidity) %>%
  group_by(complex_disease) %>%
  mutate(nr_md_comorbidities = length(unique(mendelian_disease))) %>%
  mutate(md_comorbidities = paste0(unique(mendelian_disease), collapse = ",")) %>%
  dplyr::select(-mendelian_disease) %>%
  mutate(nr_md_genes = length(causal_gene)) %>%
  mutate(md_genes = paste0(unique(causal_gene), collapse = ",")) %>%
  ungroup() %>%
  dplyr::select(complex_disease, md_comorbidities, nr_md_comorbidities, md_genes, nr_md_genes) %>%
  distinct()

## add complex disease category
cd_categories = fread("processed_data/complex_disease_category.txt")
cd_md_genes = left_join(cd_md_genes, cd_categories, by = "complex_disease") %>%
  dplyr::relocate(disease_category, .after = "complex_disease")

write.xlsx(cd_md_genes, "supplementary_tables/S1_complex_diseases_md_comorbidities_genes.xlsx", overwrite = TRUE)

## -- Candidate drugs per complex disease -- ##

drug_targets = fread("processed_data/drugbank_all_drugs_known_targets.txt")
md_genes_drugs = left_join(md_genes, drug_targets[, 1:2], by = c("causal_gene" = "drug_target"))
md_genes_drugs = na.omit(md_genes_drugs) ; rownames(md_genes_drugs) = NULL
md_genes_drugs = md_genes_drugs %>%
  dplyr::select(-causal_gene) %>%
  distinct()
cd_candidate_drugs = left_join(md_cd_comorbidity, md_genes_drugs, by = "mendelian_disease") %>%
  dplyr::select(-comorbidity, -mendelian_disease) %>%
  distinct() %>%
  arrange(complex_disease, db_id) %>%
  group_by(complex_disease) %>%
  mutate(candidate_drugs = paste0(db_id, collapse = ","),
         nr_candidate_drugs = length(db_id)) %>%
  ungroup() %>%
  dplyr::select(complex_disease, nr_candidate_drugs, candidate_drugs) %>%
  distinct()

write.xlsx(cd_candidate_drugs, "supplementary_tables/S2_complex_diseases_candidate_drugs.xlsx", overwrite = TRUE)

## -- Genetic similarity metrics per Mendelian disease - cancer pair -- ##

genetic_overlap = fread("processed_data/md_cancers_genetic_overlap.txt")
coexpression = fread("processed_data/md_cancers_coexpression.txt")

combined = full_join(genetic_overlap, coexpression, by = c("cancer", "mendelian_disease", "comorbidity"))
colnames(combined) = c("mendelian_disease", "cancer", "pvalue_onesided_geneoverlap", "comorbidity", "genetic_similarity_geneoverlap", "adj_pvalue_coexpression", "genetic_similarity_coexpression")
# fix information for Familial Dysautonomia
combined[541:550, "genetic_similarity_coexpression"] = 0

combined = combined %>%
  dplyr::select(cancer, mendelian_disease, gene_overlap_pvalue = pvalue_onesided_geneoverlap, coexpression_pvalue = adj_pvalue_coexpression) %>%
  distinct() %>%
  group_by(cancer, mendelian_disease)

write.xlsx(combined, "supplementary_tables/S3_mendelian_diseases_cancers_genetic_similarity.xlsx", overwrite = TRUE)
