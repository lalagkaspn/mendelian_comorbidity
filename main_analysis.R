
library(data.table)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ComplexHeatmap)

## Set type for image compression based on operating system
## For macOS, X11 installation is required (link: https://www.xquartz.org/)
# if Windows OS
if (Sys.info()['sysname'] == "Windows") {
  type_compression = "windows"
} 
# if Linux
if (Sys.info()['sysname'] == "Linux") {
  type_compression = "windows"
}
# if macOS
if (Sys.info()['sysname'] == "Darwin") {
  type_compression = "cairo"
} 

#### --- all comorbidities --- ####

### observed ###

## load comorbidities
md_cd_comorbidities = fread("processed_data/md_cd_comorbidities.txt") %>% 
  dplyr::select(mendelian_disease, complex_disease) %>%
  arrange(mendelian_disease)

## comorbidity matrix
md_cd_comorbidities_matrix = as.data.frame.matrix(table(md_cd_comorbidities[, 2:1]))
md_cd_comorbidities_matrix = md_cd_comorbidities_matrix %>% mutate(complex_disease = rownames(md_cd_comorbidities_matrix), .before = Achromatopsia)
rownames(md_cd_comorbidities_matrix) = NULL

## load md genes
md_genes = fread("processed_data/md_genes.txt")

## load drug - targets
db_drug_targets = fread("processed_data/drugbank_all_drugs_known_targets.txt") %>%
  dplyr::select(db_id, drug_target) %>% 
  distinct()

## calculate number of targets per drug
drugs_nr_targets = fread("processed_data/drugbank_all_drugs_known_targets.txt") %>%
  dplyr::select(db_id, drug_target) %>%
  group_by(db_id) %>% 
  mutate(total_targets = length(drug_target)) %>%
  dplyr::select(db_id, total_targets) %>% 
  distinct() %>%
  ungroup()

## load complex disease categories
complex_disease_categories = fread("raw_data/complex_disease_category.txt")

## load indicated/investigated drugs
investigated_indicated_drugs = fread("processed_data/drugs_inv_ind_per_disease.txt")
investigated_indicated_drugs = melt(investigated_indicated_drugs, "drugbank_id", colnames(investigated_indicated_drugs)[2:ncol(investigated_indicated_drugs)])
investigated_indicated_drugs = na.omit(investigated_indicated_drugs)
rownames(investigated_indicated_drugs) = NULL
investigated_indicated_drugs = investigated_indicated_drugs %>% dplyr::select(drugbank_id, complex_disease = variable) %>% distinct()
investigated_indicated_drugs$indicated_investigated = 1

## run logistic regression - observed ##
# unique Mendelian and complex diseases
unique_cd = unique(md_cd_comorbidities$complex_disease)
unique_md = unique(md_cd_comorbidities$mendelian_disease)
# drugs targeting md genes
drugs_targeting_md_genes = left_join(md_genes, db_drug_targets, by = c("causal_gene" = "drug_target")) %>%
  na.omit() %>%
  dplyr::select(db_id) %>%
  distinct()
# create logistic regression input - each row is a drug-complex disease pair 
log_input_all = data.frame(db_id = rep(drugs_targeting_md_genes$db_id, each = length(unique_cd)))
log_input_all$complex_disease = unique_cd
# add information about disease category
log_input_all = left_join(log_input_all, complex_disease_categories, by = "complex_disease")
# add information about number of targets for each drug
log_input_all = left_join(log_input_all, drugs_nr_targets, by = "db_id")  
# add information about investigated/indicated drugs within each drug-cancer pair
log_input_all = left_join(log_input_all, investigated_indicated_drugs, by = c("complex_disease" = "complex_disease", "db_id" = "drugbank_id"))
log_input_all$indicated_investigated = ifelse(is.na(log_input_all$indicated_investigated), 0, 1)
# add information about comorbidity
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

all_obs_pvalue = summary(log_reg_all)$coefficients[8, 4]
all_obs_or = exp(summary(log_reg_all)$coefficients[8, 1])
all_obs_ci = exp(confint(log_reg_all)[8, ])

### permutations ###

## define function to shuffle complex disease comorbidities for each Mendelian disease
shuffle = function(comorbidity_matrix) {
  comorbidity_matrix_shuffled = comorbidity_matrix
  
  ## shuffling
  for (i in 2:ncol(comorbidity_matrix_shuffled)) {
    # Exclude first column - complex disease name
    comorbidity_matrix_shuffled[, i] = sample(comorbidity_matrix_shuffled[, i], replace = FALSE)
  }
  
  ## convert to original format
  comorbidity_matrix_shuffled = tidyr::gather(comorbidity_matrix_shuffled, key = "mendelian_disease", value = "comorbidity", -complex_disease) %>%
    dplyr::relocate(mendelian_disease, .before = complex_disease) %>% 
    filter(comorbidity == 1) %>% 
    dplyr::select(-comorbidity) %>%
    distinct()
  
  ## return
  return(comorbidity_matrix_shuffled)
}

# permutations (save odds ratio)
all_or_perm = c()
for (permutation in 1:1000) {
  
  ## shuffle comorbidity matrix
  md_cd_comorbidities_shuffled = shuffle(md_cd_comorbidities_matrix)

  ## logistic regression input - each row is a drug-complex disease pair 
  log_input_perm = data.frame(db_id = rep(drugs_targeting_md_genes$db_id, each = length(unique_cd)))
  log_input_perm$complex_disease = unique_cd
  # add information about disease category
  log_input_perm = left_join(log_input_perm, complex_disease_categories, by = "complex_disease")
  # add information about number of targets for each drug
  log_input_perm = left_join(log_input_perm, drugs_nr_targets, by = "db_id")  
  # add information about investigated/indicated drugs within each drug-cancer pair
  log_input_perm = left_join(log_input_perm, investigated_indicated_drugs, by = c("complex_disease" = "complex_disease", 
                                                                        "db_id" = "drugbank_id"))
  # if indicated_investigated is NA --> this drug is not indicated/investigated for this cancer
  log_input_perm$indicated_investigated = ifelse(is.na(log_input_perm$indicated_investigated), 0, 1)
  
  # add information about comorbidity
  x = left_join(md_cd_comorbidities_shuffled, md_genes, by = "mendelian_disease")
  x = na.omit(x) ; rownames(x) = NULL
  x = x %>% dplyr::select(complex_disease, causal_gene) %>% distinct
  x = left_join(x, db_drug_targets, by = c("causal_gene" = "drug_target"))
  x = na.omit(x) ; rownames(x) = NULL
  x = x %>% dplyr::select(-causal_gene) %>% distinct()
  x$comorbidity = 1
  log_input_perm = left_join(log_input_perm, x, by = c("complex_disease", "db_id"))
  log_input_perm$comorbidity = ifelse(is.na(log_input_perm$comorbidity), 0, 1)
  
  ## run logistic regression
  glm_fits_comorbidity_perm = glm(indicated_investigated ~ total_targets + disease_category + comorbidity,
                                  data = log_input_perm, 
                                  family = binomial())
  log_summary_comorbidity_perm = summary(glm_fits_comorbidity_perm)$coefficients
  
  ## update permutation vector
  all_or_perm = c(all_or_perm, log_summary_comorbidity_perm[8, 1])
  
  ## track progress
  cat(permutation, "\n")
}

# calculate permutation 5%, 50% and 95% quantiles
all_perm_ci = exp(quantile(all_or_perm, c(0.05, 0.50, 0.95)))

sum(all_obs_or <= exp(all_or_perm)) / 1000

#### --- per disease category --- ####

## observed ##
## unique complex diseases
cancers = complex_disease_categories %>% filter(disease_category == "cancer") ; cancers = unique(cancers$complex_disease)
cardiovascular = complex_disease_categories %>% filter(disease_category == "cardiovascular") ; cardiovascular = unique(cardiovascular$complex_disease)
hormonal = complex_disease_categories %>% filter(disease_category == "hormonal") ; hormonal = unique(hormonal$complex_disease)
immune = complex_disease_categories %>% filter(disease_category == "immune") ; immune = unique(immune$complex_disease)
neurological = complex_disease_categories %>% filter(disease_category == "neurological") ; neurological = unique(neurological$complex_disease)
ophthalmological = complex_disease_categories %>% filter(disease_category == "ophthalmological") ; ophthalmological = unique(ophthalmological$complex_disease)

## define function to run logistic regression per disease category
logistic_regression_per_dis_cat = function(diseases_in_category, comorbidity_dataframe, nr_drug_targets) {
  
  dis_cat_specific_comorbidity = comorbidity_dataframe %>% 
    filter(complex_disease %in% diseases_in_category) %>%
    distinct()
  
  ## unique md and cd based on the comorbidity filter applied
  unique_cds_temp = unique(dis_cat_specific_comorbidity$complex_disease)
  unique_mds_temp = unique(dis_cat_specific_comorbidity$mendelian_disease)
  
  ## md genes
  md_genes_temp = md_genes %>% filter(mendelian_disease %in% unique_mds_temp)
  
  ## drugs targeting md genes
  drugs_targeting_md_genes_temp = left_join(md_genes_temp, db_drug_targets, by = c("causal_gene" = "drug_target")) %>%
    na.omit() %>%
    dplyr::select(db_id) %>%
    distinct() 
  cd_md_genes_drugs = left_join(md_genes_temp, db_drug_targets, by = c("causal_gene" = "drug_target")) %>% 
    na.omit()
  
  # create logistic regression input - each row is a drug-complex disease pair 
  log_input = data.frame(db_id = rep(drugs_targeting_md_genes_temp$db_id, each = length(unique_cds_temp)))
  log_input$complex_disease = unique_cds_temp
  # add information about number of targets for each drug
  log_input = left_join(log_input, nr_drug_targets, by = "db_id")  
  # add information about investigated/indicated drugs within each drug-cancer pair
  log_input = left_join(log_input, investigated_indicated_drugs, by = c("complex_disease" = "complex_disease", 
                                                                        "db_id" = "drugbank_id"))
  # if indicated_investigated is NA --> this drug is not indicated/investigated for this cancer
  log_input$indicated_investigated = ifelse(is.na(log_input$indicated_investigated), 0, 1)
  
  # add information about comorbidity
  x = left_join(dis_cat_specific_comorbidity, cd_md_genes_drugs, by = "mendelian_disease")
  x = na.omit(x) ; rownames(x) = NULL
  x = x %>% dplyr::select(complex_disease, db_id) %>% distinct
  x$comorbidity = 1
  log_input = left_join(log_input, x, by = c("complex_disease", "db_id"))
  log_input$comorbidity = ifelse(is.na(log_input$comorbidity), 0, 1)
  
  ## run logistic regression ##
  glm_fits = glm(indicated_investigated ~ total_targets + comorbidity,
                 data = log_input, 
                 family = binomial())
  
  return(glm_fits)
}

log_reg_cancer = logistic_regression_per_dis_cat(diseases_in_category = cancers, comorbidity_dataframe = md_cd_comorbidities, nr_drug_targets = drugs_nr_targets)
log_reg_cardiovascular = logistic_regression_per_dis_cat(diseases_in_category = cardiovascular, comorbidity_dataframe = md_cd_comorbidities, nr_drug_targets = drugs_nr_targets)
log_reg_hormonal = logistic_regression_per_dis_cat(diseases_in_category = hormonal, comorbidity_dataframe = md_cd_comorbidities, nr_drug_targets = drugs_nr_targets)
log_reg_immune = logistic_regression_per_dis_cat(diseases_in_category = immune, comorbidity_dataframe = md_cd_comorbidities, nr_drug_targets = drugs_nr_targets)
log_reg_neurological = logistic_regression_per_dis_cat(diseases_in_category = neurological, comorbidity_dataframe = md_cd_comorbidities, nr_drug_targets = drugs_nr_targets)
log_reg_ophthalmological = logistic_regression_per_dis_cat(diseases_in_category = ophthalmological, comorbidity_dataframe = md_cd_comorbidities, nr_drug_targets = drugs_nr_targets)

# observed odds ratios and pvalues
cancers_obs_pvalue = summary(log_reg_cancer)$coefficients[3, 4]
cancers_obs_or = exp(summary(log_reg_cancer)$coefficients[3, 1])
cancers_obs_ci = exp(confint(log_reg_cancer)[3, ])

cardiovascular_obs_pvalue = summary(log_reg_cardiovascular)$coefficients[3, 4]
cardiovascular_obs_or = exp(summary(log_reg_cardiovascular)$coefficients[3, 1])
cardiovascular_obs_ci = exp(confint(log_reg_cardiovascular)[3, ])

hormonal_obs_pvalue = summary(log_reg_hormonal)$coefficients[3, 4]
hormonal_obs_or = exp(summary(log_reg_hormonal)$coefficients[3, 1])
hormonal_obs_ci = exp(confint(log_reg_hormonal)[3, ])

immune_obs_pvalue = summary(log_reg_immune)$coefficients[3, 4]
immune_obs_or = exp(summary(log_reg_immune)$coefficients[3, 1])
immune_obs_ci = exp(confint(log_reg_immune)[3, ])

neurological_obs_pvalue = summary(log_reg_neurological)$coefficients[3, 4]
neurological_obs_or = exp(summary(log_reg_neurological)$coefficients[3, 1])
neurological_obs_ci = exp(confint(log_reg_neurological)[3, ])

ophthalmological_obs_pvalue = summary(log_reg_ophthalmological)$coefficients[3, 4]
ophthalmological_obs_or = exp(summary(log_reg_ophthalmological)$coefficients[3, 1])
ophthalmological_obs_ci = exp(confint(log_reg_ophthalmological)[3, ])

## permutations ##

## function to run permutations per disease category
logistic_regression_per_dis_cat_perm = function(diseases_in_category, comorbidity_dataframe, nr_drug_targets) {
  
  dis_cat_specific_comorbidity = comorbidity_dataframe %>% 
    filter(complex_disease %in% diseases_in_category) %>%
    distinct()
  
  ## comorbidity matrix
  comorbidity_matrix = as.data.frame.matrix(table(dis_cat_specific_comorbidity[, 2:1]))
  comorbidity_matrix = comorbidity_matrix %>% mutate(complex_disease = rownames(comorbidity_matrix), .before = colnames(comorbidity_matrix)[1])
  rownames(comorbidity_matrix) = NULL
  
  # shuffle comorbidity matrix (comorbidities within each Mendelian disease)
  comorbidity_matrix_shuffled = comorbidity_matrix
  for (i in 2:ncol(comorbidity_matrix_shuffled)) {
    # Exclude first column - complex disease name
    comorbidity_matrix_shuffled[, i] = sample(comorbidity_matrix_shuffled[, i], replace = FALSE)
  }
  
  # transform back to original format
  comorbidity_matrix_shuffled = tidyr::gather(comorbidity_matrix_shuffled, key = "mendelian_disease", value = "comorbidity", -complex_disease) %>%
    dplyr::relocate(mendelian_disease, .before = complex_disease) %>% 
    filter(comorbidity == 1) %>% 
    dplyr::select(-comorbidity) %>%
    distinct()
  
  # unique Mendelian and complex diseases
  unique_cds_temp = unique(comorbidity_matrix_shuffled$complex_disease)
  unique_mds_temp = unique(comorbidity_matrix_shuffled$mendelian_disease)
  
  ## md genes
  md_genes_temp = md_genes %>% filter(mendelian_disease %in% unique_mds_temp)
  
  ## drugs targeting md genes
  drugs_targeting_md_genes_temp = left_join(md_genes_temp, db_drug_targets, by = c("causal_gene" = "drug_target")) %>%
    na.omit() %>%
    dplyr::select(db_id) %>%
    distinct() 
  cd_md_genes_drugs = left_join(md_genes_temp, db_drug_targets, by = c("causal_gene" = "drug_target")) %>% 
    na.omit()
  
  # create logistic regression input - each row is a drug-complex disease pair 
  log_input = data.frame(db_id = rep(drugs_targeting_md_genes_temp$db_id, each = length(unique_cds_temp)))
  log_input$complex_disease = unique_cds_temp
  # add information about number of targets for each drug
  log_input = left_join(log_input, nr_drug_targets, by = "db_id")  
  # add information about investigated/indicated drugs within each drug-cancer pair
  log_input = left_join(log_input, investigated_indicated_drugs, by = c("complex_disease" = "complex_disease", 
                                                                        "db_id" = "drugbank_id"))
  # if indicated_investigated is NA --> this drug is not indicated/investigated for this cancer
  log_input$indicated_investigated = ifelse(is.na(log_input$indicated_investigated), 0, 1)
  
  # add information about comorbidity
  x = left_join(comorbidity_matrix_shuffled, cd_md_genes_drugs, by = "mendelian_disease")
  x = na.omit(x) ; rownames(x) = NULL
  x = x %>% dplyr::select(complex_disease, db_id) %>% distinct
  x$comorbidity = 1
  log_input = left_join(log_input, x, by = c("complex_disease", "db_id"))
  log_input$comorbidity = ifelse(is.na(log_input$comorbidity), 0, 1)
  
  ## run logistic regression ##
  glm_fits = glm(indicated_investigated ~ total_targets + comorbidity,
                 data = log_input, 
                 family = binomial())
  
  return(glm_fits)
}

# cancers
cancers_or_perm = c()
for (i in 1:1000) {
  log_reg_perm = logistic_regression_per_dis_cat_perm(diseases_in_category = cancers, comorbidity_dataframe = md_cd_comorbidities, nr_drug_targets = drugs_nr_targets)
  cancers_or_perm = c(cancers_or_perm, summary(log_reg_perm)$coefficients[3,1])
  cat(i, "\n")
}
# immune
immune_or_perm = c()
for (i in 1:1000) {
  log_reg_perm = logistic_regression_per_dis_cat_perm(diseases_in_category = immune, comorbidity_dataframe = md_cd_comorbidities, nr_drug_targets = drugs_nr_targets)
  immune_or_perm = c(immune_or_perm, summary(log_reg_perm)$coefficients[3,1])
  cat(i, "\n")
}
# neurological
neurological_or_perm = c()
for (i in 1:1000) {
  log_reg_perm = logistic_regression_per_dis_cat_perm(diseases_in_category = neurological, comorbidity_dataframe = md_cd_comorbidities, nr_drug_targets = drugs_nr_targets)
  neurological_or_perm = c(neurological_or_perm, summary(log_reg_perm)$coefficients[3,1])
  cat(i, "\n")
}
# hormonal
hormonal_or_perm = c()
for (i in 1:1000) {
  log_reg_perm = logistic_regression_per_dis_cat_perm(diseases_in_category = hormonal, comorbidity_dataframe = md_cd_comorbidities, nr_drug_targets = drugs_nr_targets)
  hormonal_or_perm = c(hormonal_or_perm, summary(log_reg_perm)$coefficients[3,1])
  cat(i, "\n")
}
# ophthalmological
ophthalmological_or_perm = c()
for (i in 1:1000) {
  log_reg_perm = logistic_regression_per_dis_cat_perm(diseases_in_category = ophthalmological, comorbidity_dataframe = md_cd_comorbidities, nr_drug_targets = drugs_nr_targets)
  ophthalmological_or_perm = c(ophthalmological_or_perm, summary(log_reg_perm)$coefficients[3,1])
  cat(i, "\n")
}
# cardiovascular
cardiovascular_or_perm = c()
for (i in 1:1000) {
  log_reg_perm = logistic_regression_per_dis_cat_perm(diseases_in_category = cardiovascular, comorbidity_dataframe = md_cd_comorbidities, nr_drug_targets = drugs_nr_targets)
  cardiovascular_or_perm = c(cardiovascular_or_perm, summary(log_reg_perm)$coefficients[3,1])
  cat(i, "\n")
}

# calculate 5%, 50% and 95% quantiles for odds ratio after permutations
# cancer
cancers_perm_ci = exp(quantile(cancers_or_perm, c(0.05, 0.50, 0.95)))
# immune
immune_perm_ci = exp(quantile(immune_or_perm, c(0.05, 0.50, 0.95)))
# neuroogical
neurological_perm_ci = exp(quantile(neurological_or_perm, c(0.05, 0.50, 0.95)))
# hormonal
hormonal_perm_ci = exp(quantile(hormonal_or_perm, c(0.05, 0.50, 0.95)))
# ophthalmological
ophthalmological_perm_ci = exp(quantile(ophthalmological_or_perm, c(0.05, 0.50, 0.95)))
# cardiovascular
cardiovascular_perm_ci = exp(quantile(cardiovascular_or_perm, c(0.05, 0.50, 0.95)))

## forestplot
data_fig2a = data.frame(category = c("Neoplasms", "Neurological", "Immune", "All"),
                        odds_ratio = c(cancers_obs_or, neurological_obs_or, immune_obs_or, all_obs_or),
                        ci_2.5 = c(cancers_obs_ci["2.5 %"], neurological_obs_ci["2.5 %"], immune_obs_ci["2.5 %"], all_obs_ci["2.5 %"]),
                        ci_97.5 = c(cancers_obs_ci["97.5 %"], neurological_obs_ci["97.5 %"], immune_obs_ci["97.5 %"], all_obs_ci["97.5 %"]))
data_fig2a$category = factor(data_fig2a$category, levels = data_fig2a$category, labels = data_fig2a$category)

fig_2a = ggplot(data_fig2a, aes(y = category, x = odds_ratio, xmin = ci_2.5, xmax = ci_97.5)) +
  geom_point(size = 3) + 
  geom_errorbarh(height = 0.2, linewidth = 1) +
  ylab("") +
  xlab("Odds Ratio") +
  scale_x_continuous(breaks = seq(0, 4, 0.5)) +
  geom_vline(xintercept = 1, color = "black", linewidth = 1, linetype = 3) +
  theme_classic() +
  theme(axis.title = element_text(size = 25, family = "Arial", colour = "black"),
        axis.text = element_text(size = 25, family = "Arial", colour = "black"),
        axis.title.x = element_text(margin = margin(t = 15)),
        legend.title = element_blank())
fig_2a
ggsave(filename = "Fig2A_oddsratio_all_per_dis_cat.tiff", 
       path = "figures/", 
       width = 10, height = 6, device = "tiff",
       dpi = 300, compression = "lzw", type = type_compression)
dev.off()

## save results
odds_ratio_per_disease_category = data.frame(disease_category = c("All", "Immune", "Neurological", "Neoplasms", "Hormonal", "Ophthalmological", "Cardiovascular"),
                                             odds_ratio_observed = c(all_obs_or, immune_obs_or, neurological_obs_or, cancers_obs_or, hormonal_obs_or, ophthalmological_obs_or, cardiovascular_obs_or),
                                             pvalue = c(all_obs_pvalue, immune_obs_pvalue, neurological_obs_pvalue, cancers_obs_pvalue, hormonal_obs_pvalue, ophthalmological_obs_pvalue, cardiovascular_obs_pvalue),
                                             ci_2.5 = c(all_obs_ci["2.5 %"], immune_obs_ci["2.5 %"], neurological_obs_ci["2.5 %"], cancers_obs_ci["2.5 %"],  hormonal_obs_ci["2.5 %"], ophthalmological_obs_ci["2.5 %"], cardiovascular_obs_ci["2.5 %"]),
                                             ci_97.5 = c(all_obs_ci["97.5 %"], immune_obs_ci["97.5 %"], neurological_obs_ci["97.5 %"], cancers_obs_ci["97.5 %"],  hormonal_obs_ci["97.5 %"], ophthalmological_obs_ci["97.5 %"], cardiovascular_obs_ci["97.5 %"]),
                                             ci_5_permutations = c(all_perm_ci["5%"], immune_perm_ci["5%"], neurological_perm_ci["5%"], cancers_perm_ci["5%"],  hormonal_perm_ci["5%"], ophthalmological_perm_ci["5%"], cardiovascular_perm_ci["5%"]),
                                             ci_50_permutations = c(all_perm_ci["50%"], immune_perm_ci["50%"], neurological_perm_ci["50%"], cancers_perm_ci["50%"],  hormonal_perm_ci["50%"], ophthalmological_perm_ci["50%"], cardiovascular_perm_ci["50%"]),
                                             ci_95_permutations = c(all_perm_ci["95%"], immune_perm_ci["95%"], neurological_perm_ci["95%"], cancers_perm_ci["95%"],  hormonal_perm_ci["95%"], ophthalmological_perm_ci["95%"], cardiovascular_perm_ci["95%"]),
                                             p_perm = c(sum(all_obs_or <= exp(all_or_perm)),
                                                        sum(immune_obs_or <= exp(immune_or_perm)),
                                                        sum(neurological_obs_or <= exp(neurological_or_perm)),
                                                        sum(cancers_obs_or <= exp(cancers_or_perm)),
                                                        sum(hormonal_obs_or <= exp(hormonal_or_perm)),
                                                        sum(ophthalmological_obs_or <= exp(ophthalmological_or_perm)),
                                                        sum(cardiovascular_obs_or <= exp(cardiovascular_or_perm))))
dir.create("results")
fwrite(odds_ratio_per_disease_category, "results/log_reg_results_per_disease_category.txt", sep = "\t", row.names = FALSE)

## Figure 2B: Example of candidate drugs from two Mendelian diseases and cancers ##
# indicated/investigated drugs
investigated_indicated_drugs_cancers =  fread("processed_data/drugs_inv_ind_per_disease.txt")
investigated_indicated_drugs_cancers = reshape2::melt(investigated_indicated_drugs_cancers, "drugbank_id", colnames(investigated_indicated_drugs_cancers)[2:ncol(investigated_indicated_drugs_cancers)])
investigated_indicated_drugs_cancers = na.omit(investigated_indicated_drugs_cancers)
rownames(investigated_indicated_drugs_cancers) = NULL
# filter for drugs ind/inv for at least one cancer
investigated_indicated_drugs_cancers = investigated_indicated_drugs_cancers %>%
  dplyr::select(drugbank_id, complex_disease = variable, phase = value) %>%
  distinct() %>%
  filter(complex_disease %in% cancers)
# cancer comorbidities
md_cancer_comorbidities = md_cd_comorbidities %>% 
  dplyr::select(mendelian_disease, complex_disease) %>%
  filter(complex_disease %in% cancers) %>%
  arrange(mendelian_disease)
unique_md_cancers = unique(md_cancer_comorbidities$mendelian_disease)
# approved drug - targets
db_approved_drugs = fread("raw_data/drug_links_approved.csv") %>% dplyr::select(db_id = "DrugBank ID", drug_name = "Name") %>% mutate(approved = 1) %>% distinct()
db_drug_targets_approve_status = left_join(db_drug_targets, db_approved_drugs, by = "db_id") ; rm(db_approved_drugs)
db_drug_targets_approve_status$approved = ifelse(is.na(db_drug_targets_approve_status$approved), 0, 1)
approved_drugs = db_drug_targets_approve_status %>% filter(approved == 1) %>%
  dplyr::select(db_id, drug_name) %>%
  distinct()

# approved drugs targeting AR
ar_genes = md_genes %>% filter(mendelian_disease == "Androgen Insensitivity Syndrome")
ar_genes = left_join(ar_genes, db_drug_targets, by = c("causal_gene" = "drug_target"))
ar_genes = left_join(ar_genes, approved_drugs, by = "db_id")
ar_genes = na.omit(ar_genes) ; rownames(ar_genes) = NULL
ar_genes = ar_genes %>% filter(db_id %in% investigated_indicated_drugs_cancers$drugbank_id)
ar_drugs = unique(ar_genes$db_id)
ar_drugs = ar_drugs[-c(12, 17, 27,28,29)]

# approved drugs targeting Retinitis Pigmentosa
ret_pigm_genes = md_genes %>% filter(mendelian_disease == "Retinitis Pigmentosa")
ret_pigm_genes = left_join(ret_pigm_genes, db_drug_targets, by = c("causal_gene" = "drug_target"))
ret_pigm_genes = left_join(ret_pigm_genes, approved_drugs, by = "db_id")
ret_pigm_genes = na.omit(ret_pigm_genes) ; rownames(ret_pigm_genes) = NULL
ret_pigm_genes = ret_pigm_genes %>% filter(db_id %in% investigated_indicated_drugs_cancers$drugbank_id)
ret_pig_drugs = unique(ret_pigm_genes$db_id)
ret_pig_drugs = ret_pig_drugs[-2]

drugs = unique(c(ar_drugs, ret_pig_drugs))
cancers = c("Female Breast Cancer", 
            "Hyperplasia of the Prostate",
            "Prostate Cancer",
            "Uterine Cancer",
            "Benign Brain Neoplasm", 
            "Malignant Brain Neoplasm",
            "Lung Cancer",
            "Melanoma",
            "Colorectal Cancer",
            "Bladder Cancer",
            "Kidney Cancer",
            "Burkitts Lymphoma", 
            "Lymphosarcoma/Reticulorsarcoma",
            "Gastric Cancer")

heatmap_matrix = matrix(ncol = length(cancers), nrow = length(drugs))
colnames(heatmap_matrix) = cancers
rownames(heatmap_matrix) = drugs

# clinical trial phase for each drug - cancer pair
for (row in 1:nrow(heatmap_matrix)) {
  for (col in 1:ncol(heatmap_matrix)) {
    
    drug = rownames(heatmap_matrix)[row]
    cd = colnames(heatmap_matrix)[col]
    
    drugs_cd = investigated_indicated_drugs_cancers %>% filter(complex_disease == cd)
    
    if (drug %in% drugs_cd$drugbank_id) {
      phase = drugs_cd %>% filter(drugbank_id == drug)
      if (phase$phase >= 4) {
        heatmap_matrix[row, col] = "Approved"
      }
      if (phase$phase == 3 | phase$phase == 3.5 | phase$phase == 2.5) {
        heatmap_matrix[row, col] = "Phase 3"
      }
      if (phase$phase == 2 | phase$phase == 1.5) {
        heatmap_matrix[row, col] = "Phase 2"
      }
      if (phase$phase == 1 | phase$phase == 0.5) {
        heatmap_matrix[row, col] = "Phase 1"
      }
      if (phase$phase == 0) {
        heatmap_matrix[row, col] = "Phase 1"
      }
    }
    if (!drug %in% drugs_cd$drugbank_id) {
      heatmap_matrix[row, col] = "Not indicated/investigated"
    }
  }
}

# information about comorbidity
comorbidities = matrix(ncol = length(cancers), nrow = length(drugs))
colnames(comorbidities) = cancers
rownames(comorbidities) = drugs

# populate comorbidity matrix
# Androgen Insensitivity Syndrome
for (row in 1:24) {
  for (col in 1:14) {
    
    cd = colnames(comorbidities)[col]
    cd_comorbidities = md_cd_comorbidities %>% filter(complex_disease == cd)
    
    if (("Androgen Insensitivity Syndrome") %in% cd_comorbidities$mendelian_disease) {
      comorbidities[row, col] = 1
    }
    if (!("Androgen Insensitivity Syndrome") %in% cd_comorbidities$mendelian_disease) {
      comorbidities[row, col] = 0
    }
  }
}
# Retinitis pigmentosa
for (row in 25:34) {
  for (col in 1:14) {
    
    cd = colnames(comorbidities)[col]
    cd_comorbidities = md_cd_comorbidities %>% filter(complex_disease == cd)
    
    if (("Retinitis Pigmentosa") %in% cd_comorbidities$mendelian_disease) {
      comorbidities[row, col] = 1
    }
    if (!("Retinitis Pigmentosa") %in% cd_comorbidities$mendelian_disease) {
      comorbidities[row, col] = 0
    }
  }
}

# rename columns and rows
colnames(heatmap_matrix) = c("Breast cancer", 
                             "Hyperplasia of prostate", 
                             "Prostate cancer",
                             "Uterine cancer", 
                             "Benign Brain neoplasm",
                             "Malignant Brain neoplasm", 
                             "Lung cancer", 
                             "Melanoma",
                             "Colorectal cancer",
                             "Bladder cancer", 
                             "Kidney cancer", 
                             "Burkitt's lymphoma",
                             "Lymphoma", 
                             "Gastric cancer")
rownames_drugs = data.frame(db_id = rownames(heatmap_matrix))
rownames_drugs = left_join(rownames_drugs, approved_drugs, by = "db_id")
rownames_drugs[17, "drug_name"] = "Cyproterone"
rownames_drugs[19, "drug_name"] = "Nandrolone"
rownames(heatmap_matrix) = rownames_drugs$drug_name

# define colors
colors = structure(c("white", "#E5E4E2", "#D3D3D3", "#A9A9A9", "#808080"), names = c("Not indicated/investigated", "Phase 1", "Phase 2", "Phase 3", "Approved"))

ht = Heatmap(heatmap_matrix,
             name = "\n\t",
             row_names_side = "left",
             column_names_side = "bottom",
             col = colors,
             row_dend_reorder = FALSE,
             cluster_rows = FALSE,
             column_dend_reorder = FALSE,
             show_column_dend = FALSE,
             cluster_columns = FALSE,
             show_heatmap_legend = FALSE,
             rect_gp = gpar(col = "black", lwd = 0.5),
             width = unit(10, "cm"),
             height = unit(15, "cm"), column_names_rot = 45,
             cell_fun = function(j, i, x, y, w, h, fill) {
               if (comorbidities[i, j] == 1) {
                 grid.points(x, y, pch = 16, size = unit(1.4, "mm"), gp = gpar(col = "black"))
                 } else if (comorbidities[i, j] != 1) {
                 grid.text("", x, y)
                 }
               },
             left_annotation = rowAnnotation("\n" = c(rep("Androgen Insensitivity Syndrome", 24), 
                                                      rep("Retinitis Pigmentosa", 10)),
                                             border = FALSE))
lgd_list = list(
  Legend(labels = "Recommended", title = "", 
         type = "points", 
         pch = 16, 
         background = "white", 
         border = TRUE, 
         legend_gp = gpar(col = 1:2)),
  Legend(labels = c("Not indicated/investigated", "Phase 1", "Phase 2", "Phase 3", "Approved"),
         title = "",
         background = "white", 
         border = TRUE, 
         legend_gp = gpar(fill = c("white", "#E5E4E2", "#D3D3D3", "#A9A9A9", "#808080")))
)

tiff("figures/Fig2B_example_cancers.tiff",
     width = 30, height = 20, units = "cm", 
     res = 1200, compression = "lzw", type = type_compression)
draw(ht, annotation_legend_list = lgd_list, annotation_legend_side = "right")
dev.off()

#### --- investigated OR approved drugs only --- ####

## indicated/investigated drugs
investigated_indicated_drugs = fread("processed_data/drugs_inv_ind_per_disease.txt")
investigated_indicated_drugs = reshape2::melt(investigated_indicated_drugs, "drugbank_id", colnames(investigated_indicated_drugs)[2:ncol(investigated_indicated_drugs)])
investigated_indicated_drugs = na.omit(investigated_indicated_drugs)
rownames(investigated_indicated_drugs) = NULL
investigated_indicated_drugs = investigated_indicated_drugs %>% dplyr::select(drugbank_id, complex_disease = variable, phase = value) %>% distinct()
investigated_indicated_drugs$phase = ifelse(investigated_indicated_drugs$phase == 1.5, 2, investigated_indicated_drugs$phase)
investigated_indicated_drugs$phase = ifelse(investigated_indicated_drugs$phase == 2.5, 3, investigated_indicated_drugs$phase)
investigated_indicated_drugs$phase = ifelse(investigated_indicated_drugs$phase == 5, 4, investigated_indicated_drugs$phase)
# phases 1, 2 or 3 (exclude unknown phase)
investigated_indicated_drugs_phase123 = investigated_indicated_drugs %>% filter(phase %in% c(1, 2, 3))
investigated_indicated_drugs_phase123$indicated_investigated = 1
# phase 1
investigated_indicated_drugs_phase1 = investigated_indicated_drugs %>% filter(phase == 1)
investigated_indicated_drugs_phase1$indicated_investigated = 1
# phase 2
investigated_indicated_drugs_phase2 = investigated_indicated_drugs %>% filter(phase == 2)
investigated_indicated_drugs_phase2$indicated_investigated = 1
# phase 3
investigated_indicated_drugs_phase3 = investigated_indicated_drugs %>% filter(phase == 3)
investigated_indicated_drugs_phase3$indicated_investigated = 1
# approved
investigated_indicated_drugs_approved = investigated_indicated_drugs %>% filter(phase == 4)
investigated_indicated_drugs_approved$indicated_investigated = 1

## logistic regression - observed ##
# create logistic regression input - each row is a drug-complex disease pair 
log_input_ct_all_phases = data.frame(db_id = rep(drugs_targeting_md_genes$db_id, each = length(unique_cd)))
log_input_ct_all_phases$complex_disease = unique_cd
# add information about disease category
log_input_ct_all_phases = left_join(log_input_ct_all_phases, complex_disease_categories, by = "complex_disease")
# add information about number of targets for each drug
log_input_ct_all_phases = left_join(log_input_ct_all_phases, drugs_nr_targets, by = "db_id")  
# add information about comorbidity
x = left_join(md_cd_comorbidities, md_genes, by = "mendelian_disease")
x = na.omit(x) ; rownames(x) = NULL
x = x %>% dplyr::select(complex_disease, causal_gene) %>% distinct
x = left_join(x, db_drug_targets, by = c("causal_gene" = "drug_target"))
x = na.omit(x) ; rownames(x) = NULL
x = x %>% dplyr::select(-causal_gene) %>% distinct()
x$comorbidity = 1
log_input_ct_all_phases = left_join(log_input_ct_all_phases, x, by = c("complex_disease", "db_id"))
log_input_ct_all_phases$comorbidity = ifelse(is.na(log_input_ct_all_phases$comorbidity), 0, 1)

# add information about investigated/indicated drugs within each drug-cancer pair
log_input_ct_all_phases = left_join(log_input_ct_all_phases, investigated_indicated_drugs_phase123[, -3], by = c("complex_disease" = "complex_disease", "db_id" = "drugbank_id")) %>%
  dplyr::rename(indicated_investigated_phase123 = indicated_investigated)
log_input_ct_all_phases$indicated_investigated_phase123 = ifelse(is.na(log_input_ct_all_phases$indicated_investigated_phase123), 0, 1)

log_input_ct_all_phases = left_join(log_input_ct_all_phases, investigated_indicated_drugs_phase1[, -3], by = c("complex_disease" = "complex_disease", "db_id" = "drugbank_id")) %>%
  dplyr::rename(indicated_investigated_phase1 = indicated_investigated)
log_input_ct_all_phases$indicated_investigated_phase1 = ifelse(is.na(log_input_ct_all_phases$indicated_investigated_phase1), 0, 1)

log_input_ct_all_phases = left_join(log_input_ct_all_phases, investigated_indicated_drugs_phase2[, -3], by = c("complex_disease" = "complex_disease", "db_id" = "drugbank_id")) %>%
  dplyr::rename(indicated_investigated_phase2 = indicated_investigated)
log_input_ct_all_phases$indicated_investigated_phase2 = ifelse(is.na(log_input_ct_all_phases$indicated_investigated_phase2), 0, 1)

log_input_ct_all_phases = left_join(log_input_ct_all_phases, investigated_indicated_drugs_phase3[, -3], by = c("complex_disease" = "complex_disease", "db_id" = "drugbank_id")) %>%
  dplyr::rename(indicated_investigated_phase3 = indicated_investigated)
log_input_ct_all_phases$indicated_investigated_phase3 = ifelse(is.na(log_input_ct_all_phases$indicated_investigated_phase3), 0, 1)

log_input_ct_all_phases = left_join(log_input_ct_all_phases, investigated_indicated_drugs_approved[, -3], by = c("complex_disease" = "complex_disease", "db_id" = "drugbank_id")) %>%
  dplyr::rename(indicated_investigated_approved = indicated_investigated)
log_input_ct_all_phases$indicated_investigated_approved = ifelse(is.na(log_input_ct_all_phases$indicated_investigated_approved), 0, 1)

# run logistic regression
# phase 123
log_reg_ct_phase123 = glm(indicated_investigated_phase123 ~ total_targets + disease_category + comorbidity,
                        data = log_input_ct_all_phases, 
                        family = binomial())
ct_phase123_obs_or = exp(summary(log_reg_ct_phase123)$coefficients[8, 1])
ct_phase123_obs_pvalue = summary(log_reg_ct_phase123)$coefficients[8, 4]
ct_phase123_obs_ci = exp(confint(log_reg_ct_phase123))[8, ]
# phase 1
log_reg_ct_phase1 = glm(indicated_investigated_phase1 ~ total_targets + disease_category + comorbidity,
                          data = log_input_ct_all_phases, 
                          family = binomial())
ct_phase1_obs_or = exp(summary(log_reg_ct_phase1)$coefficients[8, 1])
ct_phase1_obs_pvalue = summary(log_reg_ct_phase1)$coefficients[8, 4]
ct_phase1_obs_ci = exp(confint(log_reg_ct_phase1))[8, ]
# phase 2
log_reg_ct_phase2 = glm(indicated_investigated_phase2 ~ total_targets + disease_category + comorbidity,
                        data = log_input_ct_all_phases, 
                        family = binomial())
ct_phase2_obs_or = exp(summary(log_reg_ct_phase2)$coefficients[8, 1])
ct_phase2_obs_pvalue = summary(log_reg_ct_phase2)$coefficients[8, 4]
ct_phase2_obs_ci = exp(confint(log_reg_ct_phase2))[8, ]
# phase 3
log_reg_ct_phase3 = glm(indicated_investigated_phase3 ~ total_targets + disease_category + comorbidity,
                        data = log_input_ct_all_phases, 
                        family = binomial())
ct_phase3_obs_or = exp(summary(log_reg_ct_phase3)$coefficients[8, 1])
ct_phase3_obs_pvalue = summary(log_reg_ct_phase3)$coefficients[8, 4]
ct_phase3_obs_ci = exp(confint(log_reg_ct_phase3))[8, ]
# approved
log_reg_ct_approved = glm(indicated_investigated_approved ~ total_targets + disease_category + comorbidity,
                          data = log_input_ct_all_phases, 
                          family = binomial())
ct_approved_obs_or = exp(summary(log_reg_ct_approved)$coefficients[8, 1])
ct_approved_obs_pvalue = summary(log_reg_ct_approved)$coefficients[8, 4]
ct_approved_obs_ci = exp(confint(log_reg_ct_approved))[8, ]

## permutations ##

ct_phase123_perm_or = c()
ct_phase1_perm_or = c()
ct_phase2_perm_or = c()
ct_phase3_perm_or = c()
ct_approved_perm_or = c()

ct_phase123_perm_pvalue = c()
ct_phase1_perm_pvalue = c()
ct_phase2_perm_pvalue = c()
ct_phase3_perm_pvalue = c()
ct_approved_perm_pvalue = c()

for (permutation in 1:1000) {
  
  ## shuffle comorbidity matrix
  md_cd_comorbidities_shuffled = shuffle(md_cd_comorbidities_matrix)
  
  # create logistic regression input - each row is a drug-complex disease pair 
  log_input_perm_ct_all_phases = data.frame(db_id = rep(drugs_targeting_md_genes$db_id, each = length(unique_cd)))
  log_input_perm_ct_all_phases$complex_disease = unique_cd
  # add information about disease category
  log_input_perm_ct_all_phases = left_join(log_input_perm_ct_all_phases, complex_disease_categories, by = "complex_disease")
  # add information about number of targets for each drug
  log_input_perm_ct_all_phases = left_join(log_input_perm_ct_all_phases, drugs_nr_targets, by = "db_id")  
  # add information about comorbidity
  x = left_join(md_cd_comorbidities_shuffled, md_genes, by = "mendelian_disease")
  x = na.omit(x) ; rownames(x) = NULL
  x = x %>% dplyr::select(complex_disease, causal_gene) %>% distinct
  x = left_join(x, db_drug_targets, by = c("causal_gene" = "drug_target"))
  x = na.omit(x) ; rownames(x) = NULL
  x = x %>% dplyr::select(-causal_gene) %>% distinct()
  x$comorbidity = 1
  log_input_perm_ct_all_phases = left_join(log_input_perm_ct_all_phases, x, by = c("complex_disease", "db_id"))
  log_input_perm_ct_all_phases$comorbidity = ifelse(is.na(log_input_perm_ct_all_phases$comorbidity), 0, 1)
  
  # add information about investigated/indicated drugs within each drug-cancer pair
  log_input_perm_ct_all_phases = left_join(log_input_perm_ct_all_phases, investigated_indicated_drugs_phase123[, -3], by = c("complex_disease" = "complex_disease", "db_id" = "drugbank_id")) %>%
    dplyr::rename(indicated_investigated_phase123 = indicated_investigated)
  log_input_perm_ct_all_phases$indicated_investigated_phase123 = ifelse(is.na(log_input_perm_ct_all_phases$indicated_investigated_phase123), 0, 1)
  
  log_input_perm_ct_all_phases = left_join(log_input_perm_ct_all_phases, investigated_indicated_drugs_phase1[, -3], by = c("complex_disease" = "complex_disease", "db_id" = "drugbank_id")) %>%
    dplyr::rename(indicated_investigated_phase1 = indicated_investigated)
  log_input_perm_ct_all_phases$indicated_investigated_phase1 = ifelse(is.na(log_input_perm_ct_all_phases$indicated_investigated_phase1), 0, 1)
  
  log_input_perm_ct_all_phases = left_join(log_input_perm_ct_all_phases, investigated_indicated_drugs_phase2[, -3], by = c("complex_disease" = "complex_disease", "db_id" = "drugbank_id")) %>%
    dplyr::rename(indicated_investigated_phase2 = indicated_investigated)
  log_input_perm_ct_all_phases$indicated_investigated_phase2 = ifelse(is.na(log_input_perm_ct_all_phases$indicated_investigated_phase2), 0, 1)
  
  log_input_perm_ct_all_phases = left_join(log_input_perm_ct_all_phases, investigated_indicated_drugs_phase3[, -3], by = c("complex_disease" = "complex_disease", "db_id" = "drugbank_id")) %>%
    dplyr::rename(indicated_investigated_phase3 = indicated_investigated)
  log_input_perm_ct_all_phases$indicated_investigated_phase3 = ifelse(is.na(log_input_perm_ct_all_phases$indicated_investigated_phase3), 0, 1)
  
  log_input_perm_ct_all_phases = left_join(log_input_perm_ct_all_phases, investigated_indicated_drugs_approved[, -3], by = c("complex_disease" = "complex_disease", "db_id" = "drugbank_id")) %>%
    dplyr::rename(indicated_investigated_approved = indicated_investigated)
  log_input_perm_ct_all_phases$indicated_investigated_approved = ifelse(is.na(log_input_perm_ct_all_phases$indicated_investigated_approved), 0, 1)

  # run logistic regression
  # phase 123
  log_reg_perm_ct_phase123 = glm(indicated_investigated_phase123 ~ total_targets + disease_category + comorbidity,
                            data = log_input_perm_ct_all_phases, 
                            family = binomial())
  ct_phase123_perm_or = c(ct_phase123_perm_or, exp(summary(log_reg_perm_ct_phase123)$coefficients[8, 1]))
  ct_phase123_perm_pvalue = c(ct_phase123_perm_pvalue, summary(log_reg_perm_ct_phase123)$coefficients[8, 4])

  # phase 1
  log_reg_perm_ct_phase1 = glm(indicated_investigated_phase1 ~ total_targets + disease_category + comorbidity,
                          data = log_input_perm_ct_all_phases, 
                          family = binomial())
  ct_phase1_perm_or = c(ct_phase1_perm_or, exp(summary(log_reg_perm_ct_phase1)$coefficients[8, 1]))
  ct_phase1_perm_pvalue = c(ct_phase1_perm_pvalue, summary(log_reg_perm_ct_phase1)$coefficients[8, 4])
  
  # phase 2
  log_reg_perm_ct_phase2 = glm(indicated_investigated_phase2 ~ total_targets + disease_category + comorbidity,
                          data = log_input_perm_ct_all_phases, 
                          family = binomial())
  ct_phase2_perm_or = c(ct_phase2_perm_or, exp(summary(log_reg_perm_ct_phase2)$coefficients[8, 1]))
  ct_phase2_perm_pvalue = c(ct_phase2_perm_pvalue, summary(log_reg_perm_ct_phase2)$coefficients[8, 4])
  
  # phase 3
  log_reg_perm_ct_phase3 = glm(indicated_investigated_phase3 ~ total_targets + disease_category + comorbidity,
                          data = log_input_perm_ct_all_phases, 
                          family = binomial())
  ct_phase3_perm_or = c(ct_phase3_perm_or, exp(summary(log_reg_perm_ct_phase3)$coefficients[8, 1]))
  ct_phase3_perm_pvalue = c(ct_phase3_perm_pvalue, summary(log_reg_perm_ct_phase3)$coefficients[8, 4])
  
  # approved
  log_reg_perm_ct_approved = glm(indicated_investigated_approved ~ total_targets + disease_category + comorbidity,
                            data = log_input_perm_ct_all_phases, 
                            family = binomial())
  ct_approved_perm_or = c(ct_approved_perm_or, exp(summary(log_reg_perm_ct_approved)$coefficients[8, 1]))
  ct_approved_perm_pvalue = c(ct_approved_perm_pvalue, summary(log_reg_perm_ct_approved)$coefficients[8, 4])
  
  # track progress
  cat(permutation, "\n")
}

# calculate permutation p-value for phase123
sum(ct_phase123_obs_pvalue >= ct_phase123_perm_pvalue) / 1000 # 0
sum(ct_phase123_obs_or <= ct_phase123_perm_or) / 1000 # 0

# calculate permutation 5%, 50%, 95% quantiles
ct_phase123_perm_ci = quantile(ct_phase123_perm_or, c(0.05, 0.50, 0.95))
ct_phase1_perm_ci = quantile(ct_phase1_perm_or, c(0.05, 0.50, 0.95))
ct_phase2_perm_ci = quantile(ct_phase2_perm_or, c(0.05, 0.50, 0.95))
ct_phase3_perm_ci = quantile(ct_phase3_perm_or, c(0.05, 0.50, 0.95))
ct_approved_perm_ci = quantile(ct_approved_perm_or, c(0.05, 0.50, 0.95))

### forest plot ###
data_fig2c = data.frame(category = c("Approved", "Phase III", "Phase II", "Phase I"),
                        odds_ratio = c(ct_approved_obs_or, ct_phase3_obs_or, ct_phase2_obs_or, ct_phase1_obs_or),
                        ci_2.5 = c(ct_approved_obs_ci["2.5 %"], ct_phase3_obs_ci["2.5 %"], ct_phase2_obs_ci["2.5 %"], ct_phase1_obs_ci["2.5 %"]),
                        ci_97.5 = c(ct_approved_obs_ci["97.5 %"], ct_phase3_obs_ci["97.5 %"], ct_phase2_obs_ci["97.5 %"], ct_phase1_obs_ci["97.5 %"]))
data_fig2c$category = factor(data_fig2c$category, levels = data_fig2c$category, labels = data_fig2c$category)

fig_2c = ggplot(data_fig2c, aes(y = category, x = odds_ratio, xmin = ci_2.5, xmax = ci_97.5)) +
  geom_point(size = 3) + 
  geom_errorbarh(height = 0.2, linewidth = 1) +
  ylab("") +
  xlab("Odds Ratio") +
  scale_x_continuous(breaks = seq(0, 3, 0.5)) +
  geom_vline(xintercept = 1, color = "black", linewidth = 1, linetype = 3) +
  theme_classic() +
  theme(axis.title = element_text(size = 25, family = "Arial", colour = "black"),
        axis.text = element_text(size = 25, family = "Arial", colour = "black"),
        axis.title.x = element_text(margin = margin(t = 15)),
        legend.title = element_blank(),
        legend.text = element_text(size = 24, family = "Arial", colour = "black"),
        legend.key.size = unit(1, 'cm'))

fig_2c
ggsave(filename = "Fig2C_oddsratio_ct_phases.tiff", 
       path = "figures/", 
       width = 8, height = 4, device = "tiff",
       dpi = 300, compression = "lzw", type = type_compression)
dev.off()

## save results
odds_ratio_per_ct_phase = data.frame(disease_category = c("Phase I/II/III", "Phase I", "Phase II", "Phase III", "Approved"),
                                     odds_ratio_observed = c(ct_phase123_obs_or, ct_phase1_obs_or, ct_phase2_obs_or, ct_phase3_obs_or, ct_approved_obs_or),
                                     pvalue = c(ct_phase123_obs_pvalue, ct_phase1_obs_pvalue, ct_phase2_obs_pvalue, ct_phase3_obs_pvalue, ct_approved_obs_pvalue),
                                     ci_2.5 = c(ct_phase123_obs_ci["2.5 %"], ct_phase1_obs_ci["2.5 %"], ct_phase2_obs_ci["2.5 %"], ct_phase3_obs_ci["2.5 %"],  ct_approved_obs_ci["2.5 %"]),
                                     ci_97.5 = c(ct_phase123_obs_ci["97.5 %"], ct_phase1_obs_ci["97.5 %"], ct_phase2_obs_ci["97.5 %"], ct_phase3_obs_ci["97.5 %"],  ct_approved_obs_ci["97.5 %"]),
                                     ci_5_permutations = c(ct_phase123_perm_ci["5%"], ct_phase1_perm_ci["5%"], ct_phase2_perm_ci["5%"], ct_phase3_perm_ci["5%"],  ct_approved_perm_ci["5%"]),
                                     ci_50_permutations = c(ct_phase123_perm_ci["50%"], ct_phase1_perm_ci["50%"], ct_phase2_perm_ci["50%"], ct_phase3_perm_ci["50%"],  ct_approved_perm_ci["50%"]),
                                     ci_95_permutations = c(ct_phase123_perm_ci["95%"], ct_phase1_perm_ci["95%"], ct_phase2_perm_ci["95%"], ct_phase3_perm_ci["95%"],  ct_approved_perm_ci["95%"]))
fwrite(odds_ratio_per_ct_phase, "results/log_reg_results_per_clinical_trial_phase.txt", sep = "\t", row.names = FALSE)
