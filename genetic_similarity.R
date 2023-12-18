
library(dplyr)
library(data.table)
library(ggplot2)
library(openxlsx)

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

### -- comorbidity VS comorbidity & genetic similarity -- ###
## using genetic similarity data from: https://www.nature.com/articles/ncomms8033 (Melamed et al., 2015)
# load genetic similarity data
md_cancer_gensim_melamed = read.xlsx("raw_data/melamed_genetic_similarity_sup_data_4.xlsx", startRow = 9)

# keep columns containing p-values for each genetic similarity metric
md_cancer_gensim_melamed = md_cancer_gensim_melamed %>% 
  dplyr::select(mendelian_disease = MD,
                tcga_acronym = C,
                gene_enrichment, pathway_correlation, coex_CG, humannet_set, biogrid_set) %>%
  distinct()

# find mendelian-cancer pairs with at least one significant genetic similarity metric (p<0.05)
md_cancer_gensim_melamed = md_cancer_gensim_melamed %>%
  rowwise() %>%
  mutate(min_pvalue = min(gene_enrichment, pathway_correlation, coex_CG, humannet_set, biogrid_set)) %>%
  filter(min_pvalue < 0.05) %>%
  dplyr::select(mendelian_disease, tcga_acronym) %>%
  distinct()

# convert cancer abbreviations to full names
tcga_cancer_acronyms = fread("raw_data/cancers_tcga_acronyms.txt")
md_cancer_gensim_melamed = left_join(md_cancer_gensim_melamed, tcga_cancer_acronyms, by = "tcga_acronym") ; rm(tcga_cancer_acronyms)
md_cancer_gensim_melamed = md_cancer_gensim_melamed %>% 
  dplyr::select(-tcga_acronym) %>%
  mutate(genetic_similarity_melamed = 1) %>%
  distinct()

## MD-cancer comorbidities 
md_cancer_comorbidities = fread("processed_data/md_cd_comorbidities.txt") %>%
  # filter for cancers with available information for genetic similarity
  filter(complex_disease %in% md_cancer_gensim_melamed$complex_disease)
# add genetic similarity information
md_cancer_comorbidities = left_join(md_cancer_comorbidities, md_cancer_gensim_melamed, by = c("mendelian_disease", "complex_disease"))
md_cancer_comorbidities$genetic_similarity_melamed = ifelse(is.na(md_cancer_comorbidities$genetic_similarity_melamed), 0, 1)

## mendelian disease genes and drugs
md_genes = fread("processed_data/md_genes.txt")
db_drug_targets = fread("processed_data/drugbank_all_drugs_known_targets.txt") %>% 
  dplyr::select(db_id, drug_target) %>%
  distinct()
drugs_nr_targets = fread("processed_data/drugbank_all_drugs_known_targets.txt") %>%
  dplyr::select(db_id, drug_target) %>%
  group_by(db_id) %>%
  mutate(total_targets = length(unique(drug_target))) %>%
  ungroup() %>%
  dplyr::select(db_id, total_targets) %>%
  distinct()
  
md_drugs = md_genes %>%
  left_join(db_drug_targets, by = c("causal_gene" = "drug_target")) %>%
  na.omit() %>%
  dplyr::select(-causal_gene) %>%
  filter(mendelian_disease %in% md_cancer_comorbidities$mendelian_disease) %>%
  distinct()

## investigated/indicated drugs for the complex diseases
investigated_indicated_drugs = data.table::fread("processed_data/drugs_inv_ind_per_disease.txt")
investigated_indicated_drugs = reshape2::melt(investigated_indicated_drugs, "drugbank_id", colnames(investigated_indicated_drugs)[2:ncol(investigated_indicated_drugs)])
investigated_indicated_drugs = na.omit(investigated_indicated_drugs)
rownames(investigated_indicated_drugs) = NULL
investigated_indicated_drugs = investigated_indicated_drugs %>% dplyr::select(drugbank_id, complex_disease = variable) %>% distinct()
investigated_indicated_drugs$indicated_investigated = 1

### -- logistic regression -- ###

# drugs targeting MD genes
drugs_targeting_md_genes = md_drugs %>%
  dplyr::select(db_id) %>%
  distinct()

# unique cancers in the data
unique_cancers = unique(md_cancer_comorbidities$complex_disease)  

# create logistic regression input
log_input = data.frame(db_id = rep(drugs_targeting_md_genes$db_id, each = length(unique_cancers)))

log_input$complex_disease = unique_cancers

log_input = left_join(log_input, drugs_nr_targets, by = "db_id")  

log_input = left_join(log_input, investigated_indicated_drugs, by = c("complex_disease", "db_id" = "drugbank_id"))
log_input$indicated_investigated = ifelse(is.na(log_input$indicated_investigated), 0, 1)

cancer_recommended_drugs = left_join(md_cancer_comorbidities, md_drugs, by = "mendelian_disease")
cancer_recommended_drugs = cancer_recommended_drugs %>%
  dplyr::select(db_id, complex_disease, comorbidity, genetic_similarity_melamed) %>%
  distinct() %>%
  group_by(complex_disease, db_id) %>%
  mutate(comorbidity = max(comorbidity),
         genetic_similarity_melamed = max(genetic_similarity_melamed)) %>%
  ungroup() %>%
  distinct()

log_input = left_join(log_input, cancer_recommended_drugs, by = c("complex_disease", "db_id"))
log_input$comorbidity = ifelse(is.na(log_input$comorbidity), 0, log_input$comorbidity)
log_input$genetic_similarity_melamed = ifelse(is.na(log_input$genetic_similarity_melamed), 0, log_input$genetic_similarity_melamed)

log_input$comorbidity_gensim = log_input$comorbidity * log_input$genetic_similarity_melamed

# comorbidity (regardless of genetic similarity)
glm_fits_comorbidity = glm(indicated_investigated ~ total_targets + comorbidity,
                           data = log_input, 
                           family = binomial())
log_summary_comorbidity = summary(glm_fits_comorbidity)$coefficients
log_summary_comorbidity
# confidence intervals
or_comorbidity = exp(log_summary_comorbidity[3, 1])
pvalue_comorbidity = log_summary_comorbidity[3, 4]
or_comorbidity_95ci = c(exp(confint(glm_fits_comorbidity))[3, 1], exp(confint(glm_fits_comorbidity))[3, 2])
names(or_comorbidity_95ci) = c("ci_2.5", "ci_97.5")

# comorbidity AND genetic similarity
glm_fits_comorbidity_gensim = glm(indicated_investigated ~ total_targets + comorbidity_gensim,
                                  data = log_input, 
                                  family = binomial())
log_summary_comorbidity_gensim = summary(glm_fits_comorbidity_gensim)$coefficients
log_summary_comorbidity_gensim
# confidence intervals
or_comorbidity_gensim = exp(log_summary_comorbidity_gensim[3, 1])
pvalue_comorbidity_gensim = log_summary_comorbidity_gensim[3, 4]
or_comorbidity_gensim_95ci = c(exp(confint(glm_fits_comorbidity_gensim))[3, 1], exp(confint(glm_fits_comorbidity_gensim))[3, 2])
names(or_comorbidity_gensim_95ci) = c("ci_2.5", "ci_97.5")

## permutations ## 

## define functions
# define shuffle function - we shuffle the cancers associated (comorbid or comorbid and genetically similar) with each Mendelian disease
shuffle = function(matrix_to_be_shuffled) {
  
  matrix_shuffled = matrix_to_be_shuffled
  
  # shuffling
  for (i in 2:ncol(matrix_shuffled)) {
    matrix_shuffled[, i] = sample(matrix_shuffled[, i], replace = FALSE)
  }
  
  # convert to original format
  matrix_shuffled = tidyr::gather(matrix_shuffled, key = "mendelian_disease", value = "value", -cancer) %>%
    dplyr::select(mendelian_disease, cancer, value) %>%
    filter(value == 1) %>% 
    distinct()
  
  return(matrix_shuffled)
}

permutation = function(md_cancer_matrix, drugs_targeting_md_genes, drugs_nr_targets, investigated_indicated_drugs, md_drugs) {
  
  # shuffle comorbidity matrix
  md_cancer_matrix_shuffled = shuffle(md_cancer_matrix)
  
  # unique mendelian diseases and cancers in our sample
  unique_md = unique(md_cancer_matrix_shuffled$mendelian_disease)
  unique_cancers = unique(md_cancer_matrix_shuffled$cancer)
  
  # create logistic regression input - each row is a drug-cancer pair 
  log_input = data.frame(db_id = rep(drugs_targeting_md_genes$db_id, each = length(unique_cancers)))
  log_input$cancer = unique_cancers  
  # add information about number of targets for each drug
  log_input = left_join(log_input, drugs_nr_targets, by = "db_id")  
  # add information about investigated/indicated drugs within each drug-cancer pair
  log_input = left_join(log_input, investigated_indicated_drugs, by = c("cancer" = "complex_disease", 
                                                                        "db_id" = "drugbank_id"))
  # if indicated_investigated is NA --> this drug is not indicated/investigated for this cancer
  log_input$indicated_investigated = ifelse(is.na(log_input$indicated_investigated), 0, 1)
  
  # add information about comorbidity and genetic similarity support
  md_cancer_pairs_support = left_join(md_cancer_matrix_shuffled, md_drugs, by = "mendelian_disease")
  md_cancer_pairs_support = na.omit(md_cancer_pairs_support)
  md_cancer_pairs_support = md_cancer_pairs_support %>% 
    dplyr::select(db_id, cancer, value) %>%
    distinct()
  # some drug-cancer pairs might be supported by comorbidity from one MD-cancer pair but not another (same with genetic similarity)
  # for these pairs, keep the max evidence
  md_cancer_pairs_support = md_cancer_pairs_support %>% 
    group_by(db_id, cancer) %>%
    mutate(value = max(value)) %>%
    ungroup() %>%
    distinct()
  log_input = left_join(log_input, md_cancer_pairs_support, by = c("cancer", "db_id"))
  if (sum(is.na(log_input$value)) > 1) {
    log_input$value = ifelse(is.na(log_input$value), 0, 1)
  }
  
  ## run logistic regression ##
  glm_fits = glm(indicated_investigated ~ total_targets + value,
                 data = log_input, 
                 family = binomial())
  log_summary = summary(glm_fits)$coefficients
  log_summary
  # isolate needed results
  oddsratio_pvalue = c(exp(log_summary[3, 1]), log_summary[3, 4])
  names(oddsratio_pvalue) = c("odds_ratio", "pvalue")
  
  return(oddsratio_pvalue)
}

## to assess the significance of the comorbidity variable
md_cd_comorbidity_matrix = md_cancer_comorbidities %>% filter(comorbidity == 1)
md_cd_comorbidity_matrix = as.data.frame.matrix(table(md_cd_comorbidity_matrix[, 1:2]))
md_cd_comorbidity_matrix = md_cd_comorbidity_matrix %>% 
  mutate(cancer = rownames(md_cd_comorbidity_matrix)) %>%
  dplyr::select(cancer, everything())
rownames(md_cd_comorbidity_matrix) = NULL

# create empty vector to save the permutation results
or_comorbidity_perm = c()
pvalue_comorbidity_perm = c()

# run permutations
for (i in 1:1000) {
  
  # permutations
  perm_results = permutation(md_cancer_matrix = md_cd_comorbidity_matrix, 
                             drugs_targeting_md_genes, drugs_nr_targets, investigated_indicated_drugs, md_drugs)
  
  # update permutation results vectors
  or_comorbidity_perm[i] = perm_results[1]
  pvalue_comorbidity_perm[i] = perm_results[2]
  
  # track progress
  cat(i, "\n")
}

# calculate permutation p-values and odds ratio 95% CIs
# p-values
pvalue_comorbidity_perm_pvalue = sum(pvalue_comorbidity >= pvalue_comorbidity_perm) / 1000
pvalue_comorbidity_perm_pvalue # 0.007 (0.7%)
# odds ratios
or_comorbidity_perm_pvalue = sum(or_comorbidity <= or_comorbidity_perm) / 1000
or_comorbidity_perm_pvalue # 0.007 (0.7%)
or_comorbidity_perm_95ci = quantile(or_comorbidity_perm, c(0.05, 0.5, 0.95))

## to assess the significance of the comorbidity * genetic similarity
# create comorbidity matrix
md_cd_comorbidity_gensim_matrix = md_cancer_comorbidities %>% filter(comorbidity == 1 & genetic_similarity_melamed == 1)
md_cd_comorbidity_gensim_matrix = as.data.frame.matrix(table(md_cd_comorbidity_gensim_matrix[, 1:2]))
md_cd_comorbidity_gensim_matrix = md_cd_comorbidity_gensim_matrix %>% 
  mutate(cancer = rownames(md_cd_comorbidity_gensim_matrix)) %>%
  dplyr::select(cancer, everything())
rownames(md_cd_comorbidity_gensim_matrix) = NULL

# create empty vector to save the permutation results
or_comorbidity_gensim_perm = c()
pvalue_comorbidity_gensim_perm = c()

# run permutations
for (i in 1:1000) {
  
  # permutations
  perm_results = permutation(md_cancer_matrix = md_cd_comorbidity_gensim_matrix, 
                             drugs_targeting_md_genes, drugs_nr_targets, investigated_indicated_drugs, md_drugs)
  
  # update permutation results vectors
  or_comorbidity_gensim_perm[i] = perm_results[1]
  pvalue_comorbidity_gensim_perm[i] = perm_results[2]
  
  # track progress
  cat(i, "\n")
}

# calculate permutation p-values and odds ratio 95% CIs
# p-values
pvalue_comorbidity_gensim_perm_pvalue = sum(pvalue_comorbidity_gensim >= pvalue_comorbidity_gensim_perm) / 1000
pvalue_comorbidity_gensim_perm_pvalue
# odds ratios
or_comorbidity_gensim_perm_pvalue = sum(or_comorbidity_gensim <= or_comorbidity_gensim_perm) / 1000
or_comorbidity_gensim_perm_pvalue
or_comorbidity_gensim_perm_95ci = quantile(or_comorbidity_gensim_perm, c(0.05, 0.5, 0.95))

## create forestplot
data_fig4a = data.frame(category = c("Comorbidity & Genetic similarity", "Comorbidity"),
                        ci_2.5 = c(or_comorbidity_gensim_95ci["ci_2.5"], or_comorbidity_95ci["ci_2.5"]),
                        ci_97.5 = c(or_comorbidity_gensim_95ci["ci_97.5"], or_comorbidity_95ci["ci_97.5"]),
                        odds_ratio = c(or_comorbidity_gensim, or_comorbidity))
data_fig4a$category = factor(data_fig4a$category, levels = data_fig4a$category, labels = data_fig4a$category)

fig_4a = ggplot(data_fig4a, aes(y = category, x = odds_ratio, xmin = ci_2.5, xmax = ci_97.5)) +
  geom_point(size = 4) + 
  geom_errorbarh(height = 0.2, linewidth = 1.5) +
  ylab("") +
  xlab("Odds Ratio") +
  scale_x_continuous(breaks = seq(0, 8, 1), limits = c(0, 3)) +
  geom_vline(xintercept = 1, color = "black", linewidth = 1, linetype = 3) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 40, family = "Arial", color = "black", 
                                    margin = margin(t = 15)),
        axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.3),
        axis.text.y = element_text(size = 35, family = "Arial", color = "black",
                                   margin = margin(r = 15)),
        axis.text.x = element_text(size = 35, family = "Arial", color = "black"),
        legend.text = element_text(size = 30, family = "Arial", color = "black"),
        legend.title = element_blank(),
        legend.key.size = unit(2, "cm"))

fig_4a
ggsave(filename = "Fig4A_forest_plot.tiff", 
       path = "figures/",
       width = 15, height = 5, device = 'tiff',
       dpi = 300, compression = "lzw", type = type_compression)
dev.off()

####
#### We further explore the observed signal.
#### To do this, we need information about the genetic similarity between non-comorbid Mendelian disease - cancer pairs
#### However, Melamded et al. did not estimate it, as it was out of the scope of the paper
#### Therefore, we calculate the genetic similarity between any Mendelian disease - cancer pairs using two simple metrics:
####    1. Genetic overlap: tests if the same genes are mutated in a Mendelian disease - cancer pair
####    2. Co-expression: tests if the genes in a Mendelian disease - cancer pairs have similar expression across human tissues
####

#########################
##                     ##
##   Genetic overlap   ##
##                     ##
#########################

## load data
# MD-cancer comorbidities
md_cancer_comorbidities = fread("processed_data/md_cd_comorbidities.txt") %>%
  # filter for cancers with available data on TCGA
  filter(complex_disease %in% c("Female Breast Cancer", "Malignant Brain Neoplasm", "Melanoma", "Bladder Cancer", "Colorectal Cancer", "Gastric Cancer", "Kidney Cancer", "Lung Cancer", "Prostate Cancer", "Uterine Cancer"))

## Function to create list of TCGA significantly mutated genes for each cancer
sig_snv_cnv_combined = function(full_path_to_gistic2_ampl, full_path_to_gistic2_del, full_path_to_mutsig2cv, threshold_gistic2_qvalue = 0.05, threshold_gistic2_genes_in_peak = 50, threshold_mutsig2cv_qvalue = 0.05) {
  
  ## GISTIC2 amplifications
  gistic2_ampl = data.table::fread(full_path_to_gistic2_ampl)
  gistic2_ampl = gistic2_ampl[c(2,5:nrow(gistic2_ampl)), ]
  gistic2_ampl = as.data.frame(t(gistic2_ampl))
  gistic2_ampl = gistic2_ampl[-1, ]
  rownames(gistic2_ampl) = NULL
  gistic2_ampl$V1 = as.numeric(gistic2_ampl$V1)
  gistic2_ampl = gistic2_ampl %>% filter(V1 < threshold_gistic2_qvalue)
  gistic2_ampl = gistic2_ampl %>% mutate(peak = paste0("peak_", seq(1:nrow(gistic2_ampl))), .before = V1) %>% dplyr::select(-V1) %>% distinct()
  gistic2_ampl = reshape2::melt(gistic2_ampl, "peak", colnames(gistic2_ampl)[2:ncol(gistic2_ampl)])
  # remove empty cells
  gistic2_ampl = data.frame(gene = gistic2_ampl[-which(gistic2_ampl$value == ""), ])
  rownames(gistic2_ampl) = NULL
  # keep peaks with less than 50 genes --> following Rachel's paper
  gistic2_ampl = gistic2_ampl %>% 
    dplyr::select(peak = gene.peak, gene = gene.value) %>%
    distinct() %>% 
    group_by(peak) %>% 
    mutate(nr_genes_in_peak = length(gene)) %>% 
    ungroup() %>%
    filter(nr_genes_in_peak <= threshold_gistic2_genes_in_peak) %>%
    distinct() %>% 
    dplyr::select(gene) %>%
    distinct()
  # remove genes that contain "|" and then a specific transcript
  gistic2_ampl = data.frame(gene = gistic2_ampl[-grep("\\|", gistic2_ampl$gene), ])
  rownames(gistic2_ampl) = NULL
  gistic2_ampl = na.omit(gistic2_ampl)
  gistic2_ampl = unique(gistic2_ampl)
  
  cat("NOTE:",nrow(gistic2_ampl), "significantly amplified genes!\n")
  
  # GISTIC2 deletions
  gistic2_del = data.table::fread(full_path_to_gistic2_del)
  gistic2_del = gistic2_del[c(2,5:nrow(gistic2_del)), ]
  gistic2_del = as.data.frame(t(gistic2_del))
  gistic2_del = gistic2_del[-1, ]
  rownames(gistic2_del) = NULL
  gistic2_del$V1 = as.numeric(gistic2_del$V1)
  gistic2_del = gistic2_del %>% filter(V1 < threshold_gistic2_qvalue)
  gistic2_del = gistic2_del %>% mutate(peak = paste0("peak_", seq(1:nrow(gistic2_del))), .before = V1) %>% dplyr::select(-V1) %>% distinct()
  gistic2_del = reshape2::melt(gistic2_del, "peak", colnames(gistic2_del)[2:ncol(gistic2_del)])
  # remove empty cells
  gistic2_del = data.frame(gene = gistic2_del[-which(gistic2_del$value == ""), ])
  rownames(gistic2_del) = NULL
  # keep peaks with less than 50 genes
  gistic2_del = gistic2_del %>% 
    dplyr::select(peak = gene.peak, gene = gene.value) %>%
    distinct() %>% 
    group_by(peak) %>% 
    mutate(nr_genes_in_peak = length(gene)) %>% 
    ungroup() %>%
    filter(nr_genes_in_peak <= threshold_gistic2_genes_in_peak) %>%
    distinct() %>% 
    dplyr::select(gene) %>%
    distinct()
  # remove genes that contain "|" and then a specific transcript
  gistic2_del = data.frame(gene = gistic2_del[-grep("\\|", gistic2_del$gene), ])
  rownames(gistic2_del) = NULL
  gistic2_del = na.omit(gistic2_del)
  gistic2_del = unique(gistic2_del)
  
  cat("NOTE:", nrow(gistic2_del), "significantly deleted genes!\n")
  
  # Single nucleotide alterations
  mutsig2cv_sna = data.table::fread(full_path_to_mutsig2cv) %>%
    filter(q < threshold_mutsig2cv_qvalue) %>%
    dplyr::select(gene)
  
  cat("NOTE:", nrow(mutsig2cv_sna), "significantly mutated genes!\n")
  
  ## combine all of them together
  ultimate = rbind(gistic2_ampl, gistic2_del, mutsig2cv_sna)
  ultimate = unique(ultimate)
  
  cat("NOTE:", nrow(ultimate), "genes with significant SNV or CNV!")
  
  return(ultimate)
}

path_to_tcga_data = "raw_data/tcga_sig_mutated_genes/"
bladder_cancer = sig_snv_cnv_combined(full_path_to_gistic2_ampl = paste0(path_to_tcga_data, "BLCA_GISTIC2_amplifications.txt"),
                                      full_path_to_gistic2_del = paste0(path_to_tcga_data, "BLCA_GISTIC2_deletions.txt"),
                                      full_path_to_mutsig2cv = paste0(path_to_tcga_data, "BLCA_mutsig2cv.txt"), 
                                      threshold_gistic2_qvalue = 0.05, threshold_gistic2_genes_in_peak = 50, threshold_mutsig2cv_qvalue = 0.05)
female_breast_cancer = sig_snv_cnv_combined(full_path_to_gistic2_ampl = paste0(path_to_tcga_data, "BRCA_GISTIC2_amplifications.txt"),
                                            full_path_to_gistic2_del = paste0(path_to_tcga_data, "BRCA_GISTIC2_deletions.txt"),
                                            full_path_to_mutsig2cv = paste0(path_to_tcga_data, "BRCA_mutsig2cv.txt"), 
                                            threshold_gistic2_qvalue = 0.05, threshold_gistic2_genes_in_peak = 50, threshold_mutsig2cv_qvalue = 0.05)
colorectal_cancer = sig_snv_cnv_combined(full_path_to_gistic2_ampl = paste0(path_to_tcga_data, "COAD_GISTIC2_amplifications.txt"),
                                         full_path_to_gistic2_del = paste0(path_to_tcga_data, "COAD_GISTIC2_deletions.txt"),
                                         full_path_to_mutsig2cv = paste0(path_to_tcga_data, "COAD_mutsig2cv.txt"), 
                                         threshold_gistic2_qvalue = 0.05, threshold_gistic2_genes_in_peak = 50, threshold_mutsig2cv_qvalue = 0.05) # COAD is the most common gastric cancer (compared to READ)
gastric_cancer = sig_snv_cnv_combined(full_path_to_gistic2_ampl = paste0(path_to_tcga_data, "STAD_GISTIC2_amplifications.txt"),
                                      full_path_to_gistic2_del = paste0(path_to_tcga_data, "STAD_GISTIC2_deletions.txt"),
                                      full_path_to_mutsig2cv = paste0(path_to_tcga_data, "STAD_mutsig2cv.txt"), 
                                      threshold_gistic2_qvalue = 0.05, threshold_gistic2_genes_in_peak = 50, threshold_mutsig2cv_qvalue = 0.05)
kidney_cancer = sig_snv_cnv_combined(full_path_to_gistic2_ampl = paste0(path_to_tcga_data, "KIRC_GISTIC2_amplifications.txt"),
                                     full_path_to_gistic2_del = paste0(path_to_tcga_data, "KIRC_GISTIC2_deletions.txt"),
                                     full_path_to_mutsig2cv = paste0(path_to_tcga_data, "KIRC_mutsig2cv.txt"), 
                                     threshold_gistic2_qvalue = 0.05, threshold_gistic2_genes_in_peak = 50, threshold_mutsig2cv_qvalue = 0.05) # KIRC is more common kidney cancer (compared to KICH, KIRP)
lung_cancer = sig_snv_cnv_combined(full_path_to_gistic2_ampl = paste0(path_to_tcga_data, "LUAD_GISTIC2_amplifications.txt"),
                                   full_path_to_gistic2_del = paste0(path_to_tcga_data, "LUAD_GISTIC2_deletions.txt"),
                                   full_path_to_mutsig2cv = paste0(path_to_tcga_data, "LUAD_mutsig2cv.txt"), 
                                   threshold_gistic2_qvalue = 0.05, threshold_gistic2_genes_in_peak = 50, threshold_mutsig2cv_qvalue = 0.05) # LUAD is more common lung cancer (compared to LUSC)
malignant_brain_neoplasm = sig_snv_cnv_combined(full_path_to_gistic2_ampl = paste0(path_to_tcga_data, "GBM_GISTIC2_amplifications.txt"),
                                                full_path_to_gistic2_del = paste0(path_to_tcga_data, "GBM_GISTIC2_deletions.txt"),
                                                full_path_to_mutsig2cv = paste0(path_to_tcga_data, "GBM_mutsig2cv.txt"), 
                                                threshold_gistic2_qvalue = 0.05, threshold_gistic2_genes_in_peak = 50, threshold_mutsig2cv_qvalue = 0.05) # GBM is more common brain cancer (compared to LGG)
melanoma = sig_snv_cnv_combined(full_path_to_gistic2_ampl = paste0(path_to_tcga_data, "SKCM_GISTIC2_amplifications.txt"),
                                full_path_to_gistic2_del = paste0(path_to_tcga_data, "SKCM_GISTIC2_deletions.txt"),
                                full_path_to_mutsig2cv = paste0(path_to_tcga_data, "SKCM_mutsig2cv.txt"), 
                                threshold_gistic2_qvalue = 0.05, threshold_gistic2_genes_in_peak = 50, threshold_mutsig2cv_qvalue = 0.05)
prostate_cancer = sig_snv_cnv_combined(full_path_to_gistic2_ampl = paste0(path_to_tcga_data, "PRAD_GISTIC2_amplifications.txt"),
                                       full_path_to_gistic2_del = paste0(path_to_tcga_data, "PRAD_GISTIC2_deletions.txt"),
                                       full_path_to_mutsig2cv = paste0(path_to_tcga_data, "PRAD_mutsig2cv.txt"), 
                                       threshold_gistic2_qvalue = 0.05, threshold_gistic2_genes_in_peak = 50, threshold_mutsig2cv_qvalue = 0.05)
uterine_cancer = sig_snv_cnv_combined(full_path_to_gistic2_ampl = paste0(path_to_tcga_data, "UCEC_GISTIC2_amplifications.txt"),
                                      full_path_to_gistic2_del = paste0(path_to_tcga_data, "UCEC_GISTIC2_deletions.txt"),
                                      full_path_to_mutsig2cv = paste0(path_to_tcga_data, "UCEC_mutsig2cv.txt"), 
                                      threshold_gistic2_qvalue = 0.05, threshold_gistic2_genes_in_peak = 50, threshold_mutsig2cv_qvalue = 0.05)

## calculate the significance of genetic overlap for all possible MD-cancer pairs
unique_cancers
mendelian_diseases = unique(md_cancer_comorbidities$mendelian_disease) # 60 Mendelian diseases

# genes linked to the 60 Mendelian diseases
md_genes = md_genes %>%
  filter(mendelian_disease %in% mendelian_diseases)

# protein coding genes
# link: https://www.genenames.org/download/statistics-and-files/
protein_coding_genes = fread("https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt") %>%
  dplyr::select(gene_name = symbol, prev_symbol) %>%
  distinct()

# define function
genetic_overlap = function(md, cancer) {
  
  ## Mendelian disease genes
  md_genes_temp = md_genes %>% filter(mendelian_disease == md)
  md_genes_temp = unique(md_genes_temp$causal_gene)
  # filter for protein coding genes
  md_genes_temp = intersect(md_genes_temp, unique(c(protein_coding_genes$gene_name, protein_coding_genes$prev_symbol)))
  
  ## Cancer significantly altered genes
  cancer_genes_temp = get(tolower(gsub(pattern = " ", replacement = "_", cancer)))
  cancer_genes_temp = unique(cancer_genes_temp$gene)
  # filter for protein coding genes
  cancer_genes_temp = intersect(cancer_genes_temp, unique(c(protein_coding_genes$gene_name, protein_coding_genes$prev_symbol)))
  
  ## data frame to populate later
  fisher_pvalues = data.frame(mendelian_disease = NA, cancer = NA, pvalue_onesided = NA)
  
  ## create Fisher's test input
  contingency_table = matrix(nrow = 2, ncol = 2)
  colnames(contingency_table) = c("cancer gene", "non cancer gene")
  rownames(contingency_table) = c("md gene", "non md gene")
  
  ## background genes are all protein coding genes
  contingency_table[1, 1] = length(intersect(md_genes_temp, cancer_genes_temp))
  contingency_table[2, 1] = length(setdiff(cancer_genes_temp, md_genes_temp))
  contingency_table[1, 2] = length(setdiff(md_genes_temp, cancer_genes_temp))
  contingency_table[2, 2] = nrow(protein_coding_genes) - contingency_table[1, 1] - contingency_table[1, 2] - contingency_table[2, 1]
  
  ## Fisher's test
  fisher_test_onesided = fisher.test(contingency_table, alternative = "greater")
  
  ## Complete data frame
  fisher_pvalues[1, "mendelian_disease"] = md
  fisher_pvalues[1, "cancer"] = cancer
  fisher_pvalues[1, "pvalue_onesided"] = fisher_test_onesided$p.value
  
  return(fisher_pvalues)
}

genetic_overlap_results = data.frame(mendelian_disease = NA, cancer = NA, pvalue_onesided = NA)
for (mendelian_disease_temp in mendelian_diseases) {
  for (cancer_temp in unique_cancers) {
    # estimate genetic similarity
    temp = genetic_overlap(md = mendelian_disease_temp, cancer = cancer_temp)
    genetic_overlap_results = rbind(genetic_overlap_results, temp)
    
    cat(mendelian_disease_temp, "-", cancer_temp, "\n")
  }
}
rm(mendelian_disease_temp, cancer_temp, temp)
genetic_overlap_results = genetic_overlap_results[-1, ] ; rownames(genetic_overlap_results) = NULL

## annotate with comorbidity
md_cancer_comorbidities$comorbidity = 1
genetic_overlap_results = left_join(genetic_overlap_results, md_cancer_comorbidities, by = c("mendelian_disease", "cancer" = "complex_disease"))
genetic_overlap_results$comorbidity = ifelse(is.na(genetic_overlap_results$comorbidity), 0, 1)

## annotate with genetic similarity
genetic_overlap_results$genetic_similarity = ifelse(genetic_overlap_results$pvalue_onesided < 0.05, 1, 0)

## save
fwrite(genetic_overlap_results, "processed_data/md_cancers_genetic_overlap.txt", sep = "\t", row.names = FALSE)

#########################
##                     ##
##    Co-expression    ##
##                     ##
#########################

### Co-expression analysis between MD-cancer comorbidity pairs ###

## load GTEx normalized TPM per gene for all human tissues
## downloaded from Human Protein Atlas: https://www.proteinatlas.org/download/rna_tissue_gtex.tsv.zip
gtex_ntpm = fread("raw_data/rna_tissue_gtex.tsv")
gtex_ntpm = gtex_ntpm %>% 
  dplyr::select(gene_name = "Gene name", tissue = Tissue, ntpm = nTPM) %>%
  distinct()

# calculate correlations of gene expression between each cancer gene and a set of MD genes
# put genes in the columns and tissues in the rows
gtex_ntpm = reshape(gtex_ntpm, idvar = "gene_name", timevar = "tissue", direction = "wide") # ignore warnings, they have to do with the NAs
# some genes have expression measured in not all the 37 tissues --> they have NA values --> remove them
gtex_ntpm = na.omit(gtex_ntpm)
to_name_columns = gtex_ntpm$gene_name
# remove first column with gene names
gtex_ntpm = gtex_ntpm[, -1]
str(gtex_ntpm) # all values are numerics
# transpose table so correlation to be measured gene-gene
gtex_ntpm = t(gtex_ntpm)
# name columns
colnames(gtex_ntpm) = to_name_columns ; rm(to_name_columns)
# calculate pearson correlation
pearson_correlation_gtex = cor(gtex_ntpm, gtex_ntpm, method = "pearson")
# convert it to data frame
pearson_correlation_gtex = as.data.frame(pearson_correlation_gtex)
# remove NA values - occurred because expression of some genes were always 0
pearson_correlation_gtex = pearson_correlation_gtex[-which(is.na(pearson_correlation_gtex$TSPAN6)), ] # from rows
pearson_correlation_gtex = pearson_correlation_gtex[, -which(is.na(pearson_correlation_gtex[1, ]))] # from columns

## define function to calculate co-expression
coexpression = function(md, cancer) {
  
  ## Mendelian disease genes
  md_genes_temp = md_genes %>% 
    filter(mendelian_disease == md) %>%
    filter(causal_gene %in% colnames(pearson_correlation_gtex))
  md_genes_temp = unique(md_genes_temp$causal_gene)
  
  ## Cancer significantly altered genes
  cancer_genes_temp = get(tolower(gsub(pattern = " ", replacement = "_", cancer))) %>% 
    filter(gene %in% colnames(pearson_correlation_gtex))
  cancer_genes_temp = unique(cancer_genes_temp$gene)
  
  ## gene expression data
  pvalues_per_cancer_gene = data.frame(cancer_gene = NA, p_value = NA)
  
  # for each cancer gene...
  for (i in 1:length(cancer_genes_temp)) {
    
    # cancer gene 
    cancer_gene_temp = cancer_genes_temp[i]
    
    # non-MD genes
    non_md_genes = colnames(pearson_correlation_gtex)
    non_md_genes = non_md_genes[-which(non_md_genes %in% md_genes_temp)]
    
    # co-expression of cancer gene - MD genes
    cancer_md_expression_matrix_temp = pearson_correlation_gtex[cancer_gene_temp, md_genes_temp]
    cancer_md_expression_matrix_temp = unlist(c(cancer_md_expression_matrix_temp), use.names = FALSE)
    
    # co-expression of cancer gene - non-MD genes
    cancer_nonmd_expression_matrix_temp = pearson_correlation_gtex[cancer_gene_temp, non_md_genes]
    cancer_nonmd_expression_matrix_temp = unlist(c(cancer_nonmd_expression_matrix_temp), use.names = FALSE)
    
    # wilcoxon rank sum test -- testing if the co-expression between one cancer gene and a set of MD genes are higher than the co-expression of one cancer gene and a set of nonMD genes
    rank_sum_test = wilcox.test(cancer_md_expression_matrix_temp, cancer_nonmd_expression_matrix_temp, alternative = "greater")
    
    # populate data frame
    pvalues_per_cancer_gene[i, "cancer_gene"] = cancer_gene_temp
    pvalues_per_cancer_gene[i, "p_value"] = rank_sum_test$p.value
  }
  
  pvalues_per_cancer_gene$mendelian_disease = md
  pvalues_per_cancer_gene$cancer = cancer
  
  return(pvalues_per_cancer_gene)
}

coexpression_results = data.frame(cancer_gene = NA, p_value = NA, cancer = NA, mendelian_disease = NA)
for (mendelian_disease_temp in mendelian_diseases) {
  for (cancer_temp in unique_cancers) {
    # estimate genetic similarity
    temp = coexpression(md = mendelian_disease_temp, cancer = cancer_temp)
    coexpression_results = rbind(coexpression_results, temp)
    
    cat(mendelian_disease_temp, "-", cancer_temp, "\n")
  }
}

# Familial Dysautonomia has one gene but this is not found in the GTEx dataset --> couldn't be tested with this genetic similarity metric
coexpression_results = coexpression_results[-1, ]
rownames(ultimate) = NULL

## correct p_values for the number of cancer genes tested in each Mendelian disease - cancer pair
coexpression_results = coexpression_results %>% 
  group_by(cancer, mendelian_disease) %>%
  mutate(nr_cancer_genes = length(cancer_gene)) %>%
  mutate(adj_pvalue_threshold = p.adjust(p_value, method = "BH", n = unique(nr_cancer_genes))) %>%
  mutate(adj_pvalue_threshold = min(adj_pvalue_threshold)) %>%
  ungroup() %>%
  dplyr::select(mendelian_disease, cancer, adj_pvalue_threshold) %>%
  distinct()
coexpression_results$genetic_similarity = ifelse(coexpression_results$adj_pvalue_threshold < 0.05, 1, 0)
coexpression_results = coexpression_results %>%
  dplyr::select(cancer, mendelian_disease, adj_pvalue_threshold, genetic_similarity) %>%
  distinct()

## annotate with comorbidity
coexpression_results = left_join(coexpression_results, md_cancer_comorbidities, by = c("mendelian_disease", "cancer" = "complex_disease"))
coexpression_results$comorbidity = ifelse(is.na(coexpression_results$comorbidity), 0, 1)

# save
fwrite(coexpression_results, "processed_data/md_cancers_coexpression.txt", sep = "\t", row.names = FALSE)

#########################
##                     ##
##   Combine metrics   ##
##                     ##
#########################

md_cd_comorbidities_gensim_ultimate = full_join(genetic_overlap_results, coexpression_results, by = c("cancer", "mendelian_disease", "comorbidity"))
colnames(md_cd_comorbidities_gensim_ultimate) = c("mendelian_disease", "cancer", "pvalue_onesided_geneoverlap", "comorbidity", "genetic_similarity_geneoverlap", "adj_pvalue_coexpression", "genetic_similarity_coexpression")
# fix information for Familial Dysautonomia
md_cd_comorbidities_gensim_ultimate[541:550, "genetic_similarity_coexpression"] = 0

# sum up values of genetic_similarity gene overlap and co-expression columns
md_cd_comorbidities_gensim_ultimate$gensim_combined = md_cd_comorbidities_gensim_ultimate$genetic_similarity_geneoverlap + md_cd_comorbidities_gensim_ultimate$genetic_similarity_coexpression
summary(md_cd_comorbidities_gensim_ultimate$gensim_combined)

# if gensim_combined = 0 --> no evidence for genetic similarity between a MD-cancer pair
# if gensim_combined = 1 --> evidence for genetic similarity between a MD-cancer pair from only one metric
# if gensim_combined = 2 --> evidence for genetic similarity between a MD-cancer pair from both metrics
# based on the above, replace this column with 0/1 where 1 indicates evidence for genetic similarity between a MD-cancer from at least one metric
md_cd_comorbidities_gensim_ultimate$gensim_combined = ifelse(md_cd_comorbidities_gensim_ultimate$gensim_combined == 0, 0, 1)
table(md_cd_comorbidities_gensim_ultimate$gensim_combined)
md_cd_comorbidities_gensim_ultimate = md_cd_comorbidities_gensim_ultimate %>% 
  dplyr::select(cancer, mendelian_disease, comorbidity, gensim_combined) %>% 
  distinct()

#### --- visualizations --- ####

# create stacked barplot showing the number of pairs with/without comorbidity and with/without genetic similarity
table(md_cd_comorbidities_gensim_ultimate$comorbidity, md_cd_comorbidities_gensim_ultimate$gensim_combined)
# rows are comorbidity (no/yes)
# columns are genetic similarity (no/yes)
data_fig4b = data.frame(comorbidity = c("Comorbid", "Comorbid", "Non comorbid", "Non comorbid"),
                        gen_sim = c( "Yes", "No", "Yes", "No"),
                        number_of_pairs = c(table(md_cd_comorbidities_gensim_ultimate$comorbidity, md_cd_comorbidities_gensim_ultimate$gensim_combined)[2, 2],
                                            table(md_cd_comorbidities_gensim_ultimate$comorbidity, md_cd_comorbidities_gensim_ultimate$gensim_combined)[2, 1],
                                            table(md_cd_comorbidities_gensim_ultimate$comorbidity, md_cd_comorbidities_gensim_ultimate$gensim_combined)[1, 2],
                                            table(md_cd_comorbidities_gensim_ultimate$comorbidity, md_cd_comorbidities_gensim_ultimate$gensim_combined)[1, 1]))
data_fig4b$comorbidity = factor(data_fig4b$comorbidity, levels = data_fig4b$comorbidity, labels = data_fig4b$comorbidity)
data_fig4b$gen_sim = factor(data_fig4b$gen_sim, levels = data_fig4b$gen_sim, labels = data_fig4b$gen_sim)

# Stacked barplot
fig_4b = ggplot(data_fig4b, aes(x = comorbidity, y = number_of_pairs, fill=gen_sim)) + 
  geom_bar(position = "stack", stat = "identity", width = 0.5) +
  labs(fill = "Genetically similar") +
  xlab("") +
  ylab("Number of Mendelian disease - cancer pairs\n") +
  scale_y_continuous(breaks = seq(0, 320, 40)) +
  scale_fill_manual(values=c("#F9E855", "#43347A")) +
  theme_classic() +
  theme(axis.title = element_text(angle = 0, hjust = 0.5, face = "bold", 
                                  margin = margin(t = 3, unit = "cm"),
                                  size = 24, family = "Arial", color = "black"),
        axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.3),
        axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,
                                   margin = margin(l = 0.5, r = 0.2, unit = "cm"),
                                   size = 24, family = "Arial", color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, face = "bold",
                                   margin = margin(t = 0.2, unit = "cm"),
                                   size = 24, family = "Arial", color = "black"),
        legend.text = element_text(size = 20, family = "Arial", color = "black"),
        legend.title = element_text(size = 24, family = "Arial", color = "black"),
        aspect.ratio = 1.3)

fig_4b
ggsave(filename = "Fig4B_md_cancer_comorbidity_gensim_stacked_barplot.tiff", 
       path = "figures/",
       width = 12, height = 9, device = 'tiff',
       dpi = 700, compression = "lzw", type = type_compression)
dev.off()

## check if there is a significant enrichment of genetically similar Mendelian disease - cancer pairs based on our genetic similarity metric and Melamed et al. metrics
contingency_table = matrix(nrow = 2, ncol = 2)
colnames(contingency_table) = c("melamed_sig", "melamed_not_sig")
rownames(contingency_table) = c("ourmetrics_sig", "ourmetrics_notsig")
# across comorbidity pairs only (remember, Melamed et al. did not estimate genetic similarity between non-comorbid Mendelian disease - cancer pairs)
comorbid_pairs = md_cancer_comorbidities
comorbid_pairs = left_join(comorbid_pairs, md_cancer_gensim_melamed, by = c("mendelian_disease", "complex_disease"))
comorbid_pairs$genetic_similarity_melamed = ifelse(is.na(comorbid_pairs$genetic_similarity_melamed), 0, comorbid_pairs$genetic_similarity_melamed)
comorbid_pairs = left_join(comorbid_pairs, md_cd_comorbidities_gensim_ultimate, by = c("mendelian_disease", "complex_disease" = "cancer", "comorbidity"))
table(comorbid_pairs$genetic_similarity_melamed, comorbid_pairs$gensim_combined)
contingency_table[1, 1] = 64
contingency_table[2, 1] = 34
contingency_table[1, 2] = 71
contingency_table[2, 2] = 145
fisher_ovelap_melamed_our_gensim_metrics_comorbid_pairs = fisher.test(contingency_table, alternative = "greater") 
fisher_ovelap_melamed_our_gensim_metrics_comorbid_pairs
rm(contingency_table, comorbid_pairs)

## -- logistic regression analyses -- ##

## create logistic regression input ##

# unique Mendelian diseases and cancers in our sample
mendelian_diseases
unique_cancers

# unique drugs targeting the genes associated with the 60 Mendelian diseases
drugs_targeting_md_genes

# create logistic regression input - each row is a drug-cancer pair 
log_input = data.frame(db_id = rep(drugs_targeting_md_genes$db_id, each = length(unique_cancers)))
log_input$cancer = unique_cancers  
# add information about number of targets for each drug
log_input = left_join(log_input, drugs_nr_targets, by = "db_id")  
# add information about investigated/indicated drugs within each drug-cancer pair
log_input = left_join(log_input, investigated_indicated_drugs, by = c("cancer" = "complex_disease", 
                                                                      "db_id" = "drugbank_id"))
# if indicated_investigated is NA --> this drug is not indicated/investigated for this cancer
log_input$indicated_investigated = ifelse(is.na(log_input$indicated_investigated), 0, 1)

# add information about comorbidity and genetic similarity support
md_cancer_pairs_drugs_support = left_join(md_cd_comorbidities_gensim_ultimate, md_drugs, by = "mendelian_disease")
md_cancer_pairs_drugs_support = na.omit(md_cancer_pairs_drugs_support) ; rownames(md_cancer_pairs_drugs_support) = NULL
md_cancer_pairs_drugs_support = md_cancer_pairs_drugs_support %>% 
  dplyr::select(db_id, cancer, comorbidity, gensim_combined) %>%
  distinct()
# some drug-cancer pairs might be supported by comorbidity from one MD-cancer pair but not another (same with genetic similarity)
# for these pairs, keep the max evidence
md_cancer_pairs_drugs_support = md_cancer_pairs_drugs_support %>% 
  group_by(db_id, cancer) %>%
  mutate(comorbidity = max(comorbidity),
         gensim_combined = max(gensim_combined)) %>%
  ungroup() %>%
  distinct()
log_input = left_join(log_input, md_cancer_pairs_drugs_support, by = c("cancer", "db_id"))

## -- run logistic regression -- ##
## analysis per drug-cancer pairs ##

# filter for drug-cancer pairs without comorbidity support
log_input_no_comorbidity = log_input %>% filter(comorbidity == 0)
glm_fits_no_comorbidity = glm(indicated_investigated ~ total_targets + gensim_combined,
                              data = log_input_no_comorbidity, 
                              family = binomial())
summary(glm_fits_no_comorbidity)$coefficients

# filter for drug-cancer pairs with comorbidity support
log_input_comorbidity = log_input %>% filter(comorbidity == 1)
glm_fits_comorbidity = glm(indicated_investigated ~ total_targets + gensim_combined,
                           data = log_input_comorbidity, 
                           family = binomial())
summary(glm_fits_comorbidity)$coefficients

# barplot with actual number of drug candidates that are currently investigated/indicated
# comorbidity
log_input_comorbidity
sum(log_input_comorbidity$indicated_investigated) / nrow(log_input_comorbidity) # 6.8%

# comorbidity - NO genetic similarity
log_input_comorbidity_no_gensim = log_input_comorbidity %>% filter(gensim_combined == 0)
sum(log_input_comorbidity_no_gensim$indicated_investigated) / nrow(log_input_comorbidity_no_gensim) # 5.6%

# comorbidity - genetic similarity
log_input_comorbidity_gensim = log_input_comorbidity %>% filter(gensim_combined == 1)
sum(log_input_comorbidity_gensim$indicated_investigated) / nrow(log_input_comorbidity_gensim) # 7.8%

# NO comorbidity
log_input_no_comorbidity
sum(log_input_no_comorbidity$indicated_investigated) / nrow(log_input_no_comorbidity) # 4.0%

# NO comorbidity - NO genetic similarity
log_input_no_comorbidity_no_gen_sim = log_input_no_comorbidity %>% filter(gensim_combined == 0)
sum(log_input_no_comorbidity_no_gen_sim$indicated_investigated) / nrow(log_input_no_comorbidity_no_gen_sim) # 2.4%

# NO comorbidity - genetic similarity
log_input_no_comorbidity_gen_sim = log_input_no_comorbidity %>% filter(gensim_combined == 1)
sum(log_input_no_comorbidity_gen_sim$indicated_investigated) / nrow(log_input_no_comorbidity_gen_sim) # 5.5%

## percentages
data_fig4c_fraction = data.frame(category = c("Comorbidity \n& \nGenetic similarity",
                                              "Comorbidity\n\n",
                                              "Comorbidity \n& \nNo genetic similarity",
                                              "No comorbidity \n& \nGenetic similarity",
                                              "No comorbidity\n\n",
                                              "No comorbidity \n& \nNo genetic similarity"),
                                 nr_inv_ind_drugs = c(sum(log_input_comorbidity_gensim$indicated_investigated),
                                                      sum(log_input_comorbidity$indicated_investigated),
                                                      sum(log_input_comorbidity_no_gensim$indicated_investigated),
                                                      sum(log_input_no_comorbidity_gen_sim$indicated_investigated),
                                                      sum(log_input_no_comorbidity$indicated_investigated),
                                                      sum(log_input_no_comorbidity_no_gen_sim$indicated_investigated)),
                                 nr_candidate_drugs = c(nrow(log_input_comorbidity_gensim),
                                                        nrow(log_input_comorbidity),
                                                        nrow(log_input_comorbidity_no_gensim),
                                                        nrow(log_input_no_comorbidity_gen_sim),
                                                        nrow(log_input_no_comorbidity),
                                                        nrow(log_input_no_comorbidity_no_gen_sim)))
data_fig4c_fraction$fraction = data_fig4c_fraction$nr_inv_ind_drugs / data_fig4c_fraction$nr_candidate_drugs
data_fig4c_fraction$category = factor(data_fig4c_fraction$category, levels = data_fig4c_fraction$category, labels = data_fig4c_fraction$category)

fig_4c_percentages = ggplot(data_fig4c_fraction, aes(x = category, y = fraction)) +
  geom_col(width = 0.6) +
  labs(fill = "") +
  ylab("Fraction of successful candidate drugs") +
  xlab("") +
  scale_y_continuous(breaks = seq(0, 1, 0.01)) +
  theme_classic() +
  theme(panel.background = element_rect(fill = "white", 
                                        colour = "white"),
        axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5,
                                   margin = margin(l = 0.5, r = 0.2, unit = "cm"),
                                   size = 20, family = "Arial", color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, face = "bold",
                                   margin = margin(t = 0.2, unit = "cm"),
                                   size = 20, family = "Arial", color = "black"),
        axis.title = element_text(angle = 0, hjust = 0.5, face = "bold", 
                                  margin = margin(t = 3, unit = "cm"),
                                  size = 20, family = "Arial", color = "black"),
        axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.3))

fig_4c_percentages
ggsave(filename = "Fig4C_fraction_candidate_drugs_inv_ind.tiff", 
       path = "figures/",
       width = 20, height = 12, device = 'tiff',
       dpi = 700, compression = "lzw", type = type_compression)
dev.off()

## save results
odds_ratio_genetic_similarity_melamed = data.frame(disease_category = c("Comorbidity", "Comorbidity & Genetic similarity"),
                                                   odds_ratio = c(or_comorbidity, or_comorbidity_gensim),
                                                   pvalue = c(pvalue_comorbidity, pvalue_comorbidity_gensim),
                                                   ci_2.5 = c(or_comorbidity_95ci["ci_2.5"], or_comorbidity_gensim_95ci["ci_2.5"]),
                                                   ci_97.5 = c(or_comorbidity_95ci["ci_97.5"], or_comorbidity_gensim_95ci["ci_97.5"]),
                                                   ci_5_permutations = c(or_comorbidity_perm_95ci["5%"], or_comorbidity_gensim_perm_95ci["5%"]),
                                                   ci_50_permutations = c(or_comorbidity_perm_95ci["50%"], or_comorbidity_gensim_perm_95ci["50%"]),
                                                   ci_95_permutations = c(or_comorbidity_perm_95ci["95%"], or_comorbidity_gensim_perm_95ci["95%"]))
fwrite(odds_ratio_genetic_similarity_melamed, "results/log_reg_results_comorbidity_genetic_similarity_odds_ratio_fig4a.txt", sep = "\t", row.names = FALSE)
fwrite(data_fig4c_fraction, "results/fraction_candidate_drugs_inv_ind_comorbidity_genetic_similarity.txt", sep = "\t", row.names = FALSE)
