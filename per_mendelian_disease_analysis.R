### This script implements the per Mendelian disease analysis ###

library(dplyr)
library(data.table)
library(ggplot2)
library(ggsignif)

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

## -- load data -- ##

# Mendelian - complex disease comorbidities
md_cd_comorbidities = fread("processed_data/md_cd_comorbidities.txt")
unique_cd = unique(md_cd_comorbidities$complex_disease)
unique_md = unique(md_cd_comorbidities$mendelian_disease)

# complex disease categories
complex_disease_categories = data.table::fread("processed_data/complex_disease_category.txt")
md_cd_comorbidities = left_join(md_cd_comorbidities, complex_disease_categories, by = "complex_disease")

# indicated/investigated drugs for the complex diseases
investigated_indicated_drugs = fread("processed_data/drugs_inv_ind_per_disease.txt")
investigated_indicated_drugs = reshape2::melt(investigated_indicated_drugs, "drugbank_id", colnames(investigated_indicated_drugs)[2:ncol(investigated_indicated_drugs)])
investigated_indicated_drugs = na.omit(investigated_indicated_drugs) ; rownames(investigated_indicated_drugs) = NULL
investigated_indicated_drugs = investigated_indicated_drugs %>%
  dplyr::select(drugbank_id, complex_disease = variable) %>%
  distinct()

# Mendelian disease causal genes
md_genes = data.table::fread("processed_data/md_genes.txt")

# drug - targets
db_drug_targets = fread("processed_data/drugbank_all_drugs_known_targets.txt") %>%
  dplyr::select(db_id, drug_target) %>% 
  distinct()

# number of gene-targets per drug
drugs_nr_targets = fread("processed_data/drugbank_all_drugs_known_targets.txt") %>%
  dplyr::select(db_id, drug_target) %>%
  distinct() %>%
  group_by(db_id) %>% 
  mutate(total_targets = length(drug_target)) %>%
  dplyr::select(db_id, total_targets) %>% 
  ungroup() %>%
  distinct()

## -- find drugs targeting each Mendelian disease genes -- ##

# empty list - each element is a mendelian disease
md_comorbidities = list()
for (md in unique_md) {
  md_comorbidities[[md]] = md_cd_comorbidities %>% filter(mendelian_disease == md) %>% dplyr::select(-comorbidity)
} ; rm(md)

# annotate with Mendelian disease genes
for (i in 1:length(md_comorbidities)) {
  
  # add Mendelian disease genes
  md_comorbidities[[i]] = left_join(md_comorbidities[[i]], md_genes, by = "mendelian_disease")
  md_comorbidities[[i]] = md_comorbidities[[i]] %>% 
    dplyr::select(complex_disease, disease_category, md_gene = causal_gene) %>%
    distinct()
  
  # add drugs targeting the Mendelian disease genes
  md_comorbidities[[i]] = left_join(md_comorbidities[[i]], db_drug_targets, by = c("md_gene" = "drug_target"))
  if (sum(is.na(md_comorbidities[[i]]$db_id)) != 0) {
    md_comorbidities[[i]] = md_comorbidities[[i]][-which(is.na(md_comorbidities[[i]]$db_id)), ] ; rownames(md_comorbidities[[i]]) = NULL
  }
  
  # keep needed columns
  md_comorbidities[[i]] = md_comorbidities[[i]] %>% 
    dplyr::select(complex_disease, disease_category, db_id) %>%
    distinct()
  
  cat(i, "\n")
} ; rm(i)

### NOTE: 19 Mendelian diseases are causally associated with genes NOT targeted by any existing drug with know drug targets
###       Therefore, we remove them from the downstream analysis
md_comorbidities = md_comorbidities[lapply(md_comorbidities, nrow) > 0] # 71 Mendelian diseases remain

## -- logistic regression per Mendelian disease -- ##

# create list of logistic regression inputs for each Mendelian disease
log_reg_input = vector("list", length = length(md_comorbidities))
names(log_reg_input) = names(md_comorbidities)

for (i in 1:length(md_comorbidities)) {
  
  # if drugs are recommended for repurposing...
  if (nrow(md_comorbidities[[i]]) != 0) {
    
    ## candidate drugs
    candidate_drugs_temp = unique(md_comorbidities[[i]]$db_id)
    
    ## complex diseases that these drugs are recommended for repurposing
    complex_diseases_temp = unique(md_comorbidities[[i]]$complex_disease)
    
    ## complex diseases these drugs are indicated/investigated
    # not all drugs will be indicated/investigated for a disease --> keep this in mind when you build the logistic regression input data frame
    ind_inv_temp = investigated_indicated_drugs %>%
      filter(drugbank_id %in% candidate_drugs_temp) %>%
      mutate(indicated_investigated = 1)
    
    ## create logistic regression input
    
    # each row should be a drug-complex disease pair
    log_reg_input[[i]] = data.frame(drug = rep(candidate_drugs_temp, each = length(unique_cd)))
    
    # add complex diseases
    log_reg_input[[i]]$disease = unique_cd

    # add disease category
    log_reg_input[[i]] = left_join(log_reg_input[[i]], complex_disease_categories, by = c("disease" = "complex_disease"))
    
    # add number of targets
    log_reg_input[[i]] = left_join(log_reg_input[[i]], drugs_nr_targets, by = c("drug" = "db_id"))
    
    # add indicated/investigated drugs
    log_reg_input[[i]] = left_join(log_reg_input[[i]], ind_inv_temp, by = c("drug" = "drugbank_id", "disease" = "complex_disease"))
    log_reg_input[[i]]$indicated_investigated = ifelse(is.na(log_reg_input[[i]]$indicated_investigated), 0, 1)
    
    # add drug candidates information
    log_reg_input[[i]]$recommended = ifelse(log_reg_input[[i]]$disease %in% complex_diseases_temp, 1, 0)
  }
  
  ## track progress
  cat(i, "\n")
} ; rm(i, candidate_drugs_temp, complex_diseases_temp, ind_inv_temp)

## run logistic regression
log_reg_results = vector("list", length = length(log_reg_input))
names(log_reg_results) = names(log_reg_input)

for (i in 1:length(log_reg_input)) {
  
  # run logistic regression
  glm_fits_temp = glm(indicated_investigated ~  disease_category + total_targets + recommended, 
                      data = log_reg_input[[i]], 
                      family = binomial())
  log_summary_temp = summary(glm_fits_temp)$coefficients
  
  # populate list
  log_reg_results[[i]] = log_summary_temp
    
  ## track progress
  cat(i, "\n")
} ; rm(i, glm_fits_temp, log_summary_temp)

## create a data frame with the results of logistic regression for each Mendelian disease
log_reg_results_summary = data.frame(mendelian_disease = names(md_comorbidities), beta = NA, pvalue = NA)

for (i in 1:length(log_reg_results)){
  
  log_reg_results[[i]] = data.frame(log_reg_results[[i]])
  log_reg_results_summary[i, "beta"] = log_reg_results[[i]]["recommended", "Estimate"]
  log_reg_results_summary[i, "pvalue"] = log_reg_results[[i]]["recommended", "Pr...z.."]
  
} ; rm(i)

log_reg_results_summary = log_reg_results_summary %>% 
  arrange(pvalue)
log_reg_results_summary %>% filter(pvalue < 0.05) %>% nrow() # 9 significant Mendelian diseases at a nominal level (p<0.05)

## -- annotate each Mendelian disease with number of complex disease comorbidities and number of drugs targeting its associated genes -- ##

## number of comorbidities per Mendelian disease
md_nr_comorbidities = md_cd_comorbidities %>%
  dplyr::select(mendelian_disease, complex_disease) %>%
  group_by(mendelian_disease) %>%
  mutate(nr_comorbidities = length(complex_disease)) %>%
  ungroup() %>%
  dplyr::select(mendelian_disease, nr_comorbidities) %>%
  distinct()
log_reg_results_summary = left_join(log_reg_results_summary, md_nr_comorbidities, by = "mendelian_disease")

## number of drugs targeting the genes associated with each Mendelian disease
md_genes_drugs = left_join(md_genes, db_drug_targets, by = c("causal_gene" = "drug_target"))
md_genes_drugs = na.omit(md_genes_drugs) ; rownames(md_genes_drugs) = NULL
md_genes_drugs = md_genes_drugs %>%
  group_by(mendelian_disease) %>% 
  mutate(nr_drugs = length(db_id)) %>%
  ungroup() %>%
  dplyr::select(mendelian_disease, nr_drugs) %>% 
  distinct()
log_reg_results_summary = left_join(log_reg_results_summary, md_genes_drugs, by = "mendelian_disease")

## -- visualizations -- ##
log_reg_results_summary = log_reg_results_summary %>% arrange(nr_comorbidities)
log_reg_results_summary$nr_comorbidities = factor(log_reg_results_summary$nr_comorbidities, levels = log_reg_results_summary$nr_comorbidities, labels = log_reg_results_summary$nr_comorbidities)
log_reg_results_summary$beta = ifelse(exp(log_reg_results_summary$beta) > 1, 1, 0)
log_reg_results_summary$beta = factor(log_reg_results_summary$beta, levels = c(0, 1), labels = c("No", "Yes"))

fig_s3 = ggplot(log_reg_results_summary, aes(x = nr_drugs, y = -log10(pvalue))) +
  geom_point(alpha = 0.4, size = 1.5, aes(color = beta)) +
  # scale_x_discrete(breaks = seq(1, 61, 1)) +
  scale_y_continuous(breaks = seq(0, 11, 1)) +
  labs(color = "odds ratio > 1") +
  xlab("Number of drugs") +
  ylab("neg_log10_pvalue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  theme_classic() +
  theme(axis.text = element_text(size = 14, family = "Arial", color = "black"),
        axis.title = element_text(size = 18, family = "Arial", color = "black"),
        title = element_text(size = 18, family = "Arial", color = "black"), 
        legend.text = element_text(size = 12, family = "Arial", color = "black"))

fig_s3
ggsave(filename = "FigS3_per_MD_analysis_results.tiff", 
       path = "figures/",
       width = 7, height = 5, device = 'tiff',
       dpi = 700, compression = "lzw", type = type_compression)
dev.off()

# histogram drug-to-gene ratio per Mendelian disease
nr_genes_per_md = md_genes %>%
  group_by(mendelian_disease) %>%
  mutate(nr_genes = length(causal_gene)) %>%
  ungroup() %>%
  dplyr::select(mendelian_disease, nr_genes) %>%
  distinct()
log_reg_results_summary = left_join(log_reg_results_summary, nr_genes_per_md, by = "mendelian_disease")
log_reg_results_summary$drug_to_gene_ratio = log_reg_results_summary$nr_drugs / log_reg_results_summary$nr_genes

fig3a = ggplot(log_reg_results_summary, aes(x = nr_drugs)) +
  geom_histogram(color = "black", fill = "lightblue", bins = 71) +
  scale_x_continuous(breaks = c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150)) +
  scale_y_continuous(breaks = seq(0, 9)) +
  xlab("Number of drugs") +
  ylab("Number of Mendelian diseases") +
  theme_classic() +
  theme(axis.text = element_text(size = 24, family = "Arial", color = "black"),
        axis.title = element_text(size = 30, family = "Arial", color = "black"))

fig3a
ggsave(filename = "Fig3A_drugs_per_MD.tiff", 
       path = "figures/",
       width = 12, height = 8, device = 'tiff',
       dpi = 700, compression = "lzw", type = type_compression)
dev.off()

# number of drugs - rank sum test
log_reg_results_summary$sig = ifelse(log_reg_results_summary$pvalue < 0.05, 1, 0)
log_reg_results_summary$sig = factor(log_reg_results_summary$sig, levels = c(0, 1), labels = c("Non-significant \nMendelian diseases", "Significant \nMendelian diseases"))
fig_3b = ggplot(log_reg_results_summary, aes(y = nr_drugs, x = sig)) +
  geom_boxplot() +
  # geom_jitter(width = 0.1, alpha = 0.5, size = 4) +
  xlab("") +
  ylab("Number of drugs") +
  labs(title = "") +
  scale_y_continuous(breaks = seq(0, 170, 20)) +
  theme_classic() +
  geom_signif(comparisons = list(c("Significant \nMendelian diseases", "Non-significant \nMendelian diseases")), test = "wilcox.test", test.args = list(alternative = "greater"),   
              map_signif_level = FALSE, textsize = 10) +
  annotate(geom = "text", x = 1.5, y = 150, label = "WIilcoxon rank-sum test", size = 9) +
  theme(axis.text = element_text(size = 30, family = "Arial", color = "black"),
        axis.title = element_text(size = 30),
        legend.position = "none")

fig_3b
ggsave(filename = "Fig3B_per_MD_drugs.tiff", 
       path = "figures/",
       width = 12, height = 12, device = 'tiff',
       dpi = 700, compression = "lzw", type = type_compression)
dev.off()

## -- observed vs permutations -- ##

## Now, that we have the signal that significant Mendelian diseases are targeted by a higher number of drugs compared to non-significant,
## we want to exclude the possibility that our significant results from the main analysis are due to the high number of drugs rather than the comorbidity information

# Mendelian - complex disease comorbidities
md_cd_comorbidities = md_cd_comorbidities %>% 
  dplyr::select(mendelian_disease, complex_disease) %>%
  arrange(mendelian_disease)

# drugs targeting md_genes
drugs_targeting_md_genes = left_join(md_genes, db_drug_targets, by = c("causal_gene" = "drug_target")) %>%
  na.omit() %>%
  dplyr::select(db_id) %>%
  distinct()

### -- observed -- ###

## create logistic regression input
# each row is a drug-complex disease pair
log_input_obs = data.frame(db_id = rep(drugs_targeting_md_genes$db_id, each = length(unique_cd)))

# add complex diseases
log_input_obs$complex_disease = unique_cd

# add disease category
log_input_obs = left_join(log_input_obs, complex_disease_categories, by = "complex_disease")

# add candidate drugs information
x = left_join(md_cd_comorbidities, md_genes, by = "mendelian_disease")
x = na.omit(x) ; rownames(x) = NULL
x = x %>% dplyr::select(complex_disease, causal_gene) %>% distinct
x = left_join(x, db_drug_targets, by = c("causal_gene" = "drug_target"))
x = na.omit(x) ; rownames(x) = NULL
x = x %>% dplyr::select(-causal_gene) %>% distinct()
x$recommended = 1
log_input_obs = left_join(log_input_obs, x, by = c("complex_disease", "db_id"))
log_input_obs$recommended = ifelse(is.na(log_input_obs$recommended), 0, 1)

# number of drug targets
log_input_obs = left_join(log_input_obs, drugs_nr_targets, by = c("db_id" = "db_id"))

# add investigated/indicated drugs
investigated_indicated_drugs$indicated_investigated = 1
log_input_obs = left_join(log_input_obs, investigated_indicated_drugs, by = c("db_id" = "drugbank_id", "complex_disease" = "complex_disease"))
log_input_obs$indicated_investigated = ifelse(is.na(log_input_obs$indicated_investigated), 0, 1)

# logistic regression
glm_fit_obs = glm(indicated_investigated ~ total_targets + disease_category + recommended,
                  data = log_input_obs, 
                  family = binomial())
log_reg_summary_obs = summary(glm_fit_obs)$coefficients
or_obs = exp(log_reg_summary_obs[8, 1])
pvalue_obs = log_reg_summary_obs[8, 4]

## -- permutations -- ##

## NOTE: in this permutation analysis, we shuffle the drugs targeting each Mendelian disease.
## This is equivalent to shuffling the rows of drugs in our logistic regression input

## create logistic regression input
# each row is a drug-complex disease pair
log_input_perm = data.frame(db_id = rep(drugs_targeting_md_genes$db_id, each = length(unique_cd)))

# add complex diseases
log_input_perm$complex_disease = unique_cd

# add disease category
log_input_perm = left_join(log_input_perm, complex_disease_categories, by = "complex_disease")

# add candidate drugs information
x = left_join(md_cd_comorbidities, md_genes, by = "mendelian_disease")
x = na.omit(x) ; rownames(x) = NULL
x = x %>% dplyr::select(complex_disease, causal_gene) %>% distinct
x = left_join(x, db_drug_targets, by = c("causal_gene" = "drug_target"))
x = na.omit(x) ; rownames(x) = NULL
x = x %>% dplyr::select(-causal_gene) %>% distinct()
x$recommended = 1
log_input_perm = left_join(log_input_perm, x, by = c("complex_disease", "db_id"))
log_input_perm$recommended = ifelse(is.na(log_input_perm$recommended), 0, 1)

## perform permutation analysis 
# empty data frame to populate with permutation results
log_reg_results_permutation = data.frame(permutation = NA,
                                         odds_ratio = NA,
                                         pvalue = NA)

for (permutation in 1:1000) {
  
  # shuffle the drugs
  drugs_original = unique(log_input_perm$db_id)
  drugs_shuffled = sample(drugs_original, replace = FALSE)
  log_input_perm_temp = log_input_perm
  log_input_perm_temp$db_id = rep(drugs_shuffled, each = 65)
  
  # add number of drug targets
  log_input_perm_temp = left_join(log_input_perm_temp, drugs_nr_targets, by = c("db_id" = "db_id"))
  
  # add investigated/indicated drugs
  log_input_perm_temp = left_join(log_input_perm_temp, investigated_indicated_drugs, by = c("db_id" = "drugbank_id", "complex_disease" = "complex_disease"))
  log_input_perm_temp$indicated_investigated = ifelse(is.na(log_input_perm_temp$indicated_investigated), 0, 1)
  
  # run logistic regression
  glm_fit_perm_temp = glm(indicated_investigated ~ total_targets + disease_category + recommended,
                          data = log_input_perm_temp, 
                          family = binomial())
  log_reg_summary_temp = summary(glm_fit_perm_temp)$coefficients
  
  # populate the results data frame
  log_reg_results_permutation[permutation, "permutation"] = permutation
  log_reg_results_permutation[permutation, "odds_ratio"] = exp(log_reg_summary_temp[8, 1])
  log_reg_results_permutation[permutation, "pvalue"] = log_reg_summary_temp[8, 4]

  ## track progress
  cat(permutation, "\n")
}

## calculate p-value after permutations
perm_pvalue = sum(pvalue_obs >= log_reg_results_permutation$pvalue) / 1000 # 0.018 = 1.8%
sum(or_obs <= log_reg_results_permutation$odds_ratio) / 1000 # 0.014 = 1.4%

# ## visualizations of results
# ggplot(log_reg_results_permutation, aes(x = odds_ratio)) +
#   geom_histogram(color = "black", fill = "lightblue") +
#   geom_vline(xintercept = or_obs, color = "red", linewidth = 0.8) +
#   xlab("Odds ratio") +
#   ylab("count") +
#   scale_x_continuous(breaks = seq(0, 2.25, 0.25), limits = c(0.75, 2.25)) +
#   annotate(geom = "text", x = 1.95, y = 145, label = paste0(sum(or_obs <= log_reg_results_permutation$odds_ratio) / 1000*100, "%"), size = 7, family = "Arial", color = "red", fontface = "bold") +
#   theme_classic() +
#   theme(axis.text = element_text(size = 20, family = "Arial", colour = "black"),
#         axis.title = element_text(size = 22, family = "Arial", colour = "black"))

## -- gene level analysis -- ##

# histogram of number of drugs per Mendelian gene
md_genes_drugs = left_join(md_genes, db_drug_targets, by = c("causal_gene" = "drug_target"))
md_genes_drugs = md_genes_drugs %>%
  dplyr::select(causal_gene, db_id) %>%
  distinct() %>%
  group_by(causal_gene) %>% 
  mutate(nr_drugs = if_else(sum(is.na(unique(db_id))) == 1, "0", "over_0")) %>%
  mutate(nr_drugs = if_else(nr_drugs != "0", length(unique(db_id)), 0)) %>%
  ungroup() %>%
  dplyr::select(causal_gene, nr_drugs) %>% 
  distinct() %>%
  filter(nr_drugs > 0)

fi3c = ggplot(md_genes_drugs, aes(x = nr_drugs)) +
  geom_histogram(color = "black", fill = "lightblue", bins = 75) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85)) +
  scale_y_continuous(breaks = seq(0, 75, 5)) +
  xlab("Number of drugs") +
  ylab("Number of genes") +
  theme_classic() +
  theme(axis.text = element_text(size = 24, family = "Arial", color = "black"),
        axis.title = element_text(size = 30, family = "Arial", color = "black"))

fi3c
ggsave(filename = "Fig3C_drugs_per_MD_gene.tiff", 
       path = "figures/",
       width = 12, height = 8, device = 'tiff',
       dpi = 700, compression = "lzw", type = type_compression)
dev.off()

## analysis per Mendelian disease gene
# keep druggable md genes
md_genes_druggable = md_genes %>% filter(causal_gene %in% db_drug_targets$drug_target)
unique_md_genes = unique(md_genes_druggable$causal_gene)
lr_inputs = vector("list", length(unique(unique_md_genes)))
names(lr_inputs) = unique_md_genes

for (i in 1:length(unique_md_genes)) {
  
  # gene
  md_gene_temp = unique_md_genes[[i]]
  
  # drugs targeting that gene
  drugs_targeting_md_gene_temp = db_drug_targets %>% filter(drug_target == md_gene_temp)
  drugs_targeting_md_gene_temp = unique(drugs_targeting_md_gene_temp$db_id)
  
  # Mendelian diseases linked to that gene
  md_disease_temp = md_genes %>% filter(causal_gene == md_gene_temp)
  md_disease_temp = unique(md_disease_temp$mendelian_disease)
  
  # comorbidities of these Mendelian diseases
  md_comorbidities_temp = md_cd_comorbidities %>% filter(mendelian_disease %in% md_disease_temp)
  md_comorbidities_temp = unique(md_comorbidities_temp$complex_disease)
  
  # each row should be a drug-complex disease pair
  log_input = data.frame(drug = rep(drugs_targeting_md_gene_temp, each = length(unique_cd)))
  
  # add complex diseases
  log_input$disease = unique_cd
  
  # add disease category
  log_input = left_join(log_input, complex_disease_categories, by = c("disease" = "complex_disease"))
  
  # add number of targets
  log_input = left_join(log_input, drugs_nr_targets, by = c("drug" = "db_id"))
  
  # add indicated/investigated drugs
  log_input = left_join(log_input, investigated_indicated_drugs, by = c("drug" = "drugbank_id", "disease" = "complex_disease"))
  log_input$indicated_investigated = ifelse(is.na(log_input$indicated_investigated), 0, 1)
  
  # add recommended candidate drugs
  log_input$recommended = ifelse(log_input$disease %in% md_comorbidities_temp, 1, 0)
  
  lr_inputs[[i]] = log_input
  
  cat(i, "-", length(unique_md_genes), "\n")
}

## run logistic regression
log_reg_results = data.frame(md_genes = unique_md_genes, comorbidity_OR = NA, comorbidity_P = NA)

for (i in 1:nrow(log_reg_results)) {
  
  # run logistic regression
  if (length(unique(lr_inputs[[i]]$total_targets)) == 1) {
    glm_fits_temp = glm(indicated_investigated ~  disease_category + recommended, 
                        data = lr_inputs[[i]], 
                        family = binomial())
    log_summary_temp = as.data.frame(summary(glm_fits_temp)$coefficients)
    log_reg_results[i, "comorbidity_OR"] = exp(log_summary_temp[7, 1])
    log_reg_results[i, "comorbidity_P"] = log_summary_temp[7, 4]
    
  }
  
  if (length(unique(lr_inputs[[i]]$total_targets)) > 1) {
    glm_fits_temp = glm(indicated_investigated ~  disease_category + total_targets + recommended, 
                        data = lr_inputs[[i]], 
                        family = binomial())
    log_summary_temp = as.data.frame(summary(glm_fits_temp)$coefficients)
    log_reg_results[i, "comorbidity_OR"] = exp(log_summary_temp[8, 1])
    log_reg_results[i, "comorbidity_P"] = log_summary_temp[8, 4]
  }
  
  cat(i, "-", nrow(log_reg_results), "\n")
}

log_reg_results$sig = ifelse(log_reg_results$comorbidity_P < 0.05, 1, 0)

## add nr of targets for each MD gene
md_genes_druggable = left_join(md_genes_druggable, db_drug_targets, by = c("causal_gene" = "drug_target"))
md_genes_druggable = md_genes_druggable %>% 
  group_by(causal_gene) %>%
  mutate(nr_drugs = length(unique(db_id))) %>%
  ungroup() %>%
  distinct()
md_genes_nr_drugs = md_genes_druggable %>%
  dplyr::select(causal_gene, nr_drugs) %>% 
  distinct()

log_reg_results = left_join(log_reg_results, md_genes_nr_drugs, by = c("md_genes" = "causal_gene"))
log_reg_results$sig = factor(log_reg_results$sig, levels = c("0", "1"), labels = c("Non-significant \ngenes", "Significant \ngenes"))

fig_3d = ggplot(log_reg_results, aes(x = sig, y = nr_drugs)) +
  geom_boxplot() +
  # geom_point(alpha = 0.5, position = "jitter") +
  xlab("") +
  ylab("Number of drugs") +
  geom_signif(comparisons = list(c("Significant \ngenes", "Non-significant \ngenes")), test.args = list(alternative = "greater"), 
              map_signif_level = FALSE, textsize = 10) +
  annotate(geom = "text", x = 1.5, y = 85, label = "WIilcoxon rank-sum test", size = 9) +
  scale_y_continuous(breaks = seq(0, 85, 10)) +
  theme_classic() +
  theme(axis.text = element_text(size = 30, family = "Arial", color = "black"),
        axis.title = element_text(size = 30),
        legend.position = "none")

fig_3d
ggsave(filename = "Fig3D_per_gene_drugs.tiff", 
       path = "figures/",
       width = 12, height = 12, device = 'tiff',
       dpi = 700, compression = "lzw", type = type_compression)
dev.off()
