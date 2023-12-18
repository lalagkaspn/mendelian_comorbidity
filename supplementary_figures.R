## supplemenetary figures of observed and permuted ORs for the main analysis ##

library(data.table)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)

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

## -- main analysis -- ##

## all comorbidities - observed

# load comorbidities
md_cd_comorbidities = fread("processed_data/md_cd_comorbidities.txt") %>% 
  dplyr::select(mendelian_disease, complex_disease) %>%
  arrange(mendelian_disease)

nr_comorbidities_per_md = md_cd_comorbidities %>% 
  group_by(mendelian_disease) %>%
  mutate(nr_comorbidities = length(complex_disease)) %>%
  ungroup() %>%
  dplyr::select(-complex_disease) %>% 
  distinct()

# comorbidity matrix
md_cd_comorbidities_matrix = as.data.frame.matrix(table(md_cd_comorbidities[, 2:1]))
md_cd_comorbidities_matrix = md_cd_comorbidities_matrix %>% mutate(complex_disease = rownames(md_cd_comorbidities_matrix), .before = Achromatopsia)
rownames(md_cd_comorbidities_matrix) = NULL

# load md genes
md_genes = fread("processed_data/md_genes.txt")

# load drug - targets
db_drug_targets = fread("processed_data/drugbank_all_drugs_known_targets.txt") %>%
  dplyr::select(db_id, drug_target) %>% 
  distinct()

# calculate number of targets per drug
drugs_nr_targets = fread("processed_data/drugbank_all_drugs_known_targets.txt") %>%
  dplyr::select(db_id, drug_target) %>%
  group_by(db_id) %>% 
  mutate(total_targets = length(drug_target)) %>%
  dplyr::select(db_id, total_targets) %>% 
  distinct() %>%
  ungroup()

# load complex disease categories
complex_disease_categories = fread("raw_data/complex_disease_category.txt")

# load indicated/investigated drugs
investigated_indicated_drugs = fread("processed_data/drugs_inv_ind_per_disease.txt")
investigated_indicated_drugs = melt(investigated_indicated_drugs, "drugbank_id", colnames(investigated_indicated_drugs)[2:ncol(investigated_indicated_drugs)])
investigated_indicated_drugs = na.omit(investigated_indicated_drugs)
rownames(investigated_indicated_drugs) = NULL
investigated_indicated_drugs = investigated_indicated_drugs %>% dplyr::select(drugbank_id, complex_disease = variable) %>% distinct()
investigated_indicated_drugs$indicated_investigated = 1

# run logistic regression - observed ##
log_reg_observed_nrcomorbidities = data.frame(comorbidity_filter = 9:61, 
                                              observed_OR = 0,
                                              observed_P = 0,
                                              observed_OR_95ci_lower = 0,
                                              observed_OR_95ci_upper = 0,
                                              logr_cd_inlcuded = 0,
                                              logr_drugs_included = 0)
for (i in 1:nrow(log_reg_observed_nrcomorbidities)) {
  
  comorbidity_filter = log_reg_observed_nrcomorbidities[i, "comorbidity_filter"]
  mds = nr_comorbidities_per_md %>% filter(nr_comorbidities <= comorbidity_filter)
  
  md_cd_comorbidities_filt = md_cd_comorbidities %>% filter(mendelian_disease %in% mds$mendelian_disease)
  # unique Mendelian and complex diseases
  unique_cd_temp = unique(md_cd_comorbidities_filt$complex_disease)
  unique_md_temp = unique(md_cd_comorbidities_filt$mendelian_disease)
  
  md_genes_filt = md_genes %>% filter(mendelian_disease %in% mds$mendelian_disease)
  
  drugs_targeting_md_genes = left_join(md_genes_filt, db_drug_targets, by = c("causal_gene" = "drug_target")) %>%
    na.omit() %>%
    dplyr::select(db_id) %>%
    distinct()
  # create logistic regression input - each row is a drug-complex disease pair 
  log_input_all = data.frame(db_id = rep(drugs_targeting_md_genes$db_id, each = length(unique_cd_temp)))
  log_input_all$complex_disease = unique_cd_temp
  # add information about disease category
  log_input_all = left_join(log_input_all, complex_disease_categories, by = "complex_disease")
  # add information about number of targets for each drug
  log_input_all = left_join(log_input_all, drugs_nr_targets, by = "db_id")  
  # add information about investigated/indicated drugs within each drug-cancer pair
  log_input_all = left_join(log_input_all, investigated_indicated_drugs, by = c("complex_disease" = "complex_disease", "db_id" = "drugbank_id"))
  log_input_all$indicated_investigated = ifelse(is.na(log_input_all$indicated_investigated), 0, 1)
  # add information about comorbidity
  x = left_join(md_cd_comorbidities_filt, md_genes_filt, by = "mendelian_disease")
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
  log_reg_observed_nrcomorbidities[i, "logr_cd_inlcuded"] = length(unique(log_input_all$complex_disease))
  log_reg_observed_nrcomorbidities[i, "logr_drugs_included"] = length(unique(log_input_all$db_id))
  summary_temp = as.data.frame(summary(log_reg_all)$coefficients)
  summary_temp$rownames = rownames(summary_temp)
  com_temp = summary_temp %>% filter(rownames == "comorbidity")
  log_reg_observed_nrcomorbidities[i, "observed_OR"] = exp(com_temp$Estimate)
  log_reg_observed_nrcomorbidities[i, "observed_P"] = com_temp$`Pr(>|z|)`
  
  ci = as.data.frame(exp(confint(log_reg_all)))
  log_reg_observed_nrcomorbidities[i, "observed_OR_95ci_lower"] = ci["comorbidity", 1]
  log_reg_observed_nrcomorbidities[i, "observed_OR_95ci_upper"] = ci["comorbidity", 2]
}
rm(comorbidity_filter, mds, md_cd_comorbidities_filt, unique_cd_temp, unique_md_temp, x, md_genes_filt, drugs_targeting_md_genes, log_input_all, log_reg_all, com_temp, summary_temp, ci, i)

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
all_or_perm = vector("list", 1000)
names(all_or_perm) = paste0("permutation_", 1:1000)

for (permutation in 1:1000) {
  
  log_reg_perm_nrcomorbidities = data.frame(comorbidity_filter = 9:61, 
                                            perm_OR = 0,
                                            perm_P = 0)
  
  for (i in 1:nrow(log_reg_perm_nrcomorbidities)) {
    
    comorbidity_filter = log_reg_perm_nrcomorbidities[i, "comorbidity_filter"]
    mds = nr_comorbidities_per_md %>% filter(nr_comorbidities <= comorbidity_filter)
    
    ## comorbidity matrix
    md_cd_comorbidities_matrix_perm = as.data.frame.matrix(table(md_cd_comorbidities_filt[, 2:1]))
    md_cd_comorbidities_matrix_perm = md_cd_comorbidities_matrix_perm %>% mutate(complex_disease = rownames(md_cd_comorbidities_matrix_perm), .before = Achromatopsia)
    rownames(md_cd_comorbidities_matrix_perm) = NULL
    
    md_cd_comorbidities_filt = md_cd_comorbidities %>% filter(mendelian_disease %in% mds$mendelian_disease)
    
    ## shuffle comorbidity matrix
    md_cd_comorbidities_shuffled = shuffle(md_cd_comorbidities_matrix_perm)
    
    # unique Mendelian and complex diseases
    unique_cd_temp = unique(md_cd_comorbidities_shuffled$complex_disease)
    unique_md_temp = unique(md_cd_comorbidities_shuffled$mendelian_disease)
    
    md_genes_filt = md_genes %>% filter(mendelian_disease %in% md_cd_comorbidities_shuffled$mendelian_disease)
    
    drugs_targeting_md_genes = left_join(md_genes_filt, db_drug_targets, by = c("causal_gene" = "drug_target")) %>%
      na.omit() %>%
      dplyr::select(db_id) %>%
      distinct()
    # create logistic regression input - each row is a drug-complex disease pair 
    log_input_all = data.frame(db_id = rep(drugs_targeting_md_genes$db_id, each = length(unique_cd_temp)))
    log_input_all$complex_disease = unique_cd_temp
    # add information about disease category
    log_input_all = left_join(log_input_all, complex_disease_categories, by = "complex_disease")
    # add information about number of targets for each drug
    log_input_all = left_join(log_input_all, drugs_nr_targets, by = "db_id")  
    # add information about investigated/indicated drugs within each drug-cancer pair
    log_input_all = left_join(log_input_all, investigated_indicated_drugs, by = c("complex_disease" = "complex_disease", "db_id" = "drugbank_id"))
    log_input_all$indicated_investigated = ifelse(is.na(log_input_all$indicated_investigated), 0, 1)
    # add information about comorbidity
    x = left_join(md_cd_comorbidities_shuffled, md_genes_filt, by = "mendelian_disease")
    x = na.omit(x) ; rownames(x) = NULL
    x = x %>% dplyr::select(complex_disease, causal_gene) %>% distinct
    x = left_join(x, db_drug_targets, by = c("causal_gene" = "drug_target"))
    x = na.omit(x) ; rownames(x) = NULL
    x = x %>% dplyr::select(-causal_gene) %>% distinct()
    x$comorbidity = 1
    log_input_all = left_join(log_input_all, x, by = c("complex_disease", "db_id"))
    log_input_all$comorbidity = ifelse(is.na(log_input_all$comorbidity), 0, 1)
    # logistic regression
    if (length(unique(log_input_all$disease_category)) == 1) {
      log_reg_all = glm(indicated_investigated ~ total_targets + comorbidity,
                        data = log_input_all, 
                        family = binomial())
    } else {
      log_reg_all = glm(indicated_investigated ~ total_targets + disease_category + comorbidity,
                        data = log_input_all, 
                        family = binomial())
    }
    
    
    summary_temp = as.data.frame(summary(log_reg_all)$coefficients)
    summary_temp$rownames = rownames(summary_temp)
    com_temp = summary_temp %>% filter(rownames == "comorbidity")
    log_reg_perm_nrcomorbidities[i, "perm_OR"] = exp(com_temp$Estimate)
    log_reg_perm_nrcomorbidities[i, "perm_P"] = com_temp$`Pr(>|z|)`
  }
  
  all_or_perm[[permutation]] = log_reg_perm_nrcomorbidities
  
  ## track progress
  cat(permutation, "\n")
}

all_permuted_ors = data.frame(comorbidity_filter = NA,
                              perm_OR = NA,
                              perm_P = NA) 

for (i in 9:61) {
  
  comorbidity_filter_perm = i
  
  temp = data.frame(comorbidity_filter = rep(i, 1000),
                    perm_OR = NA,
                    perm_P = NA)
  for (z in 1:length(all_or_perm)) {
    
    perm_values_temp = all_or_perm[[z]] %>% filter(comorbidity_filter == comorbidity_filter_perm)
    temp[z, "perm_OR"] = perm_values_temp$perm_OR
    temp[z, "perm_P"] = perm_values_temp$perm_P
  }
  
  all_permuted_ors = rbind(all_permuted_ors, temp)
  cat(i, "\n")
} ; rm(temp, i, z, comorbidity_filter_perm, perm_values_temp)
all_permuted_ors = all_permuted_ors[-1, ] ; rownames(all_permuted_ors) = NULL

all_permuted_ors$comorbidity_filter = paste0("≤", all_permuted_ors$comorbidity_filter)
all_permuted_ors$comorbidity_filter = factor(all_permuted_ors$comorbidity_filter, levels = all_permuted_ors$comorbidity_filter, labels = all_permuted_ors$comorbidity_filter)

## For observed OR: plot OR with 95% CI
log_reg_observed_nrcomorbidities$comorbidity_filter = paste0("≤", log_reg_observed_nrcomorbidities$comorbidity_filter)
log_reg_observed_nrcomorbidities$comorbidity_filter = factor(log_reg_observed_nrcomorbidities$comorbidity_filter, levels = log_reg_observed_nrcomorbidities$comorbidity_filter, labels = log_reg_observed_nrcomorbidities$comorbidity_filter)
log_reg_observed_nrcomorbidities$status = "observed"

## For permuted ORs: plot 5th, 50th and 95th percentiles
all_permuted_ors_filt = all_permuted_ors %>% 
  group_by(comorbidity_filter) %>%
  mutate(percentile_5 = as.numeric(quantile(perm_OR, 0.05)),
         percentile_50 = as.numeric(quantile(perm_OR, 0.50)),
         percentile_95 = as.numeric(quantile(perm_OR, 0.95))) %>%
  filter(perm_OR >= percentile_5 & perm_OR <= percentile_95) %>%
  ungroup() %>%
  dplyr::select(comorbidity_filter, percentile_5, percentile_50, percentile_95) %>%
  distinct() %>%
  mutate(status = "permuted")
# rename columns for rbind later --> will change them after that
colnames(all_permuted_ors_filt) = c("comorbidity_filter", "observed_OR_95ci_lower", "observed_OR", "observed_OR_95ci_upper", "status")

all_obs_perm_combined = rbind(log_reg_observed_nrcomorbidities[, c(1,2,4,5,8)], all_permuted_ors_filt)
colnames(all_obs_perm_combined) = c("comorbidity_filter", "OR", "OR_lower_limit", "OR_upper_limit","status")
all_obs_perm_combined$comorbidity_filter = rep(9:61, 2)
all_obs_perm_combined$status = factor(all_obs_perm_combined$status, levels = all_obs_perm_combined$status, labels = all_obs_perm_combined$status)

s1_a = ggplot(all_obs_perm_combined) +
  geom_ribbon(aes(x = comorbidity_filter, ymin = OR_lower_limit, ymax = OR_upper_limit, fill = status),
              alpha = 0.3) +
  scale_fill_manual(values = c("#00BFC4","#F8766D"))+
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 2) +
  geom_line(aes(x = comorbidity_filter, y = OR, color = status), linewidth = 2, show.legend = FALSE) +
  scale_color_manual(values = c("#00BFC4","#F8766D")) +
  xlab("") +
  scale_x_continuous(breaks = seq(9, 61, 4), labels = paste0("≤", seq(9, 61, 4)), expand = c(0, 0.5)) +
  ylab("Odds Ratio") +
  labs(fill = "") +
  theme_classic() +
  theme(axis.title = element_text(size = 45, color = "black", family = "Arial"),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text.x = element_text(size = 40, color = "black", family = "Arial", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "Arial", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_text(size = 23, color = "black", family = "Arial"),
        legend.text = element_text(size = 40, family = "Arial"),
        legend.position = "top", legend.key.size = unit(2, "cm"))
s1_b = ggplot(log_reg_observed_nrcomorbidities, aes(x = comorbidity_filter, y = logr_drugs_included)) +
  geom_line(group = "comorbidity_filter", linewidth = 1.5) +
  scale_x_discrete(breaks = paste0("≤", seq(9, 61, 4)), labels = paste0("≤", seq(9, 61, 4)), expand = c(0, 0.5)) +
  scale_y_continuous(breaks = seq(0, 800, 200)) +
  xlab("Comorbidities per Mendelian disease") +
  ylab("Number of drugs") +
  labs(fill = "") +
  theme_classic() +
  theme(axis.title = element_text(size = 45, color = "black", family = "Arial"),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text.x = element_text(size = 40, color = "black", family = "Arial", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "Arial", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"))
ggarrange(s1_a, s1_b, 
          ncol = 1, nrow = 2, heights = c(3,1), align = "v")
ggsave(filename = "S1_main_analysis.png",
       path = "supplementary_figures/", 
       width = 45, height = 30, device = "png",
       dpi = 300, type = type_compression)

## -- Per disease category -- ##

# unique complex diseases
cancers = complex_disease_categories %>% filter(disease_category == "cancer") ; cancers = unique(cancers$complex_disease)
cardiovascular = complex_disease_categories %>% filter(disease_category == "cardiovascular") ; cardiovascular = unique(cardiovascular$complex_disease)
hormonal = complex_disease_categories %>% filter(disease_category == "hormonal") ; hormonal = unique(hormonal$complex_disease)
immune = complex_disease_categories %>% filter(disease_category == "immune") ; immune = unique(immune$complex_disease)
neurological = complex_disease_categories %>% filter(disease_category == "neurological") ; neurological = unique(neurological$complex_disease)
ophthalmological = complex_disease_categories %>% filter(disease_category == "ophthalmological") ; ophthalmological = unique(ophthalmological$complex_disease)

## observed OR ##
# define function
per_dis_cat_obs_or = function(cd_in_category, comorbidity_dataframe, drug_targets, nr_drug_targets) {
  
  ## filter for MD comorbidities of these complex diseases
  md_cd_comorbidities_in_category = comorbidity_dataframe %>% filter(complex_disease %in% cd_in_category)
  ## find number of comorbidities per Mendelian disease
  md_nr_comorbidities_temp = md_cd_comorbidities_in_category %>% 
    group_by(mendelian_disease) %>%
    mutate(nr_comorbidities = length(complex_disease)) %>%
    ungroup() %>%
    dplyr::select(mendelian_disease, nr_comorbidities) %>%
    distinct()
  min_comorbidities_temp = summary(md_nr_comorbidities_temp$nr_comorbidities)[[1]]
  max_comorbidities_temp = summary(md_nr_comorbidities_temp$nr_comorbidities)[[6]]
  
  # create empty data frame to populate later
  log_reg_results_temp = data.frame(comorbidity_filter = min_comorbidities_temp:max_comorbidities_temp,
                                    obs_or = 0,
                                    obs_pvalue = 0,
                                    ci_95_lower = 0,
                                    ci_95_upper = 0,
                                    nr_drugs_in_lr = 0)
  for (i in min_comorbidities_temp:max_comorbidities_temp) {
    
    # filter comorbidity matrix
    md_to_include = md_nr_comorbidities_temp %>% filter(nr_comorbidities <= i)
    # MD-CD comorbidities
    md_cd_comorbidities_filt = md_cd_comorbidities_in_category %>% filter(mendelian_disease %in% md_to_include$mendelian_disease)
    # unique md and cd based on the comorbidity filter applied
    unique_cds_temp = unique(md_cd_comorbidities_filt$complex_disease)
    unique_mds_temp = unique(md_cd_comorbidities_filt$mendelian_disease)
    # MD genes
    md_genes_temp = md_genes %>% filter(mendelian_disease %in% unique_mds_temp)
    # drugs targeting md genes
    drugs_targeting_md_genes_temp = left_join(md_genes_temp, db_drug_targets, by = c("causal_gene" = "drug_target")) %>%
      na.omit() %>%
      dplyr::select(db_id) %>%
      distinct() 
    cd_md_genes_drugs = left_join(md_genes_temp, db_drug_targets, by = c("causal_gene" = "drug_target")) %>% na.omit()
    if (nrow(cd_md_genes_drugs) == 0) {
      log_reg_results_temp[i, "obs_or"] = NA
      log_reg_results_temp[i, "obs_pvalue"] = NA
      
      log_reg_results_temp[i, "ci_95_lower"] = NA
      log_reg_results_temp[i, "ci_95_upper"] = NA
      
      log_reg_results_temp[i, "nr_drugs_in_lr"] = NA
    }
    if (nrow(cd_md_genes_drugs) > 0) {
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
      x = left_join(md_cd_comorbidities_filt, cd_md_genes_drugs, by = "mendelian_disease")
      x = na.omit(x) ; rownames(x) = NULL
      x = x %>% dplyr::select(complex_disease, db_id) %>% distinct
      x$comorbidity = 1
      log_input = left_join(log_input, x, by = c("complex_disease", "db_id"))
      log_input$comorbidity = ifelse(is.na(log_input$comorbidity), 0, 1)
      
      ## run logistic regression ##
      glm_fits = glm(indicated_investigated ~ total_targets + comorbidity,
                     data = log_input, 
                     family = binomial())
      log_reg_summary_temp = summary(glm_fits)$coefficients
      
      if (nrow(log_reg_summary_temp) >= 2){
        log_reg_results_temp[i, "obs_or"] = exp(log_reg_summary_temp["comorbidity", "Estimate"])
        log_reg_results_temp[i, "obs_pvalue"] = log_reg_summary_temp["comorbidity", "Pr(>|z|)"]
        
        if (length(unique(log_input$indicated_investigated)) == 1) {
          log_reg_results_temp[i, "ci_95_lower"] = NA
          log_reg_results_temp[i, "ci_95_upper"] = NA
          log_reg_results_temp[i, "nr_drugs_in_lr"] = length(unique(log_input$db_id))
        } else {
          log_reg_results_temp[i, "ci_95_lower"] = exp(confint(glm_fits)["comorbidity", "2.5 %"])
          log_reg_results_temp[i, "ci_95_upper"] = exp(confint(glm_fits)["comorbidity", "97.5 %"])
          log_reg_results_temp[i, "nr_drugs_in_lr"] = length(unique(log_input$db_id))
        }
      }
      if (nrow(log_reg_summary_temp) < 2){
        log_reg_results_temp[i, "obs_or"] = NA
        log_reg_results_temp[i, "obs_pvalue"] = NA
        
        log_reg_results_temp[i, "ci_95_lower"] = NA
        log_reg_results_temp[i, "ci_95_upper"] = NA
        
        log_reg_results_temp[i, "nr_drugs_in_lr"] = length(unique(log_input$db_id))
      }
      
      
    }
  }
  return(log_reg_results_temp)
}
rm(comorbidity_dataframe, cd_in_category, md_cd_comorbidities_in_category, md_nr_comorbidities_temp, min_comorbidities_temp, max_comorbidities_temp,
   i, md_to_include, md_cd_comorbidities_filt, unique_cds_temp, unique_mds_temp, md_genes_temp, drugs_targeting_md_genes_temp,
   cd_md_genes_drugs, log_input, x, glm_fits, log_reg_summary_temp, log_reg_results_temp)
# run function per disease category
cancers_obs_or = per_dis_cat_obs_or(cd_in_category = cancers, comorbidity_dataframe = md_cd_comorbidities, drug_targets = db_drug_targets, nr_drug_targets = drugs_nr_targets)
hormonal_obs_or = per_dis_cat_obs_or(cd_in_category = hormonal, comorbidity_dataframe = md_cd_comorbidities, drug_targets = db_drug_targets, nr_drug_targets = drugs_nr_targets)
cardiovascular_obs_or = per_dis_cat_obs_or(cd_in_category = cardiovascular, comorbidity_dataframe = md_cd_comorbidities, drug_targets = db_drug_targets, nr_drug_targets = drugs_nr_targets)
ophthalmological_obs_or = per_dis_cat_obs_or(cd_in_category = ophthalmological, comorbidity_dataframe = md_cd_comorbidities, drug_targets = db_drug_targets, nr_drug_targets = drugs_nr_targets)
immune_obs_or = per_dis_cat_obs_or(cd_in_category = immune, comorbidity_dataframe = md_cd_comorbidities, drug_targets = db_drug_targets, nr_drug_targets = drugs_nr_targets)
neurological_obs_or = per_dis_cat_obs_or(cd_in_category = neurological, comorbidity_dataframe = md_cd_comorbidities, drug_targets = db_drug_targets, nr_drug_targets = drugs_nr_targets)

## permuted ORs ##
# define function
per_dis_cat_perm_or = function(comorbidity_dataframe, cd_in_category, drug_targets, nr_drug_targets) {
  
  ## filter for MD comorbidities of these complex diseases
  md_cd_comorbidities_in_category = comorbidity_dataframe %>% filter(complex_disease %in% cd_in_category)
  ## find number of comorbidities per Mendelian disease
  md_nr_comorbidities_temp = md_cd_comorbidities_in_category %>% 
    group_by(mendelian_disease) %>%
    mutate(nr_comorbidities = length(complex_disease)) %>%
    ungroup() %>%
    dplyr::select(mendelian_disease, nr_comorbidities) %>%
    distinct()
  min_comorbidities_temp = summary(md_nr_comorbidities_temp$nr_comorbidities)[[1]]
  max_comorbidities_temp = summary(md_nr_comorbidities_temp$nr_comorbidities)[[6]]
  
  # create empty data frame to populate later
  perm_or_temp = vector("list", length = 1000) ; names(perm_or_temp) = paste0("permutation_", 1:1000)
  
  for (permutation in 1:length(perm_or_temp)) {
    log_reg_results_temp = data.frame(comorbidity_filter = min_comorbidities_temp:max_comorbidities_temp,
                                      perm_or = 0,
                                      perm_pvalue = 0)
    
    for (i in min_comorbidities_temp:max_comorbidities_temp) {
      md_to_include = md_nr_comorbidities_temp %>% filter(nr_comorbidities <= i)
      # MD-CD comorbidities
      md_cd_comorbidities_filt = md_cd_comorbidities_in_category %>% filter(mendelian_disease %in% md_to_include$mendelian_disease)
      # shuffle the comorbidity matrix
      comorbidity_matrix_temp = as.data.frame.matrix(table(md_cd_comorbidities_filt[, 2:1]))
      comorbidity_matrix_temp = comorbidity_matrix_temp %>% mutate(complex_disease = rownames(comorbidity_matrix_temp), .before = colnames(comorbidity_matrix_temp)[1])
      rownames(comorbidity_matrix_temp) = NULL
      comorbidity_matrix_shuffled = comorbidity_matrix_temp
      for (z in 1:ncol(comorbidity_matrix_shuffled)) {
        comorbidity_matrix_shuffled[, z] = sample(comorbidity_matrix_shuffled[, z], replace = FALSE)
      }
      comorbidity_matrix_shuffled = tidyr::gather(comorbidity_matrix_shuffled, key = "mendelian_disease", value = "comorbidity", -complex_disease) %>%
        dplyr::relocate(mendelian_disease, .before = complex_disease) %>% 
        filter(comorbidity == 1) %>% 
        dplyr::select(-comorbidity) %>%
        distinct()
      # unique md and cd based on the comorbidity filter applied
      unique_cds_temp = unique(comorbidity_matrix_shuffled$complex_disease)
      unique_mds_temp = unique(comorbidity_matrix_shuffled$mendelian_disease)
      # MD genes
      md_genes_temp = md_genes %>% filter(mendelian_disease %in% unique_mds_temp)
      # drugs targeting md genes
      drugs_targeting_md_genes_temp = left_join(md_genes_temp, db_drug_targets, by = c("causal_gene" = "drug_target")) %>%
        na.omit() %>%
        dplyr::select(db_id) %>%
        distinct() 
      cd_md_genes_drugs = left_join(md_genes_temp, db_drug_targets, by = c("causal_gene" = "drug_target")) %>% na.omit()
      if (nrow(cd_md_genes_drugs) == 0) {
        log_reg_results_temp[i, "perm_or"] = NA
        log_reg_results_temp[i, "perm_pvalue"] = NA
      }
      if (nrow(cd_md_genes_drugs) > 0) {
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
        log_reg_summary_temp = summary(glm_fits)$coefficients
        if (nrow(log_reg_summary_temp) >= 2){
          log_reg_results_temp[i, "perm_or"] = exp(log_reg_summary_temp["comorbidity", "Estimate"])
          log_reg_results_temp[i, "perm_pvalue"] = log_reg_summary_temp["comorbidity", "Pr(>|z|)"]
        }
        if (nrow(log_reg_summary_temp) < 2){
          log_reg_results_temp[i, "perm_or"] = NA
          log_reg_results_temp[i, "perm_pvalue"] = NA
        }
      }
    }
    perm_or_temp[[paste0("permutation_", permutation)]] = log_reg_results_temp
    cat(permutation, "\n")
  }
  return(perm_or_temp)
}
# run function per disease category
cancers_perm_or = per_dis_cat_perm_or(comorbidity_dataframe = md_cd_comorbidities, cd_in_category = cancers, drug_targets = db_drug_targets, nr_drug_targets = drugs_nr_targets)
hormonal_perm_or = per_dis_cat_perm_or(comorbidity_dataframe = md_cd_comorbidities, cd_in_category = hormonal, drug_targets = db_drug_targets, nr_drug_targets = drugs_nr_targets)
cardiovascular_perm_or = per_dis_cat_perm_or(comorbidity_dataframe = md_cd_comorbidities, cd_in_category = cardiovascular, drug_targets = db_drug_targets, nr_drug_targets = drugs_nr_targets)
ophthalmological_perm_or = per_dis_cat_perm_or(comorbidity_dataframe = md_cd_comorbidities, cd_in_category = ophthalmological, drug_targets = db_drug_targets, nr_drug_targets = drugs_nr_targets)
immune_perm_or = per_dis_cat_perm_or(comorbidity_dataframe = md_cd_comorbidities, cd_in_category = immune, drug_targets = db_drug_targets, nr_drug_targets = drugs_nr_targets)
neurological_perm_or = per_dis_cat_perm_or(comorbidity_dataframe = md_cd_comorbidities, cd_in_category = neurological, drug_targets = db_drug_targets, nr_drug_targets = drugs_nr_targets)

# combine all permuted ORs to one data frame
cancers_perm_or_all = data.frame(comorbidity_filter = NA,
                                 perm_OR = NA,
                                 perm_P = NA) 
for (i in 1:nrow(cancers_obs_or)) {
  comorbidity_filter_perm = i
  temp = data.frame(comorbidity_filter = rep(i, 1000),
                    perm_OR = NA,
                    perm_P = NA)
  for (z in 1:length(cancers_perm_or)) {
    perm_values_temp = cancers_perm_or[[z]] %>% filter(comorbidity_filter == comorbidity_filter_perm)
    temp[z, "perm_OR"] = perm_values_temp$perm_or
    temp[z, "perm_P"] = perm_values_temp$perm_pvalue
  }
  cancers_perm_or_all = rbind(cancers_perm_or_all, temp)
  cat(i, "\n")
} ; rm(temp, i, z, comorbidity_filter_perm, perm_values_temp)
cancers_perm_or_all = cancers_perm_or_all[-1, ] ; rownames(cancers_perm_or_all) = NULL
cancers_perm_or_all$comorbidity_filter = paste0("≤", cancers_perm_or_all$comorbidity_filter)
cancers_perm_or_all$comorbidity_filter = factor(cancers_perm_or_all$comorbidity_filter, levels = cancers_perm_or_all$comorbidity_filter, labels = cancers_perm_or_all$comorbidity_filter)

hormonal_perm_or_all = data.frame(comorbidity_filter = NA,
                                 perm_OR = NA,
                                 perm_P = NA) 
for (i in 1:nrow(hormonal_obs_or)) {
  comorbidity_filter_perm = i
  temp = data.frame(comorbidity_filter = rep(i, 1000),
                    perm_OR = NA,
                    perm_P = NA)
  for (z in 1:length(hormonal_perm_or)) {
    perm_values_temp = hormonal_perm_or[[z]] %>% filter(comorbidity_filter == comorbidity_filter_perm)
    temp[z, "perm_OR"] = perm_values_temp$perm_or
    temp[z, "perm_P"] = perm_values_temp$perm_pvalue
  }
  hormonal_perm_or_all = rbind(hormonal_perm_or_all, temp)
  cat(i, "\n")
} ; rm(temp, i, z, comorbidity_filter_perm, perm_values_temp)
hormonal_perm_or_all = hormonal_perm_or_all[-1, ] ; rownames(hormonal_perm_or_all) = NULL
hormonal_perm_or_all$comorbidity_filter = paste0("≤", hormonal_perm_or_all$comorbidity_filter)
hormonal_perm_or_all$comorbidity_filter = factor(hormonal_perm_or_all$comorbidity_filter, levels = hormonal_perm_or_all$comorbidity_filter, labels = hormonal_perm_or_all$comorbidity_filter)

cardiovascular_perm_or_all = data.frame(comorbidity_filter = NA,
                                  perm_OR = NA,
                                  perm_P = NA) 
for (i in 1:nrow(cardiovascular_obs_or)) {
  comorbidity_filter_perm = i
  temp = data.frame(comorbidity_filter = rep(i, 1000),
                    perm_OR = NA,
                    perm_P = NA)
  for (z in 1:length(cardiovascular_perm_or)) {
    perm_values_temp = cardiovascular_perm_or[[z]] %>% filter(comorbidity_filter == comorbidity_filter_perm)
    temp[z, "perm_OR"] = perm_values_temp$perm_or
    temp[z, "perm_P"] = perm_values_temp$perm_pvalue
  }
  cardiovascular_perm_or_all = rbind(cardiovascular_perm_or_all, temp)
  cat(i, "\n")
} ; rm(temp, i, z, comorbidity_filter_perm, perm_values_temp)
cardiovascular_perm_or_all = cardiovascular_perm_or_all[-1, ] ; rownames(cardiovascular_perm_or_all) = NULL
cardiovascular_perm_or_all$comorbidity_filter = paste0("≤", cardiovascular_perm_or_all$comorbidity_filter)
cardiovascular_perm_or_all$comorbidity_filter = factor(cardiovascular_perm_or_all$comorbidity_filter, levels = cardiovascular_perm_or_all$comorbidity_filter, labels = cardiovascular_perm_or_all$comorbidity_filter)

ophthalmological_perm_or_all = data.frame(comorbidity_filter = NA,
                                        perm_OR = NA,
                                        perm_P = NA) 
for (i in 1:nrow(ophthalmological_obs_or)) {
  comorbidity_filter_perm = i
  temp = data.frame(comorbidity_filter = rep(i, 1000),
                    perm_OR = NA,
                    perm_P = NA)
  for (z in 1:length(ophthalmological_perm_or)) {
    perm_values_temp = ophthalmological_perm_or[[z]] %>% filter(comorbidity_filter == comorbidity_filter_perm)
    temp[z, "perm_OR"] = perm_values_temp$perm_or
    temp[z, "perm_P"] = perm_values_temp$perm_pvalue
  }
  ophthalmological_perm_or_all = rbind(ophthalmological_perm_or_all, temp)
  cat(i, "\n")
} ; rm(temp, i, z, comorbidity_filter_perm, perm_values_temp)
ophthalmological_perm_or_all = ophthalmological_perm_or_all[-1, ] ; rownames(ophthalmological_perm_or_all) = NULL
ophthalmological_perm_or_all$comorbidity_filter = paste0("≤", ophthalmological_perm_or_all$comorbidity_filter)
ophthalmological_perm_or_all$comorbidity_filter = factor(ophthalmological_perm_or_all$comorbidity_filter, levels = ophthalmological_perm_or_all$comorbidity_filter, labels = ophthalmological_perm_or_all$comorbidity_filter)

immune_perm_or_all = data.frame(comorbidity_filter = NA,
                                          perm_OR = NA,
                                          perm_P = NA) 
for (i in 1:nrow(immune_obs_or)) {
  comorbidity_filter_perm = i
  temp = data.frame(comorbidity_filter = rep(i, 1000),
                    perm_OR = NA,
                    perm_P = NA)
  for (z in 1:length(immune_perm_or)) {
    perm_values_temp = immune_perm_or[[z]] %>% filter(comorbidity_filter == comorbidity_filter_perm)
    temp[z, "perm_OR"] = perm_values_temp$perm_or
    temp[z, "perm_P"] = perm_values_temp$perm_pvalue
  }
  immune_perm_or_all = rbind(immune_perm_or_all, temp)
  cat(i, "\n")
} ; rm(temp, i, z, comorbidity_filter_perm, perm_values_temp)
immune_perm_or_all = immune_perm_or_all[-1, ] ; rownames(immune_perm_or_all) = NULL
immune_perm_or_all$comorbidity_filter = paste0("≤", immune_perm_or_all$comorbidity_filter)
immune_perm_or_all$comorbidity_filter = factor(immune_perm_or_all$comorbidity_filter, levels = immune_perm_or_all$comorbidity_filter, labels = immune_perm_or_all$comorbidity_filter)

neurological_perm_or_all = data.frame(comorbidity_filter = NA,
                                perm_OR = NA,
                                perm_P = NA) 
for (i in 1:nrow(neurological_obs_or)) {
  comorbidity_filter_perm = i
  temp = data.frame(comorbidity_filter = rep(i, 1000),
                    perm_OR = NA,
                    perm_P = NA)
  for (z in 1:length(neurological_perm_or)) {
    perm_values_temp = neurological_perm_or[[z]] %>% filter(comorbidity_filter == comorbidity_filter_perm)
    temp[z, "perm_OR"] = perm_values_temp$perm_or
    temp[z, "perm_P"] = perm_values_temp$perm_pvalue
  }
  neurological_perm_or_all = rbind(neurological_perm_or_all, temp)
  cat(i, "\n")
} ; rm(temp, i, z, comorbidity_filter_perm, perm_values_temp)
neurological_perm_or_all = neurological_perm_or_all[-1, ] ; rownames(neurological_perm_or_all) = NULL
neurological_perm_or_all$comorbidity_filter = paste0("≤", neurological_perm_or_all$comorbidity_filter)
neurological_perm_or_all$comorbidity_filter = factor(neurological_perm_or_all$comorbidity_filter, levels = neurological_perm_or_all$comorbidity_filter, labels = neurological_perm_or_all$comorbidity_filter)

### create plots ###

## Immune ##
## For observed OR: plot OR with 95% CI
immune_obs_or = na.omit(immune_obs_or) ; rownames(immune_obs_or) = NULL
immune_obs_or$comorbidity_filter = paste0("≤", immune_obs_or$comorbidity_filter)
immune_obs_or$comorbidity_filter = factor(immune_obs_or$comorbidity_filter, levels = immune_obs_or$comorbidity_filter, labels = immune_obs_or$comorbidity_filter)
immune_obs_or$status = "observed"
## For permuted ORs: plot 5th, 50th and 95th percentiles
immune_perm_or_all = na.omit(immune_perm_or_all) ; rownames(immune_perm_or_all) = NULL
immune_perm_or_all$comorbidity_filter = factor(immune_perm_or_all$comorbidity_filter, levels = immune_perm_or_all$comorbidity_filter, labels = immune_perm_or_all$comorbidity_filter)
immune_perm_or_all_filt = immune_perm_or_all %>% 
  group_by(comorbidity_filter) %>%
  mutate(percentile_5 = as.numeric(quantile(perm_OR, 0.05)),
         percentile_50 = as.numeric(quantile(perm_OR, 0.50)),
         percentile_95 = as.numeric(quantile(perm_OR, 0.95))) %>%
  filter(perm_OR >= percentile_5 & perm_OR <= percentile_95) %>%
  ungroup() %>%
  dplyr::select(comorbidity_filter, percentile_50, percentile_5,percentile_95) %>%
  distinct() %>%
  mutate(status = "permuted")
# rename columns for rbind later --> will change them after that
colnames(immune_perm_or_all_filt) = c("comorbidity_filter", "obs_or", "ci_95_lower", "ci_95_upper", "status")

immune_obs_perm_combined = rbind(immune_obs_or[, c(1,2,4,5,7)], immune_perm_or_all_filt)
colnames(immune_obs_perm_combined) = c("comorbidity_filter", "OR", "OR_lower_limit", "OR_upper_limit","status")
immune_obs_perm_combined$comorbidity_filter = rep(2:19, 2)
immune_obs_perm_combined$status = factor(immune_obs_perm_combined$status, levels = immune_obs_perm_combined$status, labels = immune_obs_perm_combined$status)

s2_a = ggplot(immune_obs_perm_combined) +
  geom_ribbon(aes(x = comorbidity_filter, ymin = OR_lower_limit, ymax = OR_upper_limit, fill = status),
              alpha = 0.3) +
  scale_fill_manual(values = c("#00BFC4","#F8766D"))+
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 2) +
  geom_line(aes(x = comorbidity_filter, y = OR, color = status), linewidth = 2, show.legend = FALSE) +
  scale_color_manual(values = c("#00BFC4","#F8766D")) +
  scale_x_continuous(breaks = seq(2, 19, 1), labels = paste0("≤", seq(2, 19, 1)), expand = c(0, 0.5)) +
  scale_y_continuous(breaks = seq(0, 16, 2)) +
  labs(title = "Immune diseases") +
  xlab("") +
  ylab("Odds Ratio") +
  labs(fill = "") +
  theme_classic() +
  theme(plot.title = element_text(size = 50, color = "black", family = "Arial", face = "bold"),
        axis.title = element_text(size = 45, color = "black", family = "Arial"),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text.x = element_text(size = 40, color = "black", family = "Arial", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "Arial", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_text(size = 23, color = "black", family = "Arial"),
        legend.text = element_text(size = 40, family = "Arial"),
        legend.position = "top", legend.key.size = unit(2, "cm"))
s2_b = ggplot(immune_obs_or, aes(x = comorbidity_filter, y = nr_drugs_in_lr)) +
  geom_line(group = "comorbidity_filter", linewidth = 1.5) +
  scale_x_discrete(expand = c(0, 0.5)) +
  scale_y_continuous(breaks = seq(0, 800, 200), limits = c(0, 800)) +
  xlab("Comorbidities per Mendelian disease") +
  ylab("Number of drugs") +
  labs(fill = "") +
  theme_classic() +
  theme(axis.title = element_text(size = 45, color = "black", family = "Arial"),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text.x = element_text(size = 40, color = "black", family = "Arial", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "Arial", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"))
ggarrange(s2_a, s2_b, 
          ncol = 1, nrow = 2, heights = c(3,1), align = "v")
ggsave(filename = "S2_immune_diseases.png",
       path = "supplementary_figures/", 
       width = 30, height = 20, device = "png",
       dpi = 300, type = type_compression)

## neurological ##
## For observed OR: plot OR with 95% CI
neurological_obs_or = na.omit(neurological_obs_or) ; rownames(neurological_obs_or) = NULL
neurological_obs_or$comorbidity_filter = paste0("≤", neurological_obs_or$comorbidity_filter)
neurological_obs_or$comorbidity_filter = factor(neurological_obs_or$comorbidity_filter, levels = neurological_obs_or$comorbidity_filter, labels = neurological_obs_or$comorbidity_filter)
neurological_obs_or$status = "observed"
## For permuted ORs: plot 5th, 50th and 95th percentiles
neurological_perm_or_all = na.omit(neurological_perm_or_all) ; rownames(neurological_perm_or_all) = NULL
neurological_perm_or_all$comorbidity_filter = factor(neurological_perm_or_all$comorbidity_filter, levels = neurological_perm_or_all$comorbidity_filter, labels = neurological_perm_or_all$comorbidity_filter)
neurological_perm_or_all_filt = neurological_perm_or_all %>% 
  group_by(comorbidity_filter) %>%
  mutate(percentile_5 = as.numeric(quantile(perm_OR, 0.05)),
         percentile_50 = as.numeric(quantile(perm_OR, 0.50)),
         percentile_95 = as.numeric(quantile(perm_OR, 0.95))) %>%
  filter(perm_OR >= percentile_5 & perm_OR <= percentile_95) %>%
  ungroup() %>%
  dplyr::select(comorbidity_filter, percentile_50, percentile_5,percentile_95) %>%
  distinct() %>%
  mutate(status = "permuted")
# rename columns for rbind later --> will change them after that
colnames(neurological_perm_or_all_filt) = c("comorbidity_filter", "obs_or", "ci_95_lower", "ci_95_upper", "status")

neurological_obs_perm_combined = rbind(neurological_obs_or[, c(1,2,4,5,7)], neurological_perm_or_all_filt)
colnames(neurological_obs_perm_combined) = c("comorbidity_filter", "OR", "OR_lower_limit", "OR_upper_limit","status")
neurological_obs_perm_combined = neurological_obs_perm_combined %>% filter(comorbidity_filter != "≤2") # remove permutations with comorbidity filter ≤2 as all the ORs and Pvalues are 1
neurological_obs_perm_combined$comorbidity_filter = rep(3:15, 2)
neurological_obs_perm_combined$status = factor(neurological_obs_perm_combined$status, levels = neurological_obs_perm_combined$status, labels = neurological_obs_perm_combined$status)

s3_a = ggplot(neurological_obs_perm_combined) +
  geom_ribbon(aes(x = comorbidity_filter, ymin = OR_lower_limit, ymax = OR_upper_limit, fill = status),
              alpha = 0.3) +
  scale_fill_manual(values = c("#00BFC4","#F8766D"))+
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 2) +
  geom_line(aes(x = comorbidity_filter, y = OR, color = status), linewidth = 2, show.legend = FALSE) +
  scale_color_manual(values = c("#00BFC4","#F8766D")) +
  scale_x_continuous(breaks = seq(3, 15, 1), labels = paste0("≤", seq(3, 15, 1)), expand = c(0, 0.5)) +
  scale_y_continuous(breaks = seq(0, 6, 1), limits = c(0, 6)) +
  labs(title = "Neurological diseases") +
  xlab("") +
  ylab("Odds Ratio") +
  labs(fill = "") +
  theme_classic() +
  theme(plot.title = element_text(size = 50, color = "black", family = "Arial", face = "bold"),
        axis.title = element_text(size = 45, color = "black", family = "Arial"),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text.x = element_text(size = 40, color = "black", family = "Arial", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "Arial", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_text(size = 23, color = "black", family = "Arial"),
        legend.text = element_text(size = 40, family = "Arial"),
        legend.position = "top", legend.key.size = unit(2, "cm"))
s3_b = ggplot(neurological_obs_or, aes(x = comorbidity_filter, y = nr_drugs_in_lr)) +
  geom_line(group = "comorbidity_filter", linewidth = 1.5) +
  scale_x_discrete(expand = c(0, 0.5)) +
  scale_y_continuous(breaks = seq(0, 800, 200), limits = c(0, 800)) +
  xlab("Comorbidities per Mendelian disease") +
  ylab("Number of drugs") +
  labs(fill = "") +
  theme_classic() +
  theme(axis.title = element_text(size = 45, color = "black", family = "Arial"),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text.x = element_text(size = 40, color = "black", family = "Arial", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "Arial", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"))
ggarrange(s3_a, s3_b, 
          ncol = 1, nrow = 2, heights = c(3,1), align = "v")
ggsave(filename = "S3_neurological_diseases.png",
       path = "supplementary_figures/", 
       width = 30, height = 20, device = "png",
       dpi = 300, type = type_compression)

## Cancers ##
## For observed OR: plot OR with 95% CI
cancers_obs_or$comorbidity_filter = paste0("≤", cancers_obs_or$comorbidity_filter)
cancers_obs_or$comorbidity_filter = factor(cancers_obs_or$comorbidity_filter, levels = cancers_obs_or$comorbidity_filter, labels = cancers_obs_or$comorbidity_filter)
cancers_obs_or$status = "observed"
## For permuted ORs: plot 5th, 50th and 95th percentiles
cancers_perm_or_all_filt = cancers_perm_or_all %>% 
  group_by(comorbidity_filter) %>%
  mutate(percentile_5 = as.numeric(quantile(perm_OR, 0.05)),
         percentile_50 = as.numeric(quantile(perm_OR, 0.50)),
         percentile_95 = as.numeric(quantile(perm_OR, 0.95))) %>%
  filter(perm_OR >= percentile_5 & perm_OR <= percentile_95) %>%
  ungroup() %>%
  dplyr::select(comorbidity_filter, percentile_50, percentile_5,percentile_95) %>%
  distinct() %>%
  mutate(status = "permuted")
# rename columns for rbind later --> will change them after that
colnames(cancers_perm_or_all_filt) = c("comorbidity_filter", "obs_or", "ci_95_lower", "ci_95_upper", "status")

cancers_obs_perm_combined = rbind(cancers_obs_or[, c(1,2,4,5,7)], cancers_perm_or_all_filt)
colnames(cancers_obs_perm_combined) = c("comorbidity_filter", "OR", "OR_lower_limit", "OR_upper_limit","status")
cancers_obs_perm_combined$comorbidity_filter = rep(1:14, 2)
cancers_obs_perm_combined$status = factor(cancers_obs_perm_combined$status, levels = cancers_obs_perm_combined$status, labels = cancers_obs_perm_combined$status)

s4_a = ggplot(cancers_obs_perm_combined) +
  geom_ribbon(aes(x = comorbidity_filter, ymin = OR_lower_limit, ymax = OR_upper_limit, fill = status),
              alpha = 0.3) +
  scale_fill_manual(values = c("#00BFC4","#F8766D"))+
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 2) +
  geom_line(aes(x = comorbidity_filter, y = OR, color = status), linewidth = 2, show.legend = FALSE) +
  scale_color_manual(values = c("#00BFC4","#F8766D")) +
  scale_x_continuous(breaks = seq(1, 14, 1), labels = paste0("≤", seq(1, 14, 1)), expand = c(0, 0.5)) +
  scale_y_continuous(breaks = seq(0, 3, 0.5)) +
  labs(title = "Neoplasms") +
  xlab("") +
  ylab("Odds Ratio") +
  labs(fill = "") +
  theme_classic() +
  theme(plot.title = element_text(size = 50, color = "black", family = "Arial", face = "bold"),
        axis.title = element_text(size = 45, color = "black", family = "Arial"),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text.x = element_text(size = 40, color = "black", family = "Arial", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "Arial", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_text(size = 23, color = "black", family = "Arial"),
        legend.text = element_text(size = 40, family = "Arial"),
        legend.position = "top", legend.key.size = unit(2, "cm"))
s4_b = ggplot(cancers_obs_or, aes(x = comorbidity_filter, y = nr_drugs_in_lr)) +
  geom_line(group = "comorbidity_filter", linewidth = 1.5) +
  scale_x_discrete(expand = c(0, 0.5)) +
  scale_y_continuous(breaks = seq(0, 800, 200), limits = c(0, 800)) +
  xlab("Comorbidities per Mendelian disease") +
  ylab("Number of drugs") +
  labs(fill = "") +
  theme_classic() +
  theme(axis.title = element_text(size = 45, color = "black", family = "Arial"),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text.x = element_text(size = 40, color = "black", family = "Arial", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "Arial", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"))
ggarrange(s4_a, s4_b, 
          ncol = 1, nrow = 2, heights = c(3,1), align = "v")
ggsave(filename = "S4_neoplasms.png",
       path = "supplementary_figures/", 
       width = 30, height = 20, device = "png",
       dpi = 300, type = type_compression)

## hormonal ##
## For observed OR: plot OR with 95% CI
hormonal_obs_or = na.omit(hormonal_obs_or) ; rownames(hormonal_obs_or) = NULL
hormonal_obs_or$comorbidity_filter = paste0("≤", hormonal_obs_or$comorbidity_filter)
hormonal_obs_or$comorbidity_filter = factor(hormonal_obs_or$comorbidity_filter, levels = hormonal_obs_or$comorbidity_filter, labels = hormonal_obs_or$comorbidity_filter)
hormonal_obs_or$status = "observed"
## For permuted ORs: plot 5th, 50th and 95th percentiles
hormonal_perm_or_all = na.omit(hormonal_perm_or_all) ; rownames(hormonal_perm_or_all) = NULL
hormonal_perm_or_all$comorbidity_filter = factor(hormonal_perm_or_all$comorbidity_filter, levels = hormonal_perm_or_all$comorbidity_filter, labels = hormonal_perm_or_all$comorbidity_filter)
hormonal_perm_or_all_filt = hormonal_perm_or_all %>% 
  group_by(comorbidity_filter) %>%
  mutate(percentile_5 = as.numeric(quantile(perm_OR, 0.05)),
         percentile_50 = as.numeric(quantile(perm_OR, 0.50)),
         percentile_95 = as.numeric(quantile(perm_OR, 0.95))) %>%
  filter(perm_OR >= percentile_5 & perm_OR <= percentile_95) %>%
  ungroup() %>%
  dplyr::select(comorbidity_filter, percentile_50, percentile_5,percentile_95) %>%
  distinct() %>%
  mutate(status = "permuted")
# rename columns for rbind later --> will change them after that
colnames(hormonal_perm_or_all_filt) = c("comorbidity_filter", "obs_or", "ci_95_lower", "ci_95_upper", "status")

hormonal_obs_perm_combined = rbind(hormonal_obs_or[, c(1,2,4,5,7)], hormonal_perm_or_all_filt)
colnames(hormonal_obs_perm_combined) = c("comorbidity_filter", "OR", "OR_lower_limit", "OR_upper_limit","status")
hormonal_obs_perm_combined = hormonal_obs_perm_combined %>% filter(comorbidity_filter != "≤1") # remove it as 95% CI couldn't calculated for the observed OR
hormonal_obs_perm_combined$comorbidity_filter = rep(2:8, 2)
hormonal_obs_perm_combined$status = factor(hormonal_obs_perm_combined$status, levels = hormonal_obs_perm_combined$status, labels = hormonal_obs_perm_combined$status)

s5_a = ggplot(hormonal_obs_perm_combined) +
  geom_ribbon(aes(x = comorbidity_filter, ymin = OR_lower_limit, ymax = OR_upper_limit, fill = status),
              alpha = 0.3) +
  scale_fill_manual(values = c("#00BFC4","#F8766D"))+
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 2) +
  geom_line(aes(x = comorbidity_filter, y = OR, color = status), linewidth = 2, show.legend = FALSE) +
  scale_color_manual(values = c("#00BFC4","#F8766D")) +
  scale_x_continuous(breaks = seq(2, 8, 1), labels = paste0("≤", seq(2, 8, 1)), expand = c(0, 0.5)) +
  scale_y_continuous(breaks = seq(0, 12, 2), limits = c(0, 12)) +
  labs(title = "Hormonal diseases") +
  xlab("") +
  ylab("Odds Ratio") +
  labs(fill = "") +
  theme_classic() +
  theme(plot.title = element_text(size = 50, color = "black", family = "Arial", face = "bold"),
        axis.title = element_text(size = 45, color = "black", family = "Arial"),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text.x = element_text(size = 40, color = "black", family = "Arial", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "Arial", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_text(size = 23, color = "black", family = "Arial"),
        legend.text = element_text(size = 40, family = "Arial"),
        legend.position = "top", legend.key.size = unit(2, "cm"))
s5_b = ggplot(hormonal_obs_or, aes(x = comorbidity_filter, y = nr_drugs_in_lr)) +
  geom_line(group = "comorbidity_filter", linewidth = 1.5) +
  scale_x_discrete(expand = c(0, 0.5)) +
  scale_y_continuous(breaks = seq(0, 800, 200), limits = c(0, 800)) +
  xlab("Comorbidities per Mendelian disease") +
  ylab("Number of drugs") +
  labs(fill = "") +
  theme_classic() +
  theme(axis.title = element_text(size = 45, color = "black", family = "Arial"),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text.x = element_text(size = 40, color = "black", family = "Arial", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "Arial", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"))
ggarrange(s5_a, s5_b, 
          ncol = 1, nrow = 2, heights = c(3,1), align = "v")
ggsave(filename = "S5_hormonal_diseases.png",
       path = "supplementary_figures/", 
       width = 30, height = 20, device = "png",
       dpi = 300, type = type_compression)

## cardiovascular ##
## For observed OR: plot OR with 95% CI
cardiovascular_obs_or = na.omit(cardiovascular_obs_or) ; rownames(cardiovascular_obs_or) = NULL
cardiovascular_obs_or$comorbidity_filter = paste0("≤", cardiovascular_obs_or$comorbidity_filter)
cardiovascular_obs_or$comorbidity_filter = factor(cardiovascular_obs_or$comorbidity_filter, levels = cardiovascular_obs_or$comorbidity_filter, labels = cardiovascular_obs_or$comorbidity_filter)
cardiovascular_obs_or$status = "observed"
## For permuted ORs: plot 5th, 50th and 95th percentiles
cardiovascular_perm_or_all = na.omit(cardiovascular_perm_or_all) ; rownames(cardiovascular_perm_or_all) = NULL
cardiovascular_perm_or_all$comorbidity_filter = factor(cardiovascular_perm_or_all$comorbidity_filter, levels = cardiovascular_perm_or_all$comorbidity_filter, labels = cardiovascular_perm_or_all$comorbidity_filter)
cardiovascular_perm_or_all_filt = cardiovascular_perm_or_all %>% 
  group_by(comorbidity_filter) %>%
  mutate(percentile_5 = as.numeric(quantile(perm_OR, 0.05)),
         percentile_50 = as.numeric(quantile(perm_OR, 0.50)),
         percentile_95 = as.numeric(quantile(perm_OR, 0.95))) %>%
  filter(perm_OR >= percentile_5 & perm_OR <= percentile_95) %>%
  ungroup() %>%
  dplyr::select(comorbidity_filter, percentile_50, percentile_5,percentile_95) %>%
  distinct() %>%
  mutate(status = "permuted")
# rename columns for rbind later --> will change them after that
colnames(cardiovascular_perm_or_all_filt) = c("comorbidity_filter", "obs_or", "ci_95_lower", "ci_95_upper", "status")

cardiovascular_obs_perm_combined = rbind(cardiovascular_obs_or[, c(1,2,4,5,7)], cardiovascular_perm_or_all_filt)
colnames(cardiovascular_obs_perm_combined) = c("comorbidity_filter", "OR", "OR_lower_limit", "OR_upper_limit","status")
cardiovascular_obs_perm_combined$comorbidity_filter = rep(1:4, 2)
cardiovascular_obs_perm_combined$status = factor(cardiovascular_obs_perm_combined$status, levels = cardiovascular_obs_perm_combined$status, labels = cardiovascular_obs_perm_combined$status)

s6_a = ggplot(cardiovascular_obs_perm_combined) +
  geom_ribbon(aes(x = comorbidity_filter, ymin = OR_lower_limit, ymax = OR_upper_limit, fill = status),
              alpha = 0.3) +
  scale_fill_manual(values = c("#00BFC4","#F8766D"))+
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 2) +
  geom_line(aes(x = comorbidity_filter, y = OR, color = status), linewidth = 2, show.legend = FALSE) +
  scale_color_manual(values = c("#00BFC4","#F8766D")) +
  scale_x_continuous(breaks = seq(1, 4, 1), labels = paste0("≤", seq(1, 4, 1)), expand = c(0, 0.5)) +
  scale_y_continuous(breaks = seq(0, 4, 1), limits = c(0, 4)) +
  labs(title = "Cardiovascular diseases") +
  xlab("") +
  ylab("Odds Ratio") +
  labs(fill = "") +
  theme_classic() +
  theme(plot.title = element_text(size = 50, color = "black", family = "Arial", face = "bold"),
        axis.title = element_text(size = 45, color = "black", family = "Arial"),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text.x = element_text(size = 40, color = "black", family = "Arial", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "Arial", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_text(size = 23, color = "black", family = "Arial"),
        legend.text = element_text(size = 40, family = "Arial"),
        legend.position = "top", legend.key.size = unit(2, "cm"))
s6_b = ggplot(ophthalmological_obs_or, aes(x = comorbidity_filter, y = nr_drugs_in_lr)) +
  geom_line(group = "comorbidity_filter", linewidth = 1.5) +
  scale_x_discrete(expand = c(0, 0.5)) +
  scale_y_continuous(breaks = seq(0, 800, 200), limits = c(0, 800)) +
  xlab("Comorbidities per Mendelian disease") +
  ylab("Number of drugs") +
  labs(fill = "") +
  theme_classic() +
  theme(axis.title = element_text(size = 45, color = "black", family = "Arial"),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text.x = element_text(size = 40, color = "black", family = "Arial", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "Arial", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"))
ggarrange(s6_a, s6_b, 
          ncol = 1, nrow = 2, heights = c(3,1), align = "v")
ggsave(filename = "S6_cardiovascular_diseases.png",
       path = "supplementary_figures/", 
       width = 30, height = 15, device = "png",
       dpi = 300, type = type_compression)

## ophthalmological ##
## For observed OR: plot OR with 95% CI
ophthalmological_obs_or = na.omit(ophthalmological_obs_or) ; rownames(ophthalmological_obs_or) = NULL
ophthalmological_obs_or$comorbidity_filter = paste0("≤", ophthalmological_obs_or$comorbidity_filter)
ophthalmological_obs_or$comorbidity_filter = factor(ophthalmological_obs_or$comorbidity_filter, levels = ophthalmological_obs_or$comorbidity_filter, labels = ophthalmological_obs_or$comorbidity_filter)
ophthalmological_obs_or$status = "observed"
## For permuted ORs: plot 5th, 50th and 95th percentiles
ophthalmological_perm_or_all = na.omit(ophthalmological_perm_or_all) ; rownames(ophthalmological_perm_or_all) = NULL
ophthalmological_perm_or_all$comorbidity_filter = factor(ophthalmological_perm_or_all$comorbidity_filter, levels = ophthalmological_perm_or_all$comorbidity_filter, labels = ophthalmological_perm_or_all$comorbidity_filter)
ophthalmological_perm_or_all_filt = ophthalmological_perm_or_all %>% 
  group_by(comorbidity_filter) %>%
  mutate(percentile_5 = as.numeric(quantile(perm_OR, 0.05)),
         percentile_50 = as.numeric(quantile(perm_OR, 0.50)),
         percentile_95 = as.numeric(quantile(perm_OR, 0.95))) %>%
  filter(perm_OR >= percentile_5 & perm_OR <= percentile_95) %>%
  ungroup() %>%
  dplyr::select(comorbidity_filter, percentile_50, percentile_5,percentile_95) %>%
  distinct() %>%
  mutate(status = "permuted")
# rename columns for rbind later --> will change them after that
colnames(ophthalmological_perm_or_all_filt) = c("comorbidity_filter", "obs_or", "ci_95_lower", "ci_95_upper", "status")

ophthalmological_obs_perm_combined = rbind(ophthalmological_obs_or[, c(1,2,4,5,7)], ophthalmological_perm_or_all_filt)
colnames(ophthalmological_obs_perm_combined) = c("comorbidity_filter", "OR", "OR_lower_limit", "OR_upper_limit","status")
ophthalmological_obs_perm_combined$comorbidity_filter = rep(1:5, 2)
ophthalmological_obs_perm_combined$status = factor(ophthalmological_obs_perm_combined$status, levels = ophthalmological_obs_perm_combined$status, labels = ophthalmological_obs_perm_combined$status)

s7_a = ggplot(ophthalmological_obs_perm_combined) +
  geom_ribbon(aes(x = comorbidity_filter, ymin = OR_lower_limit, ymax = OR_upper_limit, fill = status),
              alpha = 0.3) +
  scale_fill_manual(values = c("#00BFC4","#F8766D"))+
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 2) +
  geom_line(aes(x = comorbidity_filter, y = OR, color = status), linewidth = 2, show.legend = FALSE) +
  scale_color_manual(values = c("#00BFC4","#F8766D")) +
  scale_x_continuous(breaks = seq(1, 5, 1), labels = paste0("≤", seq(1, 5, 1)), expand = c(0, 0.5)) +
  scale_y_continuous(breaks = seq(0, 125, 25), limits = c(0, 125)) +
  labs(title = "Ophthalmological diseases") +
  xlab("") +
  ylab("Odds Ratio") +
  labs(fill = "") +
  theme_classic() +
  theme(plot.title = element_text(size = 50, color = "black", family = "Arial", face = "bold"),
        axis.title = element_text(size = 45, color = "black", family = "Arial"),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text.x = element_text(size = 40, color = "black", family = "Arial", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "Arial", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_text(size = 23, color = "black", family = "Arial"),
        legend.text = element_text(size = 40, family = "Arial"),
        legend.position = "top", legend.key.size = unit(2, "cm"))
s7_b = ggplot(ophthalmological_obs_or, aes(x = comorbidity_filter, y = nr_drugs_in_lr)) +
  geom_line(group = "comorbidity_filter", linewidth = 1.5) +
  scale_x_discrete(expand = c(0, 0.5)) +
  scale_y_continuous(breaks = seq(0, 800, 200), limits = c(0, 800)) +
  xlab("Comorbidities per Mendelian disease") +
  ylab("Number of drugs") +
  labs(fill = "") +
  theme_classic() +
  theme(axis.title = element_text(size = 45, color = "black", family = "Arial"),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text.x = element_text(size = 40, color = "black", family = "Arial", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "Arial", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"))
ggarrange(s7_a, s7_b, 
          ncol = 1, nrow = 2, heights = c(3,1), align = "v")
ggsave(filename = "S7_ophthalmological_diseases.png",
       path = "supplementary_figures/", 
       width = 30, height = 15, device = "png",
       dpi = 300, type = type_compression)

### -- per clinical trial phase -- ###

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

phase123_obs_or = data.frame(comorbidity_filter = 10:61,
                             obs_or = 0,
                             obs_pvalue = 0,
                             ci95_lower = 0,
                             ci95_upper = 0)
phase1_obs_or = data.frame(comorbidity_filter = 10:61,
                           obs_or = 0,
                           obs_pvalue = 0,
                           ci95_lower = 0,
                           ci95_upper = 0)
phase2_obs_or = data.frame(comorbidity_filter = 10:61,
                           obs_or = 0,
                           obs_pvalue = 0,
                           ci95_lower = 0,
                           ci95_upper = 0)
phase3_obs_or = data.frame(comorbidity_filter = 10:61,
                           obs_or = 0,
                           obs_pvalue = 0,
                           ci95_lower = 0,
                           ci95_upper = 0)
approved_obs_or = data.frame(comorbidity_filter = 10:61,
                             obs_or = 0,
                             obs_pvalue = 0,
                             ci95_lower = 0,
                             ci95_upper = 0)

for (i in 1:nrow(phase123_obs_or)) {
  
  comorbidity_filter = phase123_obs_or[i, "comorbidity_filter"]
  
  mds_to_keep = nr_comorbidities_per_md %>% filter(nr_comorbidities <= comorbidity_filter)
  md_cd_comorbidities_filt = md_cd_comorbidities %>% filter(mendelian_disease %in% mds_to_keep$mendelian_disease)
  md_genes_filt = md_genes %>% filter(mendelian_disease %in% mds_to_keep$mendelian_disease)
  drugs_targeting_md_genes = left_join(md_genes_filt, db_drug_targets, by = c("causal_gene" = "drug_target")) %>%
    na.omit() %>%
    dplyr::select(db_id) %>%
    distinct()
  
  unique_cd = unique(md_cd_comorbidities_filt$complex_disease)
  
  # create logistic regression input
  log_input_ct_all_phases = data.frame(db_id = rep(drugs_targeting_md_genes$db_id, each = length(unique_cd)))
  log_input_ct_all_phases$complex_disease = unique_cd
  # add information about disease category
  log_input_ct_all_phases = left_join(log_input_ct_all_phases, complex_disease_categories, by = "complex_disease")
  # add information about number of targets for each drug
  log_input_ct_all_phases = left_join(log_input_ct_all_phases, drugs_nr_targets, by = "db_id")  
  # add information about comorbidity
  x = left_join(md_cd_comorbidities_filt, md_genes_filt, by = "mendelian_disease")
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
  phase123_obs_or[i, "obs_or"] = exp(summary(log_reg_ct_phase123)$coefficients["comorbidity", "Estimate"])
  phase123_obs_or[i, "obs_pvalue"] = summary(log_reg_ct_phase123)$coefficients["comorbidity", 4]
  phase123_obs_or[i, "ci95_lower"] = exp(confint(log_reg_ct_phase123))["comorbidity", "2.5 %"]
  phase123_obs_or[i, "ci95_upper"] = exp(confint(log_reg_ct_phase123))["comorbidity", "97.5 %"] 
  
  # phase 1
  log_reg_ct_phase1 = glm(indicated_investigated_phase1 ~ total_targets + disease_category + comorbidity,
                          data = log_input_ct_all_phases, 
                          family = binomial())
  phase1_obs_or[i, "obs_or"] = exp(summary(log_reg_ct_phase1)$coefficients["comorbidity", "Estimate"])
  phase1_obs_or[i, "obs_pvalue"] = summary(log_reg_ct_phase1)$coefficients["comorbidity", 4]
  phase1_obs_or[i, "ci95_lower"] = exp(confint(log_reg_ct_phase1))["comorbidity", "2.5 %"]
  phase1_obs_or[i, "ci95_upper"] = exp(confint(log_reg_ct_phase1))["comorbidity", "97.5 %"] 
  # phase 2
  log_reg_ct_phase2 = glm(indicated_investigated_phase2 ~ total_targets + disease_category + comorbidity,
                          data = log_input_ct_all_phases, 
                          family = binomial())
  phase2_obs_or[i, "obs_or"] = exp(summary(log_reg_ct_phase2)$coefficients["comorbidity", "Estimate"])
  phase2_obs_or[i, "obs_pvalue"] = summary(log_reg_ct_phase2)$coefficients["comorbidity", 4]
  phase2_obs_or[i, "ci95_lower"] = exp(confint(log_reg_ct_phase2))["comorbidity", "2.5 %"]
  phase2_obs_or[i, "ci95_upper"] = exp(confint(log_reg_ct_phase2))["comorbidity", "97.5 %"] 
  # phase 3
  log_reg_ct_phase3 = glm(indicated_investigated_phase3 ~ total_targets + disease_category + comorbidity,
                          data = log_input_ct_all_phases, 
                          family = binomial())
  phase3_obs_or[i, "obs_or"] = exp(summary(log_reg_ct_phase3)$coefficients["comorbidity", "Estimate"])
  phase3_obs_or[i, "obs_pvalue"] = summary(log_reg_ct_phase3)$coefficients["comorbidity", 4]
  phase3_obs_or[i, "ci95_lower"] = exp(confint(log_reg_ct_phase3))["comorbidity", "2.5 %"]
  phase3_obs_or[i, "ci95_upper"] = exp(confint(log_reg_ct_phase3))["comorbidity", "97.5 %"] 
  # approved
  log_reg_ct_approved = glm(indicated_investigated_approved ~ total_targets + disease_category + comorbidity,
                            data = log_input_ct_all_phases, 
                            family = binomial())
  approved_obs_or[i, "obs_or"] = exp(summary(log_reg_ct_approved)$coefficients["comorbidity", "Estimate"])
  approved_obs_or[i, "obs_pvalue"] = summary(log_reg_ct_approved)$coefficients["comorbidity", 4]
  approved_obs_or[i, "ci95_lower"] = exp(confint(log_reg_ct_approved))["comorbidity", "2.5 %"]
  approved_obs_or[i, "ci95_upper"] = exp(confint(log_reg_ct_approved))["comorbidity", "97.5 %"] 
  
  cat(i, "\n")
}
rm(mds_to_keep, comorbidity_filter, md_cd_comorbidities_filt, md_genes_filt, drugs_targeting_md_genes, unique_cd, log_input_ct_all_phases, x, 
   log_reg_ct_phase123, log_reg_ct_phase1, log_reg_ct_phase2, log_reg_ct_phase3, log_reg_ct_approved)


## permutations ##
ct_phase123_perm_or = vector("list", 1000) ; names(ct_phase123_perm_or) = paste0("permutation_", 1:1000)
ct_phase1_perm_or = vector("list", 1000) ; names(ct_phase1_perm_or) = paste0("permutation_", 1:1000)
ct_phase2_perm_or = vector("list", 1000) ; names(ct_phase2_perm_or) = paste0("permutation_", 1:1000)
ct_phase3_perm_or = vector("list", 1000) ; names(ct_phase3_perm_or) = paste0("permutation_", 1:1000)
approved_perm_or = vector("list", 1000) ; names(approved_perm_or) = paste0("permutation_", 1:1000)

for (permutation in 1:1000) {
  
  ct_phase123_perm_or_temp = data.frame(comorbidity_filter = 10:61,
                                        perm_or = 0)
  ct_phase1_perm_or_temp = data.frame(comorbidity_filter = 10:61,
                                      perm_or = 0)
  ct_phase2_perm_or_temp = data.frame(comorbidity_filter = 10:61,
                                      perm_or = 0)
  ct_phase3_perm_or_temp = data.frame(comorbidity_filter = 10:61,
                                      perm_or = 0)
  approved_perm_or_temp = data.frame(comorbidity_filter = 10:61,
                                     perm_or = 0)
  
  for (i in 1:nrow(ct_phase123_perm_or_temp)) {
    
    comorbidity_filter = ct_phase123_perm_or_temp[i, "comorbidity_filter"]
    mds_to_keep = nr_comorbidities_per_md %>% filter(nr_comorbidities <= comorbidity_filter)
    
    md_cd_comorbidities_filt = md_cd_comorbidities %>% filter(mendelian_disease %in% mds_to_keep$mendelian_disease)
    md_cd_comorbidities_filt_matrix = as.data.frame.matrix(table(md_cd_comorbidities_filt[, 2:1]))
    md_cd_comorbidities_filt_matrix = md_cd_comorbidities_filt_matrix %>% mutate(complex_disease = rownames(md_cd_comorbidities_filt_matrix), .before = colnames(md_cd_comorbidities_filt_matrix)[1])
    rownames(md_cd_comorbidities_matrix) = NULL
    md_cd_comorbidities_filt_shuffled = shuffle(md_cd_comorbidities_filt_matrix)
    
    md_genes_filt = md_genes %>% filter(mendelian_disease %in% mds_to_keep$mendelian_disease)
    drugs_targeting_md_genes = left_join(md_genes_filt, db_drug_targets, by = c("causal_gene" = "drug_target")) %>%
      na.omit() %>%
      dplyr::select(db_id) %>%
      distinct()
    
    unique_cd = unique(md_cd_comorbidities_filt_shuffled$complex_disease)
    
    # create logistic regression input - each row is a drug-complex disease pair 
    log_input_perm_ct_all_phases = data.frame(db_id = rep(drugs_targeting_md_genes$db_id, each = length(unique_cd)))
    log_input_perm_ct_all_phases$complex_disease = unique_cd
    # add information about disease category
    log_input_perm_ct_all_phases = left_join(log_input_perm_ct_all_phases, complex_disease_categories, by = "complex_disease")
    # add information about number of targets for each drug
    log_input_perm_ct_all_phases = left_join(log_input_perm_ct_all_phases, drugs_nr_targets, by = "db_id")  
    # add information about comorbidity
    x = left_join(md_cd_comorbidities_filt_shuffled, md_genes_filt, by = "mendelian_disease")
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
    ct_phase123_perm_or_temp[i, "perm_or"] = exp(summary(log_reg_perm_ct_phase123)$coefficients["comorbidity","Estimate"])
    
    # phase 1
    log_reg_perm_ct_phase1 = glm(indicated_investigated_phase1 ~ total_targets + disease_category + comorbidity,
                                 data = log_input_perm_ct_all_phases, 
                                 family = binomial())
    ct_phase1_perm_or_temp[i, "perm_or"] = exp(summary(log_reg_perm_ct_phase1)$coefficients["comorbidity","Estimate"])
    
    # phase 2
    log_reg_perm_ct_phase2 = glm(indicated_investigated_phase2 ~ total_targets + disease_category + comorbidity,
                                 data = log_input_perm_ct_all_phases, 
                                 family = binomial())
    ct_phase2_perm_or_temp[i, "perm_or"] = exp(summary(log_reg_perm_ct_phase2)$coefficients["comorbidity","Estimate"])
    
    # phase 3
    log_reg_perm_ct_phase3 = glm(indicated_investigated_phase3 ~ total_targets + disease_category + comorbidity,
                                 data = log_input_perm_ct_all_phases, 
                                 family = binomial())
    ct_phase3_perm_or_temp[i, "perm_or"] = exp(summary(log_reg_perm_ct_phase3)$coefficients["comorbidity","Estimate"])
    
    # approved
    log_reg_perm_ct_approved = glm(indicated_investigated_approved ~ total_targets + disease_category + comorbidity,
                                   data = log_input_perm_ct_all_phases, 
                                   family = binomial())
    approved_perm_or_temp[i, "perm_or"] = exp(summary(log_reg_perm_ct_approved)$coefficients["comorbidity","Estimate"])
    
    # track progress
    cat(permutation, "-", i, "-", nrow(ct_phase123_perm_or_temp), "\n")
  }
  
  ct_phase123_perm_or[[paste0("permutation_", permutation)]] = ct_phase123_perm_or_temp
  ct_phase1_perm_or[[paste0("permutation_", permutation)]] = ct_phase1_perm_or_temp
  ct_phase2_perm_or[[paste0("permutation_", permutation)]] = ct_phase2_perm_or_temp
  ct_phase3_perm_or[[paste0("permutation_", permutation)]] = ct_phase3_perm_or_temp
  approved_perm_or[[paste0("permutation_", permutation)]] = approved_perm_or_temp
}

phase123_perm_or_all = data.frame(comorbidity_filter = NA,
                                  perm_OR = NA)
for (i in min(phase123_obs_or$comorbidity_filter):max(phase123_obs_or$comorbidity_filter)) {
  comorbidity_filter_perm = i
  temp = data.frame(comorbidity_filter = rep(i, 1000),
                    perm_OR = NA)
  for (z in 1:length(ct_phase123_perm_or)) {
    perm_values_temp = ct_phase123_perm_or[[z]] %>% filter(comorbidity_filter == comorbidity_filter_perm)
    temp[z, "perm_OR"] = perm_values_temp$perm_or
  }
  phase123_perm_or_all = rbind(phase123_perm_or_all, temp)
  cat(i, "\n")
} ; rm(temp, i, z, comorbidity_filter_perm, perm_values_temp)
phase123_perm_or_all = phase123_perm_or_all[-1, ] ; rownames(phase123_perm_or_all) = NULL
phase123_perm_or_all$comorbidity_filter = paste0("≤", phase123_perm_or_all$comorbidity_filter)
phase123_perm_or_all$comorbidity_filter = factor(phase123_perm_or_all$comorbidity_filter, levels = phase123_perm_or_all$comorbidity_filter, labels = phase123_perm_or_all$comorbidity_filter)

phase1_perm_or_all = data.frame(comorbidity_filter = NA,
                                  perm_OR = NA)
for (i in min(phase1_obs_or$comorbidity_filter):max(phase1_obs_or$comorbidity_filter)) {
  comorbidity_filter_perm = i
  temp = data.frame(comorbidity_filter = rep(i, 1000),
                    perm_OR = NA)
  for (z in 1:length(ct_phase1_perm_or)) {
    perm_values_temp = ct_phase1_perm_or[[z]] %>% filter(comorbidity_filter == comorbidity_filter_perm)
    temp[z, "perm_OR"] = perm_values_temp$perm_or
  }
  phase1_perm_or_all = rbind(phase1_perm_or_all, temp)
  cat(i, "\n")
} ; rm(temp, i, z, comorbidity_filter_perm, perm_values_temp)
phase1_perm_or_all = phase1_perm_or_all[-1, ] ; rownames(phase1_perm_or_all) = NULL
phase1_perm_or_all$comorbidity_filter = paste0("≤", phase1_perm_or_all$comorbidity_filter)
phase1_perm_or_all$comorbidity_filter = factor(phase1_perm_or_all$comorbidity_filter, levels = phase1_perm_or_all$comorbidity_filter, labels = phase1_perm_or_all$comorbidity_filter)

phase2_perm_or_all = data.frame(comorbidity_filter = NA,
                                perm_OR = NA)
for (i in min(phase2_obs_or$comorbidity_filter):max(phase2_obs_or$comorbidity_filter)) {
  comorbidity_filter_perm = i
  temp = data.frame(comorbidity_filter = rep(i, 1000),
                    perm_OR = NA)
  for (z in 1:length(ct_phase2_perm_or)) {
    perm_values_temp = ct_phase2_perm_or[[z]] %>% filter(comorbidity_filter == comorbidity_filter_perm)
    temp[z, "perm_OR"] = perm_values_temp$perm_or
  }
  phase2_perm_or_all = rbind(phase2_perm_or_all, temp)
  cat(i, "\n")
} ; rm(temp, i, z, comorbidity_filter_perm, perm_values_temp)
phase2_perm_or_all = phase2_perm_or_all[-1, ] ; rownames(phase2_perm_or_all) = NULL
phase2_perm_or_all$comorbidity_filter = paste0("≤", phase2_perm_or_all$comorbidity_filter)
phase2_perm_or_all$comorbidity_filter = factor(phase2_perm_or_all$comorbidity_filter, levels = phase2_perm_or_all$comorbidity_filter, labels = phase2_perm_or_all$comorbidity_filter)

phase3_perm_or_all = data.frame(comorbidity_filter = NA,
                                perm_OR = NA)
for (i in min(phase3_obs_or$comorbidity_filter):max(phase3_obs_or$comorbidity_filter)) {
  comorbidity_filter_perm = i
  temp = data.frame(comorbidity_filter = rep(i, 1000),
                    perm_OR = NA)
  for (z in 1:length(ct_phase3_perm_or)) {
    perm_values_temp = ct_phase3_perm_or[[z]] %>% filter(comorbidity_filter == comorbidity_filter_perm)
    temp[z, "perm_OR"] = perm_values_temp$perm_or
  }
  phase3_perm_or_all = rbind(phase3_perm_or_all, temp)
  cat(i, "\n")
} ; rm(temp, i, z, comorbidity_filter_perm, perm_values_temp)
phase3_perm_or_all = phase3_perm_or_all[-1, ] ; rownames(phase3_perm_or_all) = NULL
phase3_perm_or_all$comorbidity_filter = paste0("≤", phase3_perm_or_all$comorbidity_filter)
phase3_perm_or_all$comorbidity_filter = factor(phase3_perm_or_all$comorbidity_filter, levels = phase3_perm_or_all$comorbidity_filter, labels = phase3_perm_or_all$comorbidity_filter)

approved_perm_or_all = data.frame(comorbidity_filter = NA,
                                perm_OR = NA)
for (i in min(approved_obs_or$comorbidity_filter):max(approved_obs_or$comorbidity_filter)) {
  comorbidity_filter_perm = i
  temp = data.frame(comorbidity_filter = rep(i, 1000),
                    perm_OR = NA)
  for (z in 1:length(approved_perm_or)) {
    perm_values_temp = approved_perm_or[[z]] %>% filter(comorbidity_filter == comorbidity_filter_perm)
    temp[z, "perm_OR"] = perm_values_temp$perm_or
  }
  approved_perm_or_all = rbind(approved_perm_or_all, temp)
  cat(i, "\n")
} ; rm(temp, i, z, comorbidity_filter_perm, perm_values_temp)
approved_perm_or_all = approved_perm_or_all[-1, ] ; rownames(approved_perm_or_all) = NULL
approved_perm_or_all$comorbidity_filter = paste0("≤", approved_perm_or_all$comorbidity_filter)
approved_perm_or_all$comorbidity_filter = factor(approved_perm_or_all$comorbidity_filter, levels = approved_perm_or_all$comorbidity_filter, labels = approved_perm_or_all$comorbidity_filter)

## Visualizations ##
# Phase123
## For observed OR: plot OR with 95% CI
phase123_obs_or = na.omit(phase123_obs_or) ; rownames(phase123_obs_or) = NULL
phase123_obs_or$comorbidity_filter = paste0("≤", phase123_obs_or$comorbidity_filter)
phase123_obs_or$comorbidity_filter = factor(phase123_obs_or$comorbidity_filter, levels = phase123_obs_or$comorbidity_filter, labels = phase123_obs_or$comorbidity_filter)
phase123_obs_or$status = "observed"
## For permuted ORs: plot 5th, 50th and 95th percentiles
phase123_perm_or_all = na.omit(phase123_perm_or_all) ; rownames(phase123_perm_or_all) = NULL
phase123_perm_or_all$comorbidity_filter = factor(phase123_perm_or_all$comorbidity_filter, levels = phase123_perm_or_all$comorbidity_filter, labels = phase123_perm_or_all$comorbidity_filter)
phase123_perm_or_all_filt = phase123_perm_or_all %>% 
  group_by(comorbidity_filter) %>%
  mutate(percentile_5 = as.numeric(quantile(perm_OR, 0.05)),
         percentile_50 = as.numeric(quantile(perm_OR, 0.50)),
         percentile_95 = as.numeric(quantile(perm_OR, 0.95))) %>%
  filter(perm_OR >= percentile_5 & perm_OR <= percentile_95) %>%
  ungroup() %>%
  dplyr::select(comorbidity_filter, percentile_50, percentile_5,percentile_95) %>%
  distinct() %>%
  mutate(status = "permuted")
# rename columns for rbind later --> will change them after that
colnames(phase123_perm_or_all_filt) = c("comorbidity_filter", "obs_or", "ci95_lower", "ci95_upper", "status")

phase123_obs_perm_combined = rbind(phase123_obs_or[, c(1,2,4,5,6)], phase123_perm_or_all_filt)
colnames(phase123_obs_perm_combined) = c("comorbidity_filter", "OR", "OR_lower_limit", "OR_upper_limit", "status")
phase123_obs_perm_combined$comorbidity_filter = rep(10:61, 2)
phase123_obs_perm_combined$status = factor(phase123_obs_perm_combined$status, levels = phase123_obs_perm_combined$status, labels = phase123_obs_perm_combined$status)

s8_a = ggplot(phase123_obs_perm_combined) +
  geom_ribbon(aes(x = comorbidity_filter, ymin = OR_lower_limit, ymax = OR_upper_limit, fill = status),
              alpha = 0.3) +
  scale_fill_manual(values = c("#00BFC4","#F8766D"))+
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 2) +
  geom_line(aes(x = comorbidity_filter, y = OR, color = status), linewidth = 2, show.legend = FALSE) +
  scale_color_manual(values = c("#00BFC4","#F8766D")) +
  xlab("") +
  scale_x_continuous(breaks = seq(10, 61, 3), labels = paste0("≤", seq(10, 61, 3)), expand = c(0, 0.5)) +
  ylab("Odds Ratio") +
  labs(title = "Phase I/II/III", fill = "") +
  theme_classic() +
  theme(plot.title = element_text(size = 70, color = "black", family = "Arial", face = "bold"),
        axis.title = element_text(size = 45, color = "black", family = "Arial"),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text.x = element_text(size = 40, color = "black", family = "Arial", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "Arial", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_text(size = 23, color = "black", family = "Arial"),
        legend.text = element_text(size = 40, family = "Arial"),
        legend.position = "top", legend.key.size = unit(2, "cm"))
s8_b = ggplot(log_reg_observed_nrcomorbidities[2:53, ], aes(x = comorbidity_filter, y = logr_drugs_included)) +
  geom_line(group = "comorbidity_filter", linewidth = 1.5) +
  scale_x_discrete(breaks = paste0("≤", seq(10, 61, 3)), , expand = c(0, 0.5)) +
  scale_y_continuous(breaks = seq(0, 800, 200)) +
  xlab("Comorbidities per Mendelian disease") +
  ylab("Number of drugs") +
  labs(fill = "") +
  theme_classic() +
  theme(axis.title = element_text(size = 45, color = "black", family = "Arial"),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text.x = element_text(size = 40, color = "black", family = "Arial", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "Arial", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"))
ggarrange(s8_a, s8_b, 
          ncol = 1, nrow = 2, heights = c(3,1), align = "v")
ggsave(filename = "S8_phase_I_II_III.png",
       path = "supplementary_figures/", 
       width = 45, height = 30, device = "png",
       dpi = 300, type = type_compression)

# Phase1
## For observed OR: plot OR with 95% CI
phase1_obs_or = na.omit(phase1_obs_or) ; rownames(phase1_obs_or) = NULL
phase1_obs_or$comorbidity_filter = paste0("≤", phase1_obs_or$comorbidity_filter)
phase1_obs_or$comorbidity_filter = factor(phase1_obs_or$comorbidity_filter, levels = phase1_obs_or$comorbidity_filter, labels = phase1_obs_or$comorbidity_filter)
phase1_obs_or$status = "observed"
## For permuted ORs: plot 5th, 50th and 95th percentiles
phase1_perm_or_all = na.omit(phase1_perm_or_all) ; rownames(phase1_perm_or_all) = NULL
phase1_perm_or_all$comorbidity_filter = factor(phase1_perm_or_all$comorbidity_filter, levels = phase1_perm_or_all$comorbidity_filter, labels = phase1_perm_or_all$comorbidity_filter)
phase1_perm_or_all_filt = phase1_perm_or_all %>% 
  group_by(comorbidity_filter) %>%
  mutate(percentile_5 = as.numeric(quantile(perm_OR, 0.05)),
         percentile_50 = as.numeric(quantile(perm_OR, 0.50)),
         percentile_95 = as.numeric(quantile(perm_OR, 0.95))) %>%
  filter(perm_OR >= percentile_5 & perm_OR <= percentile_95) %>%
  ungroup() %>%
  dplyr::select(comorbidity_filter, percentile_50, percentile_5,percentile_95) %>%
  distinct() %>%
  mutate(status = "permuted")
# rename columns for rbind later --> will change them after that
colnames(phase1_perm_or_all_filt) = c("comorbidity_filter", "obs_or", "ci95_lower", "ci95_upper", "status")

phase1_obs_perm_combined = rbind(phase1_obs_or[, c(1,2,4,5,6)], phase1_perm_or_all_filt)
colnames(phase1_obs_perm_combined) = c("comorbidity_filter", "OR", "OR_lower_limit", "OR_upper_limit", "status")
phase1_obs_perm_combined$comorbidity_filter = rep(14:61, 2)
phase1_obs_perm_combined$status = factor(phase1_obs_perm_combined$status, levels = phase1_obs_perm_combined$status, labels = phase1_obs_perm_combined$status)

s9_a = ggplot(phase1_obs_perm_combined) +
  geom_ribbon(aes(x = comorbidity_filter, ymin = OR_lower_limit, ymax = OR_upper_limit, fill = status),
              alpha = 0.3) +
  scale_fill_manual(values = c("#00BFC4","#F8766D"))+
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 2) +
  geom_line(aes(x = comorbidity_filter, y = OR, color = status), linewidth = 2, show.legend = FALSE) +
  scale_color_manual(values = c("#00BFC4","#F8766D")) +
  xlab("") +
  scale_x_continuous(breaks = c(seq(14, 61, 3), 61), labels = paste0("≤", c(seq(14, 61, 3), 61)), expand = c(0, 0.5)) +
  ylab("Odds Ratio") +
  labs(title = "Phase I", fill = "") +
  theme_classic() +
  theme(plot.title = element_text(size = 70, color = "black", family = "Arial", face = "bold"),
        axis.title = element_text(size = 45, color = "black", family = "Arial"),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text.x = element_text(size = 40, color = "black", family = "Arial", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "Arial", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_text(size = 23, color = "black", family = "Arial"),
        legend.text = element_text(size = 40, family = "Arial"),
        legend.position = "top", legend.key.size = unit(2, "cm"))
s9_b = ggplot(log_reg_observed_nrcomorbidities[6:53, ], aes(x = comorbidity_filter, y = logr_drugs_included)) +
  geom_line(group = "comorbidity_filter", linewidth = 1.5) +
  scale_x_discrete(breaks = paste0("≤", c(seq(14, 61, 3), 61)), expand = c(0, 0.5)) +
  scale_y_continuous(breaks = seq(0, 800, 200)) +
  xlab("Comorbidities per Mendelian disease") +
  ylab("Number of drugs") +
  labs(fill = "") +
  theme_classic() +
  theme(axis.title = element_text(size = 45, color = "black", family = "Arial"),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text.x = element_text(size = 40, color = "black", family = "Arial", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "Arial", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"))
ggarrange(s9_a, s9_b, 
          ncol = 1, nrow = 2, heights = c(3,1), align = "v")
ggsave(filename = "S9_phase_I.png",
       path = "supplementary_figures/", 
       width = 45, height = 30, device = "png",
       dpi = 300, type = type_compression)

# phase2
## For observed OR: plot OR with 95% CI
phase2_obs_or = na.omit(phase2_obs_or) ; rownames(phase2_obs_or) = NULL
phase2_obs_or$comorbidity_filter = paste0("≤", phase2_obs_or$comorbidity_filter)
phase2_obs_or$comorbidity_filter = factor(phase2_obs_or$comorbidity_filter, levels = phase2_obs_or$comorbidity_filter, labels = phase2_obs_or$comorbidity_filter)
phase2_obs_or$status = "observed"
## For permuted ORs: plot 5th, 50th and 95th percentiles
phase2_perm_or_all = na.omit(phase2_perm_or_all) ; rownames(phase2_perm_or_all) = NULL
phase2_perm_or_all$comorbidity_filter = factor(phase2_perm_or_all$comorbidity_filter, levels = phase2_perm_or_all$comorbidity_filter, labels = phase2_perm_or_all$comorbidity_filter)
phase2_perm_or_all_filt = phase2_perm_or_all %>% 
  group_by(comorbidity_filter) %>%
  mutate(percentile_5 = as.numeric(quantile(perm_OR, 0.05)),
         percentile_50 = as.numeric(quantile(perm_OR, 0.50)),
         percentile_95 = as.numeric(quantile(perm_OR, 0.95))) %>%
  filter(perm_OR >= percentile_5 & perm_OR <= percentile_95) %>%
  ungroup() %>%
  dplyr::select(comorbidity_filter, percentile_50, percentile_5,percentile_95) %>%
  distinct() %>%
  mutate(status = "permuted")
# rename columns for rbind later --> will change them after that
colnames(phase2_perm_or_all_filt) = c("comorbidity_filter", "obs_or", "ci95_lower", "ci95_upper", "status")

phase2_obs_perm_combined = rbind(phase2_obs_or[, c(1,2,4,5,6)], phase2_perm_or_all_filt)
colnames(phase2_obs_perm_combined) = c("comorbidity_filter", "OR", "OR_lower_limit", "OR_upper_limit", "status")
phase2_obs_perm_combined$comorbidity_filter = rep(10:61, 2)
phase2_obs_perm_combined$status = factor(phase2_obs_perm_combined$status, levels = phase2_obs_perm_combined$status, labels = phase2_obs_perm_combined$status)

s10_a = ggplot(phase2_obs_perm_combined) +
  geom_ribbon(aes(x = comorbidity_filter, ymin = OR_lower_limit, ymax = OR_upper_limit, fill = status),
              alpha = 0.3) +
  scale_fill_manual(values = c("#00BFC4","#F8766D"))+
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 2) +
  geom_line(aes(x = comorbidity_filter, y = OR, color = status), linewidth = 2, show.legend = FALSE) +
  scale_color_manual(values = c("#00BFC4","#F8766D")) +
  xlab("") +
  scale_x_continuous(breaks = seq(10, 61, 3), labels = paste0("≤", seq(10, 61, 3)), expand = c(0, 0.5)) +
  ylab("Odds Ratio") +
  labs(title = "Phase II", fill = "") +
  theme_classic() +
  theme(plot.title = element_text(size = 70, color = "black", family = "Arial", face = "bold"),
        axis.title = element_text(size = 45, color = "black", family = "Arial"),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text.x = element_text(size = 40, color = "black", family = "Arial", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "Arial", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_text(size = 23, color = "black", family = "Arial"),
        legend.text = element_text(size = 40, family = "Arial"),
        legend.position = "top", legend.key.size = unit(2, "cm"))
s10_b = ggplot(log_reg_observed_nrcomorbidities[2:53, ], aes(x = comorbidity_filter, y = logr_drugs_included)) +
  geom_line(group = "comorbidity_filter", linewidth = 1.5) +
  scale_x_discrete(breaks = paste0("≤", seq(10, 61, 3)), expand = c(0, 0.5)) +
  scale_y_continuous(breaks = seq(0, 800, 200)) +
  xlab("Comorbidities per Mendelian disease") +
  ylab("Number of drugs") +
  labs(fill = "") +
  theme_classic() +
  theme(axis.title = element_text(size = 45, color = "black", family = "Arial"),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text.x = element_text(size = 40, color = "black", family = "Arial", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "Arial", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"))
ggarrange(s10_a, s10_b, 
          ncol = 1, nrow = 2, heights = c(3,1), align = "v")
ggsave(filename = "S10_phase_II.png",
       path = "supplementary_figures/", 
       width = 45, height = 30, device = "png",
       dpi = 300, type = type_compression)

# phase3
## For observed OR: plot OR with 95% CI
phase3_obs_or = na.omit(phase3_obs_or) ; rownames(phase3_obs_or) = NULL
phase3_obs_or$comorbidity_filter = paste0("≤", phase3_obs_or$comorbidity_filter)
phase3_obs_or$comorbidity_filter = factor(phase3_obs_or$comorbidity_filter, levels = phase3_obs_or$comorbidity_filter, labels = phase3_obs_or$comorbidity_filter)
phase3_obs_or$status = "observed"
## For permuted ORs: plot 5th, 50th and 95th percentiles
phase3_perm_or_all = na.omit(phase3_perm_or_all) ; rownames(phase3_perm_or_all) = NULL
phase3_perm_or_all$comorbidity_filter = factor(phase3_perm_or_all$comorbidity_filter, levels = phase3_perm_or_all$comorbidity_filter, labels = phase3_perm_or_all$comorbidity_filter)
phase3_perm_or_all_filt = phase3_perm_or_all %>% 
  group_by(comorbidity_filter) %>%
  mutate(percentile_5 = as.numeric(quantile(perm_OR, 0.05)),
         percentile_50 = as.numeric(quantile(perm_OR, 0.50)),
         percentile_95 = as.numeric(quantile(perm_OR, 0.95))) %>%
  filter(perm_OR >= percentile_5 & perm_OR <= percentile_95) %>%
  ungroup() %>%
  dplyr::select(comorbidity_filter, percentile_50, percentile_5,percentile_95) %>%
  distinct() %>%
  mutate(status = "permuted")
# rename columns for rbind later --> will change them after that
colnames(phase3_perm_or_all_filt) = c("comorbidity_filter", "obs_or", "ci95_lower", "ci95_upper", "status")

phase3_obs_perm_combined = rbind(phase3_obs_or[, c(1,2,4,5,6)], phase3_perm_or_all_filt)
colnames(phase3_obs_perm_combined) = c("comorbidity_filter", "OR", "OR_lower_limit", "OR_upper_limit", "status")
phase3_obs_perm_combined$comorbidity_filter = rep(10:61, 2)
phase3_obs_perm_combined$status = factor(phase3_obs_perm_combined$status, levels = phase3_obs_perm_combined$status, labels = phase3_obs_perm_combined$status)

s11_a = ggplot(phase3_obs_perm_combined) +
  geom_ribbon(aes(x = comorbidity_filter, ymin = OR_lower_limit, ymax = OR_upper_limit, fill = status),
              alpha = 0.3) +
  scale_fill_manual(values = c("#00BFC4","#F8766D"))+
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 2) +
  geom_line(aes(x = comorbidity_filter, y = OR, color = status), linewidth = 2, show.legend = FALSE) +
  scale_color_manual(values = c("#00BFC4","#F8766D")) +
  xlab("") +
  scale_x_continuous(breaks = seq(10, 61, 3), labels = paste0("≤", seq(10, 61, 3)), expand = c(0, 0.5)) +
  scale_y_continuous(breaks = seq(0, 13, 1), limits = c(0, 13)) +
  ylab("Odds Ratio") +
  labs(title = "Phase III", fill = "") +
  theme_classic() +
  theme(plot.title = element_text(size = 70, color = "black", family = "Arial", face = "bold"),
        axis.title = element_text(size = 45, color = "black", family = "Arial"),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text.x = element_text(size = 40, color = "black", family = "Arial", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "Arial", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_text(size = 23, color = "black", family = "Arial"),
        legend.text = element_text(size = 40, family = "Arial"),
        legend.position = "top", legend.key.size = unit(2, "cm"))
s11_b = ggplot(log_reg_observed_nrcomorbidities[2:53, ], aes(x = comorbidity_filter, y = logr_drugs_included)) +
  geom_line(group = "comorbidity_filter", linewidth = 1.5) +
  scale_x_discrete(breaks = paste0("≤", seq(10, 61, 3)), expand = c(0, 0.5)) +
  scale_y_continuous(breaks = seq(0, 800, 200)) +
  xlab("Comorbidities per Mendelian disease") +
  ylab("Number of drugs") +
  labs(fill = "") +
  theme_classic() +
  theme(axis.title = element_text(size = 45, color = "black", family = "Arial"),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text.x = element_text(size = 40, color = "black", family = "Arial", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "Arial", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"))
ggarrange(s11_a, s11_b, 
          ncol = 1, nrow = 2, heights = c(3,1), align = "v")
ggsave(filename = "S11_phase_III.png",
       path = "supplementary_figures/", 
       width = 45, height = 30, device = "png",
       dpi = 300, type = type_compression)

# approved
## For observed OR: plot OR with 95% CI
approved_obs_or = na.omit(approved_obs_or) ; rownames(approved_obs_or) = NULL
approved_obs_or$comorbidity_filter = paste0("≤", approved_obs_or$comorbidity_filter)
approved_obs_or$comorbidity_filter = factor(approved_obs_or$comorbidity_filter, levels = approved_obs_or$comorbidity_filter, labels = approved_obs_or$comorbidity_filter)
approved_obs_or$status = "observed"
## For permuted ORs: plot 5th, 50th and 95th percentiles
approved_perm_or_all = na.omit(approved_perm_or_all) ; rownames(approved_perm_or_all) = NULL
approved_perm_or_all$comorbidity_filter = factor(approved_perm_or_all$comorbidity_filter, levels = approved_perm_or_all$comorbidity_filter, labels = approved_perm_or_all$comorbidity_filter)
approved_perm_or_all_filt = approved_perm_or_all %>% 
  group_by(comorbidity_filter) %>%
  mutate(percentile_5 = as.numeric(quantile(perm_OR, 0.05)),
         percentile_50 = as.numeric(quantile(perm_OR, 0.50)),
         percentile_95 = as.numeric(quantile(perm_OR, 0.95))) %>%
  filter(perm_OR >= percentile_5 & perm_OR <= percentile_95) %>%
  ungroup() %>%
  dplyr::select(comorbidity_filter, percentile_50, percentile_5,percentile_95) %>%
  distinct() %>%
  mutate(status = "permuted")
# rename columns for rbind later --> will change them after that
colnames(approved_perm_or_all_filt) = c("comorbidity_filter", "obs_or", "ci95_lower", "ci95_upper", "status")

approved_obs_perm_combined = rbind(approved_obs_or[, c(1,2,4,5,6)], approved_perm_or_all_filt)
colnames(approved_obs_perm_combined) = c("comorbidity_filter", "OR", "OR_lower_limit", "OR_upper_limit", "status")
approved_obs_perm_combined$comorbidity_filter = rep(10:61, 2)
approved_obs_perm_combined$status = factor(approved_obs_perm_combined$status, levels = approved_obs_perm_combined$status, labels = approved_obs_perm_combined$status)

s12_a = ggplot(approved_obs_perm_combined) +
  geom_ribbon(aes(x = comorbidity_filter, ymin = OR_lower_limit, ymax = OR_upper_limit, fill = status),
              alpha = 0.3) +
  scale_fill_manual(values = c("#00BFC4","#F8766D"))+
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 2) +
  geom_line(aes(x = comorbidity_filter, y = OR, color = status), linewidth = 2, show.legend = FALSE) +
  scale_color_manual(values = c("#00BFC4","#F8766D")) +
  xlab("") +
  scale_x_continuous(breaks = seq(10, 61, 3), labels = paste0("≤", seq(10, 61, 3)), expand = c(0, 0.5)) +
  scale_y_continuous(breaks = seq(0, 3, 0.5), limits = c(0, 3)) +
  ylab("Odds Ratio") +
  labs(title = "Approved", fill = "") +
  theme_classic() +
  theme(plot.title = element_text(size = 70, color = "black", family = "Arial", face = "bold"),
        axis.title = element_text(size = 45, color = "black", family = "Arial"),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text.x = element_text(size = 40, color = "black", family = "Arial", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "Arial", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_text(size = 23, color = "black", family = "Arial"),
        legend.text = element_text(size = 40, family = "Arial"),
        legend.position = "top", legend.key.size = unit(2, "cm"))
s12_b = ggplot(log_reg_observed_nrcomorbidities[2:53, ], aes(x = comorbidity_filter, y = logr_drugs_included)) +
  geom_line(group = "comorbidity_filter", linewidth = 1.5) +
  scale_x_discrete(breaks = paste0("≤", seq(10, 61, 3)), expand = c(0, 0.5)) +
  scale_y_continuous(breaks = seq(0, 800, 200)) +
  xlab("Comorbidities per Mendelian disease") +
  ylab("Number of drugs") +
  labs(fill = "") +
  theme_classic() +
  theme(axis.title = element_text(size = 45, color = "black", family = "Arial"),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text.x = element_text(size = 40, color = "black", family = "Arial", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "Arial", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"))
ggarrange(s12_a, s12_b, 
          ncol = 1, nrow = 2, heights = c(3,1), align = "v")
ggsave(filename = "S12_approved.png",
       path = "supplementary_figures/", 
       width = 45, height = 30, device = "png",
       dpi = 300, type = type_compression)

rm(list = ls())

### comorbidity VS comorbidity + genetic similarity ###

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

md_nr_cancer_comorbidities = md_cancer_comorbidities %>%
  group_by(mendelian_disease) %>%
  mutate(nr_comorbidities = length(complex_disease)) %>%
  ungroup() %>%
  dplyr::select(mendelian_disease, nr_comorbidities) %>%
  distinct()

lr_results_cancers_comorbidity = data.frame(comorbidity_filter = min(md_nr_cancer_comorbidities$nr_comorbidities):max(md_nr_cancer_comorbidities$nr_comorbidities),
                                            or_obs_comorbidity = 0,
                                            or_obs_comorbidity_lower = 0,
                                            or_obs_comorbidity_upper = 0,
                                            or_obs_comorbidity_gensim = 0,
                                            or_obs_comorbidity_gensim_lower = 0,
                                            or_obs_comorbidity_gensim_upper = 0,
                                            nr_drugs = 0)

for (i in 1:nrow(lr_results_cancers_comorbidity)) {
  
  comorbidity_filter = lr_results_cancers_comorbidity[i, "comorbidity_filter"]
  
  md_to_include = md_nr_cancer_comorbidities %>% filter(nr_comorbidities <= comorbidity_filter)
  md_cancer_comorbidities_filt = md_cancer_comorbidities %>% filter(mendelian_disease %in% md_to_include$mendelian_disease)
  
  # drugs targeting MD genes
  drugs_targeting_md_genes = md_drugs %>%
    filter(mendelian_disease %in% md_cancer_comorbidities_filt$mendelian_disease) %>%
    dplyr::select(db_id) %>%
    distinct()
  
  # unique cancers in the data
  unique_cancers_temp = unique(md_cancer_comorbidities_filt$complex_disease)  
  
  # create logistic regression input
  log_input = data.frame(db_id = rep(drugs_targeting_md_genes$db_id, each = length(unique_cancers_temp)))
  
  log_input$complex_disease = unique_cancers_temp
  
  log_input = left_join(log_input, drugs_nr_targets, by = "db_id")  
  
  log_input = left_join(log_input, investigated_indicated_drugs, by = c("complex_disease", "db_id" = "drugbank_id"))
  log_input$indicated_investigated = ifelse(is.na(log_input$indicated_investigated), 0, 1)
  
  cancer_recommended_drugs = left_join(md_cancer_comorbidities_filt, md_drugs, by = "mendelian_disease")
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

  # comorbidity AND genetic similarity
  glm_fits_comorbidity_gensim = glm(indicated_investigated ~ total_targets + comorbidity_gensim,
                                    data = log_input, 
                                    family = binomial())
  log_summary_comorbidity_gensim = summary(glm_fits_comorbidity_gensim)$coefficients

  lr_results_cancers_comorbidity[i, "or_obs_comorbidity"] = exp(log_summary_comorbidity["comorbidity", "Estimate"])
  lr_results_cancers_comorbidity[i, "or_obs_comorbidity_lower"] = exp(confint(glm_fits_comorbidity))["comorbidity", "2.5 %"]
  lr_results_cancers_comorbidity[i, "or_obs_comorbidity_upper"] = exp(confint(glm_fits_comorbidity))["comorbidity", "97.5 %"]
  
  lr_results_cancers_comorbidity[i, "or_obs_comorbidity_gensim"] = exp(log_summary_comorbidity_gensim["comorbidity_gensim", "Estimate"])
  lr_results_cancers_comorbidity[i, "or_obs_comorbidity_gensim_lower"] = exp(confint(glm_fits_comorbidity_gensim))["comorbidity_gensim", "2.5 %"]
  lr_results_cancers_comorbidity[i, "or_obs_comorbidity_gensim_upper"] = exp(confint(glm_fits_comorbidity_gensim))["comorbidity_gensim", "97.5 %"]
  
  lr_results_cancers_comorbidity[i, "nr_drugs"] = length(unique(log_input$db_id))
  
  cat(i, "\n")
}
rm(i, comorbidity_filter, md_to_include, md_cancer_comorbidities_filt, drugs_targeting_md_genes, unique_cancers_temp, log_input, glm_fits_comorbidity, glm_fits_comorbidity_gensim, log_summary_comorbidity, log_summary_comorbidity_gensim)

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

## comorbidity
lr_results_cancers_comorbidity_permutations = data.frame(comorbidity_filter = 0,
                                                         or_comorbidity_perm = 0)
# run permutations
for (i in 1:1000) {
  for (z in 1:10) {
  
    md_to_include = md_nr_cancer_comorbidities %>% filter(nr_comorbidities <= z)
    md_cancer_comorbidities_filt = md_cancer_comorbidities %>% filter(mendelian_disease %in% md_to_include$mendelian_disease)
    
    md_cancer_comorbidities_filt_matrix = md_cancer_comorbidities_filt %>% filter(comorbidity == 1)
    md_cancer_comorbidities_filt_matrix = as.data.frame.matrix(table(md_cancer_comorbidities_filt_matrix[, 1:2]))
    md_cancer_comorbidities_filt_matrix = md_cancer_comorbidities_filt_matrix %>%
      mutate(cancer = rownames(md_cancer_comorbidities_filt_matrix)) %>%
      dplyr::select(cancer, everything())
    rownames(md_cancer_comorbidities_filt_matrix) = NULL

    # drugs targeting MD genes
    drugs_targeting_md_genes_filt = md_drugs %>%
      filter(mendelian_disease %in% md_cancer_comorbidities_filt$mendelian_disease) %>%
      dplyr::select(db_id) %>%
      distinct()
    
    # permutations
    perm_results = permutation(md_cancer_matrix = md_cancer_comorbidities_filt_matrix, 
                               drugs_targeting_md_genes_filt, drugs_nr_targets, investigated_indicated_drugs, md_drugs)
    temp = data.frame(comorbidity_filter = z,
                      or_comorbidity_perm = as.numeric(perm_results[1]))
    lr_results_cancers_comorbidity_permutations = rbind(lr_results_cancers_comorbidity_permutations, temp)
  }
  # track progress
  cat(i, "\n")
}
lr_results_cancers_comorbidity_permutations = lr_results_cancers_comorbidity_permutations[-1, ] ; rownames(lr_results_cancers_comorbidity_permutations) = NULL

## comorbidity + genetic similarity
lr_results_cancers_comorbidity_gensim_permutations = data.frame(comorbidity_filter = 0,
                                                                or_comorbidity_gensim_perm = 0)
# run permutations
for (i in 1:1000) {
  for (z in 1:10) {
    
    md_to_include = md_nr_cancer_comorbidities %>% filter(nr_comorbidities <= z)
    md_cancer_comorbidities_filt = md_cancer_comorbidities %>% filter(mendelian_disease %in% md_to_include$mendelian_disease)
    
    md_cancer_comorbidities_filt_matrix = md_cancer_comorbidities_filt %>% filter(comorbidity == 1 & genetic_similarity_melamed == 1)
    md_cancer_comorbidities_filt_matrix = as.data.frame.matrix(table(md_cancer_comorbidities_filt_matrix[, 1:2]))
    md_cancer_comorbidities_filt_matrix = md_cancer_comorbidities_filt_matrix %>%
      mutate(cancer = rownames(md_cancer_comorbidities_filt_matrix)) %>%
      dplyr::select(cancer, everything())
    rownames(md_cancer_comorbidities_filt_matrix) = NULL
    
    # drugs targeting MD genes
    drugs_targeting_md_genes_filt = md_drugs %>%
      filter(mendelian_disease %in% md_cancer_comorbidities_filt$mendelian_disease) %>%
      dplyr::select(db_id) %>%
      distinct()
    
    # permutations
    perm_results = permutation(md_cancer_matrix = md_cancer_comorbidities_filt_matrix, 
                               drugs_targeting_md_genes_filt, drugs_nr_targets, investigated_indicated_drugs, md_drugs)
    temp = data.frame(comorbidity_filter = z,
                      or_comorbidity_gensim_perm = as.numeric(perm_results[1]))
    lr_results_cancers_comorbidity_gensim_permutations = rbind(lr_results_cancers_comorbidity_gensim_permutations, temp)
  }
  # track progress
  cat(i, "\n")
}
lr_results_cancers_comorbidity_gensim_permutations = lr_results_cancers_comorbidity_gensim_permutations[-1, ] ; rownames(lr_results_cancers_comorbidity_gensim_permutations) = NULL

## visualizations

# comorbidity
## For observed OR: plot OR with 95% CI
lr_results_cancers_comorbidity$comorbidity_filter = paste0("≤", lr_results_cancers_comorbidity$comorbidity_filter)
lr_results_cancers_comorbidity$comorbidity_filter = factor(lr_results_cancers_comorbidity$comorbidity_filter, levels = lr_results_cancers_comorbidity$comorbidity_filter, labels = lr_results_cancers_comorbidity$comorbidity_filter)
lr_results_cancers_comorbidity$status = "observed"
## For permuted ORs: plot 5th, 50th and 95th percentiles
lr_results_cancers_comorbidity_permutations = lr_results_cancers_comorbidity_permutations %>% arrange(comorbidity_filter)
lr_results_cancers_comorbidity_permutations$comorbidity_filter = factor(lr_results_cancers_comorbidity_permutations$comorbidity_filter, levels = lr_results_cancers_comorbidity_permutations$comorbidity_filter, labels = lr_results_cancers_comorbidity_permutations$comorbidity_filter)
lr_results_cancers_comorbidity_permutations_filt = lr_results_cancers_comorbidity_permutations %>% 
  group_by(comorbidity_filter) %>%
  mutate(percentile_5 = as.numeric(quantile(or_comorbidity_perm, 0.05)),
         percentile_50 = as.numeric(quantile(or_comorbidity_perm, 0.50)),
         percentile_95 = as.numeric(quantile(or_comorbidity_perm, 0.95))) %>%
  filter(or_comorbidity_perm >= percentile_5 & or_comorbidity_perm <= percentile_95) %>%
  ungroup() %>%
  dplyr::select(comorbidity_filter, percentile_50, percentile_5,percentile_95) %>%
  distinct() %>%
  mutate(status = "permuted")
# rename columns for rbind later --> will change them after that
colnames(lr_results_cancers_comorbidity_permutations_filt) = c("comorbidity_filter", "or_obs_comorbidity", "or_obs_comorbidity_lower", "or_obs_comorbidity_upper", "status")

comorbidity_obs_perm_combined = rbind(lr_results_cancers_comorbidity[, c(1,2,3,4,8)], lr_results_cancers_comorbidity_permutations_filt)
colnames(comorbidity_obs_perm_combined) = c("comorbidity_filter", "OR", "OR_lower_limit", "OR_upper_limit", "status")
comorbidity_obs_perm_combined$comorbidity_filter = rep(1:10, 2)
comorbidity_obs_perm_combined$status = factor(comorbidity_obs_perm_combined$status, levels = comorbidity_obs_perm_combined$status, labels = comorbidity_obs_perm_combined$status)

s13_a = ggplot(comorbidity_obs_perm_combined) +
  geom_ribbon(aes(x = comorbidity_filter, ymin = OR_lower_limit, ymax = OR_upper_limit, fill = status),
              alpha = 0.3) +
  scale_fill_manual(values = c("#00BFC4","#F8766D"))+
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 2) +
  geom_line(aes(x = comorbidity_filter, y = OR, color = status), linewidth = 2, show.legend = FALSE) +
  scale_color_manual(values = c("#00BFC4","#F8766D")) +
  xlab("") +
  scale_x_continuous(breaks = seq(1, 10, 1), labels = paste0("≤", seq(1, 10, 1)), expand = c(0, 0.5)) +
  scale_y_continuous(breaks = seq(0, 5, 1), limits = c(0, 5)) +
  ylab("Odds Ratio") +
  labs(title = "Comorbidity", fill = "") +
  theme_classic() +
  theme(plot.title = element_text(size = 70, color = "black", family = "Arial", face = "bold"),
        axis.title = element_text(size = 45, color = "black", family = "Arial"),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text.x = element_text(size = 30, color = "black", family = "Arial", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "Arial", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_text(size = 23, color = "black", family = "Arial"),
        legend.text = element_text(size = 40, family = "Arial"),
        legend.position = "top", legend.key.size = unit(2, "cm"))
s13_b = ggplot(lr_results_cancers_comorbidity, aes(x = comorbidity_filter, y = nr_drugs)) +
  geom_line(group = "comorbidity_filter", linewidth = 1.5) +
  scale_x_discrete(expand = c(0, 0.5)) +
  scale_y_continuous(breaks = seq(0, 700, 200), limits = c(0, 700)) +
  xlab("Comorbidities per Mendelian disease") +
  ylab("Number of drugs") +
  labs(fill = "") +
  theme_classic() +
  theme(axis.title = element_text(size = 45, color = "black", family = "Arial"),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text.x = element_text(size = 30, color = "black", family = "Arial", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "Arial", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"))
ggarrange(s13_a, s13_b,
          ncol = 1, nrow = 2, heights = c(2,1), align = "v")
ggsave(filename = "S13_MD_cancers_comorbidity.png",
       path = "supplementary_figures/",
       width = 25, height = 15, device = "png",
       dpi = 300, type = type_compression)

# comorbidity + genetic similarity
## For permuted ORs: plot 5th, 50th and 95th percentiles
lr_results_cancers_comorbidity_gensim_permutations = lr_results_cancers_comorbidity_gensim_permutations %>% arrange(comorbidity_filter)
lr_results_cancers_comorbidity_gensim_permutations$comorbidity_filter = factor(lr_results_cancers_comorbidity_gensim_permutations$comorbidity_filter, levels = lr_results_cancers_comorbidity_gensim_permutations$comorbidity_filter, labels = lr_results_cancers_comorbidity_gensim_permutations$comorbidity_filter)
lr_results_cancers_comorbidity_gensim_permutations_filt = lr_results_cancers_comorbidity_gensim_permutations %>% 
  group_by(comorbidity_filter) %>%
  mutate(percentile_5 = as.numeric(quantile(or_comorbidity_gensim_perm, 0.05)),
         percentile_50 = as.numeric(quantile(or_comorbidity_gensim_perm, 0.50)),
         percentile_95 = as.numeric(quantile(or_comorbidity_gensim_perm, 0.95))) %>%
  filter(or_comorbidity_gensim_perm >= percentile_5 & or_comorbidity_gensim_perm <= percentile_95) %>%
  ungroup() %>%
  dplyr::select(comorbidity_filter, percentile_50, percentile_5,percentile_95) %>%
  distinct() %>%
  mutate(status = "permuted")
# rename columns for rbind later --> will change them after that
colnames(lr_results_cancers_comorbidity_gensim_permutations_filt) = c("comorbidity_filter", "or_obs_comorbidity_gensim", "or_obs_comorbidity_gensim_lower", "or_obs_comorbidity_gensim_upper", "status")

comorbidity_gensim_obs_perm_combined = rbind(lr_results_cancers_comorbidity[, c(1,5,6,7,8)], lr_results_cancers_comorbidity_gensim_permutations_filt)
colnames(comorbidity_gensim_obs_perm_combined) = c("comorbidity_filter", "OR", "OR_lower_limit", "OR_upper_limit", "status")
comorbidity_gensim_obs_perm_combined$comorbidity_filter = rep(1:10, 2)
comorbidity_gensim_obs_perm_combined$status = factor(comorbidity_gensim_obs_perm_combined$status, levels = comorbidity_gensim_obs_perm_combined$status, labels = comorbidity_gensim_obs_perm_combined$status)

s14_a = ggplot(comorbidity_gensim_obs_perm_combined) +
  geom_ribbon(aes(x = comorbidity_filter, ymin = OR_lower_limit, ymax = OR_upper_limit, fill = status),
              alpha = 0.3) +
  scale_fill_manual(values = c("#00BFC4","#F8766D"))+
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", linewidth = 2) +
  geom_line(aes(x = comorbidity_filter, y = OR, color = status), linewidth = 2, show.legend = FALSE) +
  scale_color_manual(values = c("#00BFC4","#F8766D")) +
  xlab("") +
  scale_x_continuous(breaks = seq(1, 10, 1), labels = paste0("≤", seq(1, 10, 1)), expand = c(0, 0.5)) +
  scale_y_continuous(breaks = seq(0, 6, 1), limits = c(0, 6)) +
  ylab("Odds Ratio") +
  labs(title = "Comorbidity + Genetic similarity", fill = "") +
  theme_classic() +
  theme(plot.title = element_text(size = 70, color = "black", family = "Arial", face = "bold"),
        axis.title = element_text(size = 45, color = "black", family = "Arial"),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text.x = element_text(size = 30, color = "black", family = "Arial", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "Arial", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"),
        legend.title = element_text(size = 23, color = "black", family = "Arial"),
        legend.text = element_text(size = 40, family = "Arial"),
        legend.position = "top", legend.key.size = unit(2, "cm"))
s14_b = ggplot(lr_results_cancers_comorbidity, aes(x = comorbidity_filter, y = nr_drugs)) +
  geom_line(group = "comorbidity_filter", linewidth = 1.5) +
  scale_x_discrete(expand = c(0, 0.5)) +
  scale_y_continuous(breaks = seq(0, 700, 200), limits = c(0, 700)) +
  xlab("Comorbidities per Mendelian disease") +
  ylab("Number of drugs") +
  labs(fill = "") +
  theme_classic() +
  theme(axis.title = element_text(size = 45, color = "black", family = "Arial"),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text.x = element_text(size = 30, color = "black", family = "Arial", margin = margin(t = 10)),
        axis.text.y = element_text(size = 40, color = "black", family = "Arial", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 2),
        axis.ticks = element_line(linewidth = 2),
        axis.ticks.length = unit(0.3, "cm"))
ggarrange(s14_a, s14_b,
          ncol = 1, nrow = 2, heights = c(2,1), align = "v")
ggsave(filename = "S14_MD_cancers_comorbidity_genetic_similarity.png",
       path = "supplementary_figures/",
       width = 25, height = 15, device = "png",
       dpi = 300, type = type_compression)
