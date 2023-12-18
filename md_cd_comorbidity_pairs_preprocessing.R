### This script processes the supplementary files from Blair et al. (S3, S4) ###

library(data.table)
library(readxl)
library(dplyr)
library(stringr)
library(reshape2)

## --- Mendelian and complex disease comorbidities --- ##

## load blair et al supplementary file Table 4
## Link: https://www.sciencedirect.com/science/article/pii/S0092867413010246?via%3Dihub
md_cd_comorbidities = read_xls("raw_data/blair_mmc4_md_cd_comorbidities.xls", col_names = TRUE, skip = 1)
md_cd_comorbidities = md_cd_comorbidities %>%
  dplyr::select(complex_disease = "Complex Disease",
                mendelian_disease = "Mendelian Disease") %>%
  mutate(comorbidity = 1) %>%
  arrange(complex_disease, mendelian_disease) %>%
  distinct() 

## remove 5 Mendelian diseases (and their comorbidities) that are due to chromosomal abnormalities --> their causal genes are not obvious
md_cd_comorbidities = md_cd_comorbidities %>%
  filter(!mendelian_disease %in% c("Down Syndrome", "Edward Syndrome", "Klinefelter Syndrome", "Patau Syndrome", "Turner_Syndrome"))
length(unique(md_cd_comorbidities$mendelian_disease)) # 90 unique Mendelian diseases
length(unique(md_cd_comorbidities$complex_disease)) # 65 unique complex diseases

## rename complex diseases for consistency among files
md_cd_comorbidities$complex_disease = gsub("Addison Disease", "Addisons Disease", md_cd_comorbidities$complex_disease, fixed = TRUE)
md_cd_comorbidities$complex_disease = gsub("Alzheimer's Disease", "Alzheimers Disease", md_cd_comorbidities$complex_disease, fixed = TRUE)
md_cd_comorbidities$complex_disease = gsub("Breast Cancer", "Female Breast Cancer", md_cd_comorbidities$complex_disease, fixed = TRUE)
md_cd_comorbidities$complex_disease = gsub("Crohn Disease", "Crohns Disease", md_cd_comorbidities$complex_disease, fixed = TRUE)
md_cd_comorbidities$complex_disease = gsub("Cushing Syndrome", "Cushings Syndrome", md_cd_comorbidities$complex_disease, fixed = TRUE)
md_cd_comorbidities$complex_disease = gsub("Diabetes Mellitus Type II", "Diabetes Mellitus Type 2", md_cd_comorbidities$complex_disease, fixed = TRUE)
md_cd_comorbidities$complex_disease = gsub("Diabetes Mellitus Type I", "Diabetes Mellitus Type 1", md_cd_comorbidities$complex_disease, fixed = TRUE)
md_cd_comorbidities$complex_disease = gsub("Lichen", "Liche", md_cd_comorbidities$complex_disease, fixed = TRUE)
md_cd_comorbidities$complex_disease = gsub("Viral Infection", "Unspecified Viral Infection (Common Cold)", md_cd_comorbidities$complex_disease, fixed = TRUE)

## --- Mendelian disease causal genes --- ##

## load blair et al supplementary file Table 3
## Link: https://www.sciencedirect.com/science/article/pii/S0092867413010246?via%3Dihub
md_genes = read_xls("raw_data/blair_mmc3_md_genes.xls", col_names = TRUE)

## keep only needed information
md_genes = md_genes %>%
  dplyr::select(mendelian_disease = "Summary Name",
                causal_genes = "Genes") %>%
  arrange(mendelian_disease, causal_genes) %>%
  distinct() 

## as above, remove 5 Mendelian diseases that are due to chromosomal abnormalities
md_genes = md_genes %>% 
  filter(!mendelian_disease %in% c("Downs Syndrome", "Edwards Syndrome", "Klinefelters Syndrome", "Pataus Syndrome", "Turners Syndrome"))

## correct Mendelian disease names to match the md_cd_comorbidities file
md_genes$mendelian_disease = gsub("Fragile X Syndrome", "Fragile-X Syndrome", md_genes$mendelian_disease)
md_genes$mendelian_disease = gsub("Bartters Syndrome", "Bartter Syndrome", md_genes$mendelian_disease)
md_genes$mendelian_disease = gsub("Congenital Hirschsprungs Disease", "Congenital Hirschsprung Disease", md_genes$mendelian_disease)
md_genes$mendelian_disease = gsub("Friedreichs Ataxia", "Friedreich Ataxia", md_genes$mendelian_disease)
md_genes$mendelian_disease = gsub("Li-Fraumeni and Related Syndromes", "Li Fraumeni and Related Syndromes", md_genes$mendelian_disease)
md_genes$mendelian_disease = gsub("Spinocerebellar ataxia", "Spinocerebellar Ataxia", md_genes$mendelian_disease)

## each row --> a Mendelian disease - causal gene pair
md_genes = cbind(md_genes[, 1], str_split_fixed(md_genes$causal_genes, ",", n = Inf))
md_genes = reshape2::melt(md_genes, "mendelian_disease", colnames(md_genes)[2:ncol(md_genes)])
md_genes = md_genes %>% 
  dplyr::select(mendelian_disease, causal_gene = value) %>%
  distinct()
md_genes = md_genes[-which(md_genes$causal_gene == ""), ] ; rownames(md_genes) = NULL
md_genes = md_genes %>% arrange(mendelian_disease, causal_gene) %>% distinct()

## save files
dir.create(path = "processed_data")
fwrite(md_cd_comorbidities, "processed_data/md_cd_comorbidities.txt", sep = "\t", row.names = FALSE)
fwrite(md_genes, "processed_data/md_genes.txt", sep = "\t", row.names = FALSE)

## --- visualizations --- ##

## histogram of the number of complex disease comorbidities per Mendelian disease
md_nr_cd_comorbidities = md_cd_comorbidities %>% 
  group_by(mendelian_disease) %>%
  mutate(nr_cd_comorbidities = length(complex_disease)) %>%
  ungroup() %>%
  dplyr::select(mendelian_disease, nr_cd_comorbidities) %>%
  distinct()

fig_1b = ggplot(md_nr_cd_comorbidities, aes(x = nr_cd_comorbidities)) +
  geom_histogram(color = "black", fill = "gray", bins = 61) +
  xlab("Number of comorbities per Mendelian disease") +
  ylab("Count") +
  scale_x_continuous(breaks = seq(0, 60, 5)) +
  scale_y_continuous(breaks = seq(0, 8, 1)) +
  theme_classic() +
  theme(axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.3), 
        axis.title = element_text(angle = 0, hjust = 0.5,
                                  margin = margin(t = 3, unit = "cm"),
                                  size = 16, family = "Arial", colour = "black"),
        axis.text = element_text(angle = 0, hjust = 0.5, vjust = 0.5,
                                 margin = margin(l = 0.5, r = 0.2, unit = "cm"),
                                 size = 16, family = "Arial", colour = "black"),
        legend.title = element_blank())

fig_1b
ggsave(filename = "Fig1B_comorbidities_per_mendelian_disease.tiff", 
       path = "figures/", 
       width = 7, height = 5, device = "tiff",
       dpi = 300, compression = "lzw", type = type_compression)
dev.off()

## number of causal genes per Mendelian disease
md_nr_genes = md_genes %>% 
  group_by(mendelian_disease) %>%
  mutate(nr_genes = length(causal_gene)) %>%
  ungroup() %>%
  dplyr::select(mendelian_disease, nr_genes) %>%
  distinct() %>%
  arrange(nr_genes)
md_nr_genes$mendelian_disease = factor(md_nr_genes$mendelian_disease, levels = md_nr_genes$mendelian_disease, labels = md_nr_genes$mendelian_disease)

ggplot(md_nr_genes, aes(x = nr_genes, y = mendelian_disease)) +
  geom_col(fill = "gray", color = "black") +
  xlab("Number of causal genes") +
  ylab("") +
  scale_x_continuous(breaks = seq(0, 51, 3)) +
  theme(axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.3), 
        axis.title = element_text(angle = 0, hjust = 0.5,
                                  margin = margin(t = 3, unit = "cm"),
                                  size = 30, family = "Arial", colour = "black"),
        axis.text = element_text(angle = 0, hjust = 0.5, vjust = 0.5,
                                 margin = margin(l = 0.5, r = 0.2, unit = "cm"),
                                 size = 30, family = "Arial", colour = "black"))
