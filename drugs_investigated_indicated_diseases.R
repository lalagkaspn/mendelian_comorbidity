### This script finds indicated/investigated drugs for the 65 complex diseases in our sample ###

library(data.table)
library(dplyr)
library(jsonlite)
library(httr)
library(ggplot2)

## To run this script, you need some files that can be downloaded only after you obtain corresponding licenses:
## 1. UMLS - MRCONSO.RFF file
## After you obtain license from UMLS, please downnload the following file:
## https://www.nlm.nih.gov/research/umls/licensedcontent/umlsknowledgesources.html?_gl=1*httt76*_ga*MTU2ODcwNTcxMC4xNjk2NjA3NzQ3*_ga_7147EPK006*MTY5ODY4MDk0Mi4yLjEuMTY5ODY4MDk1MS4wLjAuMA..*_ga_P1FPTH9PL4*MTY5ODY4MDk0Mi4yLjEuMTY5ODY4MDk1MS4wLjAuMA..
## 2. DrugBank - drug information
## After you obtain license from DrugBank, please downnload the following file:
## - Drug targets for approved drugs in DrugBank: https://go.drugbank.com/releases/latest#external-links --> "External Drug Links" --> Drug Group "Approved"
## NOTE: download and place both files in the "raw_data" folder
## Moreover, you need an apiKey for the UMLS database. After obtaining it, define the following variable:
umls_apiKey = "8881fe0b-f54d-44b4-a248-984d6de51c5f"

##################
##              ##
## investigated ##
##              ##
##################

## clinical trials data
## downloaded from the AACT: https://aact.ctti-clinicaltrials.org/
## data of download: November 4, 2022
## due to large file size, it was uploaded to Zenodo
ct = fread("https://zenodo.org/records/10058991/files/AACT_studies.txt?download=1")

# filter for interventional studies
ct_filt = ct %>% filter(study_type == "Interventional")

# filter for interventions that are biological, combination product, drug, genetic
ct_interventions = fread("https://zenodo.org/records/10058991/files/AACT_interventions.txt?download=1")
ct_interventions_filt = ct_interventions %>% filter(intervention_type %in% c("Biological", "Combination Product", "Drug", "Genetic"))
ct_filt = ct_filt %>% filter(nct_id %in% ct_interventions_filt$nct_id)
ct_filt = ct_filt %>% dplyr::select(nct_id, overall_status, last_known_status, phase)

# filter for interventions and conditions with matched mesh terms
ct_interventions_mesh = fread("https://zenodo.org/records/10058991/files/AACT_browse_interventions.txt?download=1")
ct_conditions_mesh = fread("https://zenodo.org/records/10058991/files/AACT_browse_conditions.txt?download=1")
ct_to_keep = intersect(ct_interventions_mesh$nct_id, ct_conditions_mesh$nct_id)
ct_filt = ct_filt %>% filter(nct_id %in% ct_to_keep)
ct_interventions_mesh = ct_interventions_mesh %>% filter(nct_id %in% ct_to_keep) %>% dplyr::select(-id)
ct_conditions_mesh = ct_conditions_mesh %>% filter(nct_id %in% ct_to_keep) %>% dplyr::select(-id)
rm(ct, ct_interventions, ct_interventions_filt, ct_to_keep)

# filter for mesh-list
ct_conditions_mesh = ct_conditions_mesh %>% filter(mesh_type == "mesh-list") %>% dplyr::select(nct_id, downcase_mesh_term) %>% distinct()
ct_interventions_mesh = ct_interventions_mesh %>% filter(mesh_type == "mesh-list") %>% dplyr::select(nct_id, downcase_mesh_term) %>% distinct()
ct_filt = left_join(ct_filt, ct_conditions_mesh, by = "nct_id")
ct_filt = left_join(ct_filt, ct_interventions_mesh, by = "nct_id")
ct_filt = ct_filt %>% dplyr::rename(condition_downcase_mesh = downcase_mesh_term.x, intervention_downcase_mesh = downcase_mesh_term.y)

ct_conditions_unique = ct_conditions_mesh %>% 
  filter(nct_id %in% ct_filt$nct_id) %>%
  dplyr::select(downcase_mesh_term) %>% 
  distinct() %>% 
  arrange(downcase_mesh_term)
rm(ct_conditions_mesh)

ct_interventions_unique = ct_interventions_mesh %>% 
  filter(nct_id %in% ct_filt$nct_id) %>%
  dplyr::select(downcase_mesh_term) %>% 
  distinct() %>% 
  arrange(downcase_mesh_term)
rm(ct_interventions_mesh)

length(unique(ct_filt$nct_id)) # 109,430 unique clinical trials remained
length(unique(ct_conditions_unique$downcase_mesh_term)) # 3,267 unique conditions (diseases)
length(unique(ct_interventions_unique$downcase_mesh_term)) # 3,376 unique interventions (drugs)

## conditions: add MeSH codes
# load MRCONSO.RRF (UMLS)
rrf = data.table::fread("raw_data/MRCONSO.RRF", sep = "|", quote = "")
rrf = rrf[, -19]
# column names: https://www.ncbi.nlm.nih.gov/books/NBK9685/table/ch03.T.concept_names_and_sources_file_mr/?report=objectonly
colnames(rrf) = c("CUI", "LAT", "TS", "LUI", "STT", "SUI", "ISPREF", "AUI", "SAUI", "SCUI", "SDUI", "SAB", "TTY", "CODE" , "STR", "SRL", "SUPPRESS", "CVF")
# filter for english terms only
rrf = rrf %>% filter(LAT == "ENG" & SAB == "MSH")
rrf$STR_lowercase = tolower(rrf$STR)
# add mesh_codes
ct_conditions_unique = left_join(ct_conditions_unique, rrf[, c("STR_lowercase", "CODE")], by = c("downcase_mesh_term" = "STR_lowercase"))
ct_conditions_unique = ct_conditions_unique %>% dplyr::rename(mesh_code = CODE)
# remove subheading mesh codes (starting with Q)
ct_conditions_unique = ct_conditions_unique[-which(ct_conditions_unique$mesh_code == "Q000506"), ]
rownames(ct_conditions_unique) = NULL

## filter clinical trials for conditions that are among the 65 complex diseases in Blair et al (using mesh_codes)
# blair complex diseases - MeSh terms were manually matched to their ICD10 codes
complex_disease_mesh = fread("raw_data/cd_icd10_mesh_manually_matched.txt")

# filter clinical trial conditions for Blair complex diseases
ct_conditions_unique_blair = ct_conditions_unique %>% filter(mesh_code %in% complex_disease_mesh$mesh_code)

# filter clinical trials for Blair complex diseases
ct_filt_blair = ct_filt %>% filter(condition_downcase_mesh %in% ct_conditions_unique_blair$downcase_mesh_term) %>% distinct()

# match clinical trial conditions to MeSH codes and disease categories
ct_filt_blair = left_join(ct_filt_blair, ct_conditions_unique_blair, by = c("condition_downcase_mesh" = "downcase_mesh_term"))
ct_filt_blair = left_join(ct_filt_blair, complex_disease_mesh[, c(1,3)], by = "mesh_code")

# add complex disease categories
complex_disease_categories = fread("raw_data/complex_disease_category.txt")
ct_filt_blair = left_join(ct_filt_blair, complex_disease_categories, by = "complex_disease") ; rm(complex_disease_categories)
ct_filt_blair = ct_filt_blair %>% dplyr::select(nct_id:phase, complex_disease, mesh_code, disease_category, intervention_downcase_mesh) %>% distinct()

length(unique(ct_filt_blair$nct_id)) # 36,624 clinical trials
length(unique(ct_filt_blair$complex_disease)) # 65 complex diseases
length(unique(ct_filt_blair$intervention_downcase_mesh)) # 2,362 interventions (drugs)

## interventions: add MeSH codes
ct_interventions_unique_blair = ct_filt_blair %>% dplyr::select(intervention_downcase_mesh) %>% distinct()
ct_interventions_unique_blair = left_join(ct_interventions_unique_blair, rrf[, c("STR_lowercase", "CODE")], by = c("intervention_downcase_mesh" = "STR_lowercase"))
ct_interventions_unique_blair = ct_interventions_unique_blair %>% dplyr::rename(mesh_code = CODE) %>% distinct()

## interventions: MeSH codes --> DrugBank IDs
## UMLS API crosswalk function
pageSize = 10000
# empty data frame to populate later with DrugBank IDs
ct_interventions_unique_blair_mesh_db = data.frame(downcase_mesh = NA, mesh_code = NA, db_id = NA)
for (i in 1:nrow(ct_interventions_unique_blair)) {
  
  # prepare UMLS API call
  source_vobaculary = "MSH"
  target_vocabulary = "DRUGBANK"
  mesh_code = ct_interventions_unique_blair[i, mesh_code]
  call = paste0("https://uts-ws.nlm.nih.gov/rest/crosswalk/current/source/", source_vobaculary, "/", mesh_code, "?targetSource=", target_vocabulary ,"&apiKey=", umls_apiKey, "&pageSize=", pageSize)
  
  # call UMLS API
  res = GET(call)
  
  # res is a JSON output
  data = fromJSON(rawToChar(res$content))
  
  # get drugbank_ids
  db_id = data$result$ui

  # if no drugbank_ids returned
  if (length(db_id) == 0) {
    no_match = data.frame(downcase_mesh = ct_interventions_unique_blair[i, intervention_downcase_mesh], mesh_code = mesh_code, db_id = NA)
    ct_interventions_unique_blair_mesh_db = rbind(ct_interventions_unique_blair_mesh_db, no_match)
  }

  # if only 1 drugbank_ids returned
  if (length(db_id) == 1) {
    match = data.frame(downcase_mesh = ct_interventions_unique_blair[i, intervention_downcase_mesh], mesh_code = mesh_code, db_id = db_id)
    ct_interventions_unique_blair_mesh_db = rbind(ct_interventions_unique_blair_mesh_db, match)
  }

  # if >1 drugbank_ids returned
  if (length(db_id) > 1) {
    for (z in 1:length(db_id)) {
      match = data.frame(downcase_mesh = ct_interventions_unique_blair[i, intervention_downcase_mesh], mesh_code = mesh_code, db_id = db_id[z])
      ct_interventions_unique_blair_mesh_db = rbind(ct_interventions_unique_blair_mesh_db, match)
    }
  }

  # track progress
  cat(i, "/", nrow(ct_interventions_unique_blair), "\n")
}
rm(data, match, no_match, res, call, db_id, i, mesh_code, pageSize, source_vobaculary, target_vocabulary, z)
ct_interventions_unique_blair_mesh_db = ct_interventions_unique_blair_mesh_db[-1, ] ; rownames(ct_interventions_unique_blair_mesh_db) = NULL
# interventions that didn't map to drugbank_id
ct_interventions_no_db = ct_interventions_unique_blair_mesh_db[which(is.na(ct_interventions_unique_blair_mesh_db$db_id)), ]
ct_interventions_no_db = ct_interventions_no_db %>% dplyr::select(-db_id) %>% distinct()
length(unique(ct_interventions_no_db$mesh_code)) # 569 mesh_codes not matched to db_ids
# remove unmatched
ct_interventions_unique_blair_mesh_db = ct_interventions_unique_blair_mesh_db[-which(is.na(ct_interventions_unique_blair_mesh_db$db_id)), ] ; rownames(ct_interventions_unique_blair_mesh_db) = NULL
length(unique(ct_interventions_unique_blair_mesh_db$mesh_code)) # 1,793 mesh_codes matched to db_ids

## map unmatched MeSh codes to RxNORM CUIs
unmatched_mesh_to_rxcui = data.frame(downcase_mesh = NA, mesh_code = NA, rxcui = NA)
pageSize = 10000
for (i in 1:nrow(ct_interventions_no_db)) {
 
  # prepare UMLS API call
  source_vobaculary = "MSH"
  target_vocabulary = "RXNORM"
  mesh_code = ct_interventions_no_db[i, "mesh_code"]
  call = paste0("https://uts-ws.nlm.nih.gov/rest/crosswalk/current/source/", source_vobaculary, "/", mesh_code, "?targetSource=", target_vocabulary ,"&apiKey=", umls_apiKey, "&pageSize=", pageSize)
 
  # call UMLS API
  res = GET(call)
  
  # res is a JSON output
  data = fromJSON(rawToChar(res$content))
  
  # Get rxcuis
  rxcui = data$result$ui

  # if no rxcui returned
  if (length(rxcui) == 0) {
    no_match = data.frame(downcase_mesh = ct_interventions_no_db[i, "downcase_mesh"], mesh_code = mesh_code, rxcui = NA)
    unmatched_mesh_to_rxcui = rbind(unmatched_mesh_to_rxcui, no_match)
  }

  # if only 1 rxcui returned
  if (length(rxcui) == 1) {
    match = data.frame(downcase_mesh = ct_interventions_no_db[i, "downcase_mesh"], mesh_code = mesh_code, rxcui = rxcui)
    unmatched_mesh_to_rxcui = rbind(unmatched_mesh_to_rxcui, match)
  }

  # if >1 rxcui returned
  if (length(rxcui) > 1) {
    for (z in 1:length(rxcui)) {
      match = data.frame(downcase_mesh = ct_interventions_no_db[i, "downcase_mesh"], mesh_code = mesh_code, rxcui = rxcui[z])
      unmatched_mesh_to_rxcui = rbind(unmatched_mesh_to_rxcui, match)
    }
  }

  # track progress
  cat(i, "/", nrow(ct_interventions_no_db), "\n")
}
rm(call, i, mesh_code, pageSize, rxcui, source_vobaculary, target_vocabulary, z, data, res, match, no_match)
unmatched_mesh_to_rxcui = unmatched_mesh_to_rxcui[-1, ]

# find those that didn't match to RxCUIs
unmatched_mesh_to_rxcui_no_match = unmatched_mesh_to_rxcui[which(is.na(unmatched_mesh_to_rxcui$rxcui)), ]
length(unique(unmatched_mesh_to_rxcui_no_match$mesh_code)) # 398 mesh_codes didn't match to rxcuis
rm(unmatched_mesh_to_rxcui_no_match)

# remove unmatched
unmatched_mesh_to_rxcui = na.omit(unmatched_mesh_to_rxcui) ; rownames(unmatched_mesh_to_rxcui) = NULL
length(unique(unmatched_mesh_to_rxcui$mesh_code)) # 171 mesh_codes matched to rxcuis

## The MeSh codes that matched to rxcuis might be:
# - Drug combinations --> will split to the corresponding ingredients ("part_of")
# - Specific forms of drugs (e.g., liposomal doxorubicin) --> will match to ingredient name ("form_of")
# get for each rxcui the "part_of" and "form_of"
pageSize = 10000
rxnorm_expanded = data.frame(downcase_mesh = NA, mesh_code = NA, rxcui = NA, relation = NA, rxcui_new = NA, rx_name = NA)
for (i in 1:nrow(unmatched_mesh_to_rxcui)) {
  
  rxcui = unmatched_mesh_to_rxcui[i, "rxcui"]
  call = paste0("https://uts-ws.nlm.nih.gov/rest/content/current/source/RXNORM/", rxcui, "/relations?includeAdditionalRelationLabels=has_part,form_of&apiKey=", umls_apiKey)
  
  # call UMLS API
  res = GET(call)
  
  # res is a JSON output
  data = fromJSON(rawToChar(res$content))

  # match
  if (length(data) == 4) {
    relation = data[["result"]][["additionalRelationLabel"]]
    rxcui_new = data[["result"]][["relatedId"]]
    rx_name = data[["result"]][["relatedIdName"]]
    
    if (length(relation) == 1) {
      y = data.frame(downcase_mesh = unmatched_mesh_to_rxcui[i, "downcase_mesh"],
                     mesh_code = unmatched_mesh_to_rxcui[i, "mesh_code"],
                     rxcui = unmatched_mesh_to_rxcui[i, "rxcui"],
                     relation = relation,
                     rxcui_new = basename(rxcui_new),
                     rx_name = rx_name)
      rxnorm_expanded = rbind(rxnorm_expanded, y)
    }
    
    if (length(relation) > 1) {
      for (a in 1:length(relation)) {
        y = data.frame(downcase_mesh = unmatched_mesh_to_rxcui[i, "downcase_mesh"],
                       mesh_code = unmatched_mesh_to_rxcui[i, "mesh_code"],
                       rxcui = unmatched_mesh_to_rxcui[i, "rxcui"],
                       relation = relation[a],
                       rxcui_new = basename(rxcui_new[a]),
                       rx_name = rx_name[a])
        rxnorm_expanded = rbind(rxnorm_expanded, y)
      }
    }
  }
  
  # no match
  if (length(data) == 3) {
    y = data.frame(downcase_mesh = unmatched_mesh_to_rxcui[i, "downcase_mesh"],
                   mesh_code = unmatched_mesh_to_rxcui[i, "mesh_code"],
                   rxcui = unmatched_mesh_to_rxcui[i, "rxcui"],
                   relation = NA,
                   rxcui_new = NA,
                   rx_name = NA)
    rxnorm_expanded = rbind(rxnorm_expanded, y)
  }

  cat(i, "/", nrow(unmatched_mesh_to_rxcui), "\n")
}
rm(call, i, pageSize, rxcui, a, data, res, y, relation, rx_name, rxcui_new)
rxnorm_expanded = rxnorm_expanded[-1, ]
rownames(rxnorm_expanded) = NULL

## not matched
rxnorm_not_expanded = rxnorm_expanded[which(is.na(rxnorm_expanded$rxcui_new)), ]
rxnorm_not_expanded = rxnorm_not_expanded %>% dplyr::select(downcase_mesh, mesh_code) %>% distinct()

## matched
rxnorm_expanded = na.omit(rxnorm_expanded)
rownames(rxnorm_expanded) = NULL
rxnorm_not_expanded = rxnorm_not_expanded %>% filter(!mesh_code %in% rxnorm_expanded$mesh_code) %>% distinct()
length(unique(rxnorm_not_expanded$mesh_code)) # 75 mesh_codes matched to rxcuis didn't match to form_of or part_of
length(unique(rxnorm_expanded$mesh_code)) # 96 mesh_codes matched to rxcuis matched to form_of or part_of
rm(rxnorm_not_expanded)

## match RxCUIs to DB_IDs
pageSize = 10000
rxnorm_expanded$db_id = "NA"
for (i in 1:nrow(rxnorm_expanded)) {
  
  # prepare UMLS API call
  source_vobaculary = "RXNORM"
  target_vocabulary = "DRUGBANK"
  rxcui = rxnorm_expanded[i, "rxcui_new"]
  call = paste0("https://uts-ws.nlm.nih.gov/rest/crosswalk/current/source/", source_vobaculary, "/", rxcui, "?targetSource=", target_vocabulary ,"&apiKey=", umls_apiKey, "&pageSize=", pageSize)
  
  # call UMLS API
  res = GET(call)
  
  # res is a JSON output
  data = fromJSON(rawToChar(res$content))
  
  # get drugbank_ids
  db_id = data$result$ui
  
  # if no db_id returned
  if (length(db_id) == 0) {
    next
  }
  
  # if only 1 db_id returned
  if (length(db_id) == 1) { 
    rxnorm_expanded[i, "db_id"] = db_id
  }
  
  # if >1 db_id returned
  if (length(db_id) > 1) {
    rxnorm_expanded[i, "db_id"] = paste0(db_id, collapse = ",")
  }
  
  # track progress
  cat(i, "/", nrow(rxnorm_expanded), "\n")
}
rm(call, i, pageSize, db_id, data, res, rxcui, source_vobaculary, target_vocabulary)

# remove not matched
rxnorm_expanded_no_dbid = rxnorm_expanded[which(rxnorm_expanded$db_id == "NA"), ]
rxnorm_expanded_no_dbid = unique(rxnorm_expanded_no_dbid[, 1:2])
rxnorm_expanded = rxnorm_expanded[-which(rxnorm_expanded$db_id == "NA"), ]
rxnorm_expanded_no_dbid = rxnorm_expanded_no_dbid %>% filter(!mesh_code %in% rxnorm_expanded$mesh_code)
length(unique(rxnorm_expanded_no_dbid$downcase_mesh)) # 14 mesh_codes matched to "part_of", "form_of" but didn't match to drugbank_ids
length(unique(rxnorm_expanded$mesh_code)) # 82 mesh codes matched to "part_of", "form_of" and then to drugbank_ids
rm(rxnorm_expanded_no_dbid)

## merge files
length(unique(ct_interventions_unique_blair_mesh_db$mesh_code)) # 1,793 mesh_codes matched to db_id by MeSH --> DB_IDs
length(unique(rxnorm_expanded$mesh_code)) # 82 mesh_codes matched to db_id by MeSH --> RxCUIs --> "part_of", "form_of" --> DB_IDs

rxnorm_expanded = rxnorm_expanded %>% dplyr::select(downcase_mesh, mesh_code, db_id) %>% distinct()
# split rows that have more than one db_id
x = stringr::str_split_fixed(rxnorm_expanded$db_id, ",", n = Inf)
rxnorm_expanded = rxnorm_expanded[, 1:2]
rxnorm_expanded = cbind(rxnorm_expanded, x)
rxnorm_expanded = reshape2::melt(rxnorm_expanded, c("downcase_mesh", "mesh_code"), colnames(rxnorm_expanded)[3:4])
rm(x)
rxnorm_expanded = rxnorm_expanded[, -3]
rxnorm_expanded = rxnorm_expanded[-which(rxnorm_expanded$value == ""), ]
colnames(rxnorm_expanded) = c("downcase_mesh", "mesh_code", "db_id")
rxnorm_expanded = unique(rxnorm_expanded)
rownames(rxnorm_expanded) = NULL

ct_interventions_unique_blair_mesh_db = rbind(ct_interventions_unique_blair_mesh_db, rxnorm_expanded)
rm(rxnorm_expanded)
ct_interventions_unique_blair_mesh_db = unique(ct_interventions_unique_blair_mesh_db)
length(unique(ct_interventions_unique_blair_mesh_db$mesh_code)) # 1,875 mesh_codes matched to drugbank_ids

## update clinical trial files with DB_IDs for the interventions
ct_filt_blair = left_join(ct_filt_blair, ct_interventions_unique_blair_mesh_db[, c(1,3)], by = c("intervention_downcase_mesh" = "downcase_mesh"))

## remove clinical trials without DB_IDs matched
ct_filt_blair = na.omit(ct_filt_blair) ; rownames(ct_filt_blair) = NULL

length(unique(ct_filt_blair$nct_id)) # 35,156 clinical trials
length(unique(ct_filt_blair$complex_disease)) # 65 complex diseases
length(unique(ct_filt_blair$intervention_downcase_mesh)) # 1,875 inteventions (drugs)
rm(ct_conditions_unique, ct_conditions_unique_blair, ct_filt, ct_interventions_unique, ct_interventions_unique_blair, ct_interventions_unique_blair_mesh_db, rrf)

## Rename clinical trial phases
ct_filt_blair[which(ct_filt_blair$phase == "Not Applicable"), "phase"] = 0
ct_filt_blair[which(ct_filt_blair$phase == "Early Phase 1"), "phase"] = 1
ct_filt_blair[which(ct_filt_blair$phase == "Phase 1"), "phase"] = 1
ct_filt_blair[which(ct_filt_blair$phase == "Phase 1/Phase 2"), "phase"] = 1.5
ct_filt_blair[which(ct_filt_blair$phase == "Phase 2"), "phase"] = 2
ct_filt_blair[which(ct_filt_blair$phase == "Phase 2/Phase 3"), "phase"] = 2.5
ct_filt_blair[which(ct_filt_blair$phase == "Phase 3"), "phase"] = 3
ct_filt_blair[which(ct_filt_blair$phase == "Phase 4"), "phase"] = 4
table(ct_filt_blair$phase)

# save file
fwrite(ct_filt_blair, "processed_data/ct_blair.txt", sep = "\t", row.names = FALSE)

###############
##           ##
## indicated ##
##           ##
###############

## load approved DrugBank drugs
db_drugs_approved = data.table::fread("raw_data/drug_links_approved.csv")
db_drugs_approved = db_drugs_approved %>% dplyr::select(drugbank_id = "DrugBank ID", drug_name = Name) %>% distinct()

## RxNorm ##

# annotate approved drugs with indications from RxNORM
rxnorm_drugbank_indications = data.frame(rela = NA, relaSource = NA, minConcept.rxcui = NA, minConcept.name = NA, minConcept.tty = NA, rxclassMinConceptItem.classId = NA, rxclassMinConceptItem.className = NA, rxclassMinConceptItem.classType = NA, drugbank_id = NA)
# replace all "%" with "_" in order to avoid conflicts with RxNorm API
db_drugs_approved$rxnorm_api = gsub("%", "_", db_drugs_approved$drug_name)
# for each drug: add all the information stored in RxNorm (including indications)
for (i in 1:nrow(db_drugs_approved)) {
  
  # prepare API call
  drug_name = db_drugs_approved[i, rxnorm_api]
  drug_db_id = db_drugs_approved[i, drugbank_id]
  drug_name_http = gsub(" ", "%20", drug_name)
  call = paste0("https://rxnav.nlm.nih.gov/REST/rxclass/class/byDrugName?drugName=", drug_name_http)
  
  # call API
  get_cui = GET(url = call)
  
  # convert content to text
  get_cui = content(get_cui, "text", encoding = "UTF-8")
  
  # parse data in JSON
  get_cui_json = fromJSON(get_cui, flatten = TRUE)
  
  # convert into dataframe
  get_cui_json_dataframe = get_cui_json[["rxclassDrugInfoList"]][["rxclassDrugInfo"]]
  
  # if no match
  if (is.null(get_cui_json_dataframe) == TRUE) {
    get_cui_json_dataframe = get_cui_json_dataframe[, 1:8]
    null = data.frame(rela = NA, relaSource = NA, minConcept.rxcui = NA, minConcept.name = drug_name, minConcept.tty = NA, rxclassMinConceptItem.classId = NA,
                      rxclassMinConceptItem.className = NA, rxclassMinConceptItem.classType = NA, drugbank_id = NA)
    rxnorm_drugbank_indications = rbind(rxnorm_drugbank_indications, null)
  }
  
  # if match
  if (is.null(get_cui_json_dataframe) == FALSE) {
    get_cui_json_dataframe = get_cui_json_dataframe[, 1:8]
    get_cui_json_dataframe$drugbank_id = drug_db_id
    rxnorm_drugbank_indications = rbind(rxnorm_drugbank_indications, get_cui_json_dataframe)
  }
  
  # trac progress
  cat(i, "/", nrow(db_drugs_approved), "\n")
}
rm(drug_name, drug_db_id, drug_name_http, call, get_cui, i, get_cui_json, get_cui_json_dataframe, null)
rxnorm_drugbank_indications = rxnorm_drugbank_indications[-1, ] ; rownames(rxnorm_drugbank_indications) = NULL

# keep only INGREDIENTS and DISEASE
# https://lhncbc.nlm.nih.gov/RxNav/applications/RxClassIntro.html
# https://lhncbc.nlm.nih.gov/RxNav/applications/RxNavViews.html#label:appendix
rxnorm_drugbank_indications_ingredients_disease = rxnorm_drugbank_indications %>% 
  filter(minConcept.tty == "IN" & rxclassMinConceptItem.classType == "DISEASE" & rela %in% c("may_treat", "may_prevent")) %>%
  dplyr::select(drugbank_id, indication = rxclassMinConceptItem.className, relaSource) %>%
  distinct()
rm(rxnorm_drugbank_indications)

# add MeSH codes
# UMLS MRCONSO.RRF file
rrf = data.table::fread("raw_data/MRCONSO.RRF", sep = "|", quote = "")
# NA column -> remove
rrf = rrf[, -19]
# Column names according to: https://www.ncbi.nlm.nih.gov/books/NBK9685/table/ch03.T.concept_names_and_sources_file_mr/?report=objectonly
colnames(rrf) = c("CUI", "LAT", "TS", "LUI", "STT", "SUI", "ISPREF", "AUI", "SAUI", "SCUI", "SDUI", "SAB", "TTY", "CODE" , "STR", "SRL", "SUPPRESS", "CVF")
rrf = rrf %>% filter(LAT == "ENG")

rxnorm_indications_mesh_codes = rrf %>%
  filter(SAB == "MSH" & STR %in% rxnorm_drugbank_indications_ingredients_disease$indication) %>%
  dplyr::select(CUI, SAB, CODE, STR) %>% 
  distinct()
rxnorm_drugbank_indications_ingredients_disease = left_join(rxnorm_drugbank_indications_ingredients_disease, rxnorm_indications_mesh_codes, by = c("indication" = "STR"))
rxnorm_drugbank_indications_ingredients_disease = rxnorm_drugbank_indications_ingredients_disease %>%
  dplyr::select(drugbank_id, mesh_indication = indication, mesh_code = CODE) %>%
  distinct()

## repoDB ##

## downloaded from: https://unmtid-shinyapps.net/shiny/repodb/

# load repoDB data
repodb = data.table::fread("raw_data/repoDB_Drugs_indications_2017.csv")
repodb_approved = repodb %>% 
  filter(status == "Approved") %>% 
  dplyr::select(db_id = drug_id, indication = ind_name, indication_cui = ind_id) %>%
  mutate(source = "UMLS_CUI") %>%
  distinct()
rm(repodb)

## 8 of the CUIs in repodb_approved were merged into other CUIs in later UMLS versions (2022AB)
# C0334634 (2022AA) --> C4721414
repodb_approved$indication_cui = gsub(pattern = "C0334634", replacement = "C4721414", repodb_approved$indication_cui)
# C0178299 (2018AA) --> C0037278
repodb_approved$indication_cui = gsub(pattern = "C0178299", replacement = "C0037278", repodb_approved$indication_cui)
# C1279247 (2020AA) --> C1279264
repodb_approved$indication_cui = gsub(pattern = "C1279247", replacement = "C1279264", repodb_approved$indication_cui)
# C1960397 (2018AA) --> C1709527
repodb_approved$indication_cui = gsub(pattern = "C1960397", replacement = "C1709527", repodb_approved$indication_cui)
# C0002440 (2018AB) --> C0152500
repodb_approved$indication_cui = gsub(pattern = "C0002440", replacement = "C0152500", repodb_approved$indication_cui)
# C0342826 (2022AA) --> C1843116
repodb_approved$indication_cui = gsub(pattern = "C0342826", replacement = "C1843116", repodb_approved$indication_cui)
# C0342827 (2022AA) --> C1856127
repodb_approved$indication_cui = gsub(pattern = "C0342827", replacement = "C1856127", repodb_approved$indication_cui)
# C2584877 (2021AA) --> C2750514
repodb_approved$indication_cui = gsub(pattern = "C2584877", replacement = "C2750514", repodb_approved$indication_cui)

# convert UMLS_CUIs to MeSH
umls_cuis_to_mesh = rrf %>% 
  filter(CUI %in% repodb_approved$indication_cui & SAB == "MSH") %>%
  dplyr::select(CUI, SAB, CODE, STR) %>%
  distinct()
repodb_approved = left_join(repodb_approved, umls_cuis_to_mesh, by = c("indication_cui" = "CUI"))

# remove unmatched umls_cuis
repodb_approved = na.omit(repodb_approved)
repodb_approved = repodb_approved %>% dplyr::select(drugbank_id = db_id, mesh_indication = STR, mesh_code = CODE)
rm(umls_cuis_to_mesh)

# filter for DrugBank approved drugs
repodb_approved = repodb_approved %>% filter(drugbank_id %in% db_drugs_approved$drugbank_id)
length(unique(repodb_approved$drugbank_id)) # 1,321 drugs
length(unique(repodb_approved$mesh_code)) # 734 indications

## merge RxNORM + repoDB ##
db_drugs_indications_mesh = rbind(rxnorm_drugbank_indications_ingredients_disease, repodb_approved)
rm(rxnorm_drugbank_indications_ingredients_disease, repodb_approved)
db_drugs_indications_mesh = unique(db_drugs_indications_mesh)

length(unique(db_drugs_indications_mesh$drugbank_id)) # 2,382 drugs
length(unique(db_drugs_indications_mesh$mesh_code)) # 1,411 indications


## Filter for Blair complex diseases ##

db_drugs_indications_blair = db_drugs_indications_mesh %>% filter(mesh_code %in% complex_disease_mesh$mesh_code) %>% distinct()
db_drugs_indications_blair = left_join(db_drugs_indications_blair, complex_disease_mesh, by = "mesh_code")
# add disease category
complex_disease_categories = fread("raw_data/complex_disease_category.txt")
db_drugs_indications_blair = left_join(db_drugs_indications_blair, complex_disease_categories, by = "complex_disease")
db_drugs_indications_blair = db_drugs_indications_blair %>% dplyr::select(drugbank_id, complex_disease, disease_category)
length(unique(db_drugs_indications_blair$drugbank_id)) # 939 drugs
length(unique(db_drugs_indications_blair$complex_disease)) # 58 complex diseases

db_drugs_indications_blair = unique(db_drugs_indications_blair)

# complex diseases with no indicated drugs
setdiff(complex_disease_mesh$complex_disease, db_drugs_indications_blair$complex_disease)
# Cataract
# General Frontotemporal Dementia
# Hypotony of the Eye
# Kawasaki Disease
# Picks Disease
# Purulent Endophthalmitis
# Skin Cyst

rm(db_drugs_approved, db_drugs_indications_mesh, rrf, rxnorm_indications_mesh_codes)

#######################
##                   ##
## ultimate matrices ##
##                   ##
#######################

## add disease categories to complex_disease_mesh data frame
complex_disease_mesh = left_join(complex_disease_mesh, complex_disease_categories, by = "complex_disease") ; rm(complex_disease_categories)

## Create ultimate matrix ##

# indicated/investigated drug vector
drugs_indicated_investigated = unique(c(ct_filt_blair$db_id, db_drugs_indications_blair$drugbank_id)) # 2,155 drugs
drugs_indicated_investigated = sort(drugs_indicated_investigated)

# complex diseases
complex_diseases_vector = unique(complex_disease_mesh$complex_disease)
complex_diseases_vector = sort(complex_diseases_vector)

ultimate_per_disease = matrix(ncol = length(complex_diseases_vector),
                              nrow = length(drugs_indicated_investigated))
colnames(ultimate_per_disease) = complex_diseases_vector
rownames(ultimate_per_disease) = drugs_indicated_investigated

for (i in 1:nrow(ultimate_per_disease)) {

  # drug
  drug = rownames(ultimate_per_disease)[i]

  # if drug indicated for a complex disease
  if (sum(drug %in% db_drugs_indications_blair$drugbank_id) > 0) {
    drug_conditions_indicated = db_drugs_indications_blair %>% filter(drugbank_id == drug) %>% dplyr::select(complex_disease) %>% distinct()
    for (z in 1:nrow(drug_conditions_indicated)) {
      ultimate_per_disease[i, drug_conditions_indicated[z, "complex_disease"]] = 5
    }
  }

  # if drug investigated for a complex disease
  if (sum(drug %in% ct_filt_blair$db_id) > 0) {
    drug_conditions_investigated = ct_filt_blair %>% filter(db_id == drug) %>% dplyr::select(phase, complex_disease) %>% distinct()
    
    for (b in 1:nrow(drug_conditions_investigated)) {

      disease = drug_conditions_investigated[b, complex_disease]

      # if the drug is already indicated, then next
      if (setequal(ultimate_per_disease[i, disease], 5) == TRUE) {
        next
      }

      # if no previous information about this drug-complex_disease pair, add clinical trial phase
      if (is.na(ultimate_per_disease[i, disease]) == TRUE) {
        ultimate_per_disease[i, disease] = drug_conditions_investigated[b, phase]
        next
      }

      # if previous information about this drug-complex_disease pair, keep the most advanced clinical trial phase
      if (is.na(ultimate_per_disease[i, disease]) == FALSE) {
        if (ultimate_per_disease[i, disease] >= drug_conditions_investigated[b, phase]) {
          next
        }
        if (ultimate_per_disease[i, disease] < drug_conditions_investigated[b, phase]) {
          ultimate_per_disease[i, disease] = drug_conditions_investigated[b, phase]
        }
      }
    }
  }
  cat(i, "/", nrow(ultimate_per_disease), "\n")
}
rm(drug_conditions_indicated, drug_conditions_investigated, b, disease, drug, i, z)

ultimate_per_disease = as.data.frame(ultimate_per_disease)
ultimate_per_disease = ultimate_per_disease %>% mutate(drugbank_id = rownames(ultimate_per_disease), .before = Acne)
rownames(ultimate_per_disease) = NULL

# save file
data.table::fwrite(ultimate_per_disease, "processed_data/drugs_inv_ind_per_disease.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

####################
##                ##
## visualizations ##
##                ##
####################

## reshape
ultimate_per_disease = reshape2::melt(ultimate_per_disease, "drugbank_id", colnames(ultimate_per_disease)[2:ncol(ultimate_per_disease)])
ultimate_per_disease = na.omit(ultimate_per_disease)
colnames(ultimate_per_disease) = c("drugbank_id", "complex_disease", "phase")

## correct clinical trial phases
ultimate_per_disease$phase = gsub(1.5, 2, ultimate_per_disease$phase)
ultimate_per_disease$phase = gsub(2.5, 3, ultimate_per_disease$phase)
ultimate_per_disease$phase = gsub(5, 4, ultimate_per_disease$phase) # 4 is considered to be indicated drugs
table(ultimate_per_disease$phase)
# "Hypotony of the eye" has one clinical trial but it has unknown phase --> for illustration purposes, convert this to 1
ultimate_per_disease[which(ultimate_per_disease$complex_disease == "Hypotony of the Eye"), "phase"] = 1

## remove clinical trials with unknown phases
ultimate_per_disease = ultimate_per_disease[-which(ultimate_per_disease$phase == 0), ] ; rownames(ultimate_per_disease) = NULL

## calculate number of drugs per clinical trial phase
ultimate_per_disease = ultimate_per_disease %>% 
  group_by(complex_disease) %>%
  mutate(nr_drugs = length(drugbank_id),
         phase_1 = sum(phase == 1),
         phase_2 = sum(phase == 2),
         phase_3 = sum(phase == 3),
         approved = sum(phase == 4)) %>%
  ungroup() %>%
  dplyr::select(complex_disease, phase_1, phase_2, phase_3, approved) %>%
  distinct()

## load complex disease categories
cd_categories = fread("raw_data/complex_disease_category.txt")

## gather data
ultimate_per_disease_gathered = tidyr::gather(ultimate_per_disease, status, nr_drugs, -complex_disease)
ultimate_per_disease_gathered$status = factor(ultimate_per_disease_gathered$status, levels = c("approved", "phase_3", "phase_2", "phase_1"), labels = c("Approved", "Phase 3", "Phase 2", "Phase 1"))

## add complex disease category
ultimate_per_disease_gathered = left_join(ultimate_per_disease_gathered, cd_categories, by = "complex_disease")

## calculate total number of investigated/indicated drugs per complex disease
ultimate_per_disease_gathered = ultimate_per_disease_gathered %>% 
  group_by(complex_disease) %>% 
  mutate(total_drugs = sum(nr_drugs)) %>% 
  ungroup() %>%
  distinct() %>%
  arrange(disease_category, total_drugs, status)
ultimate_per_disease_gathered$complex_disease = factor(ultimate_per_disease_gathered$complex_disease, levels = ultimate_per_disease_gathered$complex_disease, labels = ultimate_per_disease_gathered$complex_disease)
ultimate_per_disease_gathered$disease_category = factor(ultimate_per_disease_gathered$disease_category, levels = unique(ultimate_per_disease_gathered$disease_category), labels = unique(ultimate_per_disease_gathered$disease_category))

myColors <- c("brown", "blue3", "cyan3", "orange", "purple", "chartreuse3")
names(myColors) <- levels(ultimate_per_disease_gathered$disease_category)

scale_custom <- list(
  scale_fill_manual(values=c("#808080", "#A9A9A9", "#D3D3D3", "#E5E4E2")),
  scale_color_manual(aesthetics = "complex_disease", values = myColors)
)

fig1c = ggplot(ultimate_per_disease_gathered, aes(y = complex_disease, x = nr_drugs, fill = status)) +
  geom_col() +
  xlab("Number of drugs\n") +
  ylab("") +
  labs(fill = "") +
  scale_x_continuous(breaks = seq(0, 600, 50)) +
  scale_fill_manual(values=c("#808080", "#A9A9A9", "#D3D3D3", "#E5E4E2")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 25, family = "Arial", color = "black"),
        axis.text.y = element_text(size = 30, family = "Arial", color = c(rep("brown", 14),
                                                                          rep("blue3", 4),
                                                                          rep("cyan3", 8),
                                                                          rep("orange", 19),
                                                                          rep("purple", 15),
                                                                          rep("chartreuse3", 5))),
        axis.title.x = element_text(size = 30, family = "Arial", color = "black"),
        legend.text = element_text(size = 40)) +
  guides(fill = guide_legend(byrow = TRUE))

fig1c
ggsave(filename = "Fig1C_investigated_indicated_drugs_per_complex_disease.tiff", 
       path = "figures/", 
       width = 20, height = 24, device = "tiff",
       dpi = 700, compression = "lzw", type = type_compression)
dev.off()
