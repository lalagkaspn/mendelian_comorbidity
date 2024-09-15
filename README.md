# Shared etiology of Mendelian and complex disease supports drug discovery

## Lalagkas, P.N., Melamed, R.D. Shared etiology of Mendelian and complex disease supports drug discovery. BMC Med Genomics 17, 228 (2024). https://doi.org/10.1186/s12920-024-01988-3

Link to the [article](https://bmcmedgenomics.biomedcentral.com/articles/10.1186/s12920-024-01988-3).

In this work, we analyze a set of *2,908* clinically associated (comorbid) pairs of *90* Mendelian diseases and *65* complex diseases. Our goal is to leverage the well-known biology of Mendelian diseases to improve the treatment of complex diseases. To this end, we use the aforementioned clinical associations to match each complex disease to a set of Mendelian disease causal genes. We hypothesize that the drugs targeting these genes are potential candidate drugs for the complex diseases. We validate our list of candidate drugs using thousands of clinical trials and current drug indications. Additionally, we present evidence that certain Mendelian diseases hold a greater promise as resources for complex disease drug discovery. We also assess the genetic similarity between 60 Mendelian diseases and 10 cancers, as another evidence of shared biology (other than the comorbidity), and propose that the combination of comorbidity and genetic similarity can prioritize candidate drugs with a higher probability of success.

There are eight scripts in this repository. The first three scripts preprocess and create the necessary data to run the last five scripts.

- `md_cd_comorbidity_pairs_preprocessing.R` is the script to preprocess the clinically associated pairs of Mendelian and complex diseases
- `drugs_investigated_indicated_diseases.R` is the script to get all the drugs that are investigated or indicated for the 65 complex diseases
- `drug_targets_drugbank.R` is the script to preprocess DrugBank files to get all the drugs with known drug targets
- `main_analysis.R` is the script to reproduce the results and figures of the second section of the Results (*Mendelian disease comorbidity identifies drugs under current investigation or indication*)
- `per_mendelian_disease_analysis.R` is the script to reproduce the results and figures of the third section of the Results (*Prioritizing Mendelian diseases targeted by high number of drugs*)
- `genetic_similarity.R` is the script to reproduce the results and figures of the fourth section of the Results (*Combining comorbidity with genetic similarity enhances drug predictions*)
- `supplementary_tables.R` is the script to reproduce the supplementary tables
- `supplementary_figures.R` is the script to reproduce the supplementary figures

All scripts are written in **R**.

NOTE: all the publicly available data needed to run the above scripts is made available in this repository. For the non-publicly available data, information about how to download it after getting the corresponding licenses is included in each script.

If you have any questions, please reach out to [panagiotis.lalagkas@gmail.com](mailto:panagiotis.lalagkas@gmail.com)


