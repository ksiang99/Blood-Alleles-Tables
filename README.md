# **Python Packages**
`Python 3.11.10` with `Tabula 2.9.3`, `Pandas 2.2.3`, `Request 2.32.3` and `Numpy 1.26.4` are used to run the Python scripts.

# **Files**
**`01_download_ISBT_Table.py`** 
* Downloads ISBT Blood Alleles Tables using URLS provided in `ISBT_links.txt`.

**`02_extract_ISBT_table.py`**
* Extracts and cleans the **Phenotype**, **Allele Name** and **Nucleotide Change** columns in ISBT Blood Alleles Tables.
* Uses HGVS_Notation.txt to create **Blood Group**, **Gene** and **Chromosome** columns.
* Appends HGVS transcript notation in front of **Nucleotide Change** column.

**`03_ABO_add.py`**
* Adds missing nucleotide change to B alleles in the ABO blood group.

**`04_Erythrogene_table.py`**
* Generates a table with the same structure as the output of `02_extract_ISBT_table.py` using `erythrogene_alleles.xlsx`.

**`05_get_coords.py`** 
* Uses [VariantValidator](https://rest.variantvalidator.org/) to:
  * Map the **Nucleotide Change** column in extracted ISBT and Erythrogene tables to GRCh37 and GRCh38 coordinates and VCF position.
  * Create two new columns **GRCh37 Alt Allele** and **GRCh38 Alt Allele** to represent the genotype value of the nucleotide change.

**`06_edit_genotype_value.py`**
* Adjust VCF positions for `ABO:c.261G`.
* Corrects genotype values for biallelic/multiallelic sites.

**`07_blood_allele_table.py`**
* Produces two final output files:
  * `Blood_Allele_Table.tsv`
  * `Blood_Allele_Table_Separated.tsv` (Exploded version)
* Uses the following files:
  * `ISBT_variants_to_remove.tsv` and `Erythrogene_variants_to_remove.tsv` to remove alleles with nucleotide change that cannot be found or automapped by VariantValidator, have variant reference not aligning with reference sequence or involves Exon/Intron deletion.
  * `ISBT_variants_to_overwrite.tsv` provides the mapping details of nucleotide change that cannot be found using **`05_get_coords.py`** but works when using [VariantValidator Web Interface](https://rest.variantvalidator.org/)

**`08_Infer.R`**
* Uses genotype data stored in `chr[]_df_genotype.Rdata` and allele information from `Blood_Allele_Table_Separated.tsv` to infer potential phenotypes at each variant position.

**`09_summary_table.R`**
* Condenses the results from `08_Infer.R` into a summary table with blood group as columns and samples as rows.
* List the potential phenotypes for each blood group and sample.

**`10_final_summary_table.R`**
* From `09_summary_table.R`, the phenotype with the most number of Nucleotide Change associate is inferred as the "correct" phenotype for the sample.

**`11_phenotype_count.R`**
* From `10_final_summary_table.R`, calculate the positive phenotype case for ALL population and superpopulations.

**`12_infer_pred.R`**
* Works similiar to `08_Infer.R` but a position in df_genotype is replaced with predicted genotype values from a Machine Learning Model.

**`13_summary_table_pred.R`**
* Condenses each result from `12_infer_pred.R` into a summary table with blood group as columns and samples as rows.
* List the potential phenotypes for the blood group which the predicted genotype values of the position is associated to.

**`14_final_summary_table_pred.R`**
* Works similiar to `10_final_summary_table.R`

**`15_phenotype_count.R`**
* Works similiar to `11_phenotype_count.R`
