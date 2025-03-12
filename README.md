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
* Adjust VCF positions for `ABO:c.261G` Provides the correct VCF position
* Corrects genotype values for biallelic/multiallelic sites.

**`07_blood_allele_table.py`**
* Produces two final output files:
  * `Blood_Allele_Table.tsv`
  * `Blood_Allele_Table_Separated.tsv` (Exploded version of `Blood_Allele_Table.tsv`)
* Uses the following files:
  * `ISBT_variants_to_remove.tsv` and `Erythrogene_variants_to_remove.tsv` to remove alleles with nucleotide change that cannot be found or automapped by VariantValidator, have variant reference not aligning with reference sequence or involves Exon/Intron deletion.
  * `ISBT_variants_to_overwrite.tsv` provides the mapping details of nucleotide change that cannot be found using **`05_get_coords.py`** but works when using [VariantValidator Web Interface](https://rest.variantvalidator.org/)

**`08_Infer.R`**
* Uses genotype data from 1000G VCF files and `Blood_Allele_Table_Separated.tsv` to infer potential phenotypes for each position.
