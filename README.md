# **Files**
**`01_download_ISBT_Table.py`:** Downloads ISBT Blood Alleles Tables using `ISBT_links.txt`

**`02_extract_ISBT_table.py`:** Extract and clean `Phenotype`, `Allele Name` and `Nucleotide Change` columns in ISBT Blood Alleles Tables. Uses HGVS_Notation.txt to add `Blood Group`, `Gene` and `Chromosome` columns, and add HGVS transcript notation in front of `Nucleotide Change`.

**`03_ABO_add.py`**: Add missing nucleotide change in B alleles.

**`04_Erythrogene_table.py`**: Generate the same table structure as output of `02_extract_ISBT_table.py` for erythrogene alleles.
