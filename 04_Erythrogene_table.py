import pandas as pd

## Main Script ##

# Read files
file = "erythrogene_alleles.xlsx"
HGVS_df = pd.read_csv('HGVS_Notation.txt', sep='\t')
temp_df = pd.read_excel(file)

# Provide column order
df = pd.DataFrame(columns=["Blood System", "Phenotype", "Allele Name", "Gene", "Chromosome", 
                           "Nucleotide Change", "GRCh38 Coordinates", "GRCh37 Coordinates",
                           "GRCh38 VCF Position", "GRCh37 VCF Position"])

# Assign values to df columns using temp_df and HGVS_df
df["Gene"] = temp_df["Gene"]
df["Nucleotide Change"] = temp_df["Nucleotide Change"]
df["Phenotype"] = "-"
df["Allele Name"] = temp_df["Allele"] + " (" + temp_df["Gene"] + "." + (temp_df.groupby("Gene").cumcount() + 1).astype(str).str.zfill(2) + ")"
df["Blood System"] = df["Gene"].map(HGVS_df.set_index('Gene')['Blood Group'])
df["Chromosome"] = df["Gene"].map(HGVS_df.set_index('Gene')['Chromosome'])

# Each phenotype may be associated with > 1 SNP. Split the string such that one row has only 1 nucleotide change.
df['Nucleotide Change'] = df['Nucleotide Change'].str.split(';')
df = df.explode('Nucleotide Change')
df['Nucleotide Change'] = df['Nucleotide Change'].str.strip()

# Add HGVS notation and "c." in front in "Nucleotide Change" column
df["Nucleotide Change"] = "c." + df["Nucleotide Change"]
df['Transcript'] = df['Gene'].map(HGVS_df.set_index('Gene')['Transcript'])
df['Nucleotide Change'] = df.apply(
    lambda row: row['Transcript'] + ":" + row['Nucleotide Change'] if pd.notna(row['Transcript']) and row['Nucleotide Change'].startswith('c') else row['Nucleotide Change'],
    axis=1)
df.drop(columns='Transcript', inplace=True)

# Save dataframe in excel and tsv format
df.to_excel("Extracted Tables/Erythrogene_Table.xlsx", index=False)
df.to_csv("Extracted Tables/Erythrogene_Table.tsv", sep='\t', index=False)

print("Done")