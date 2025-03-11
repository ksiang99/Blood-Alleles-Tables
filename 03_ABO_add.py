# Add nucleotide changes not stated in ABO B alleles

import pandas as pd

file = "Extracted ISBT Tables/001_ABO.tsv"

df = pd.read_csv(file, sep='\t')

copy_NC_lst = df.loc[df["Allele Name"] == "ABO*B.01", "Nucleotide Change"].tolist()
copy_NC_str = ";".join(copy_NC_lst)
start_ind = df.loc[df["Allele Name"] == "ABO*B.02"].index[0]
end_ind = df.loc[df["Allele Name"] == "ABO*BEL.05"].index[0]

for ind in range(start_ind, end_ind + 1):
    df.loc[ind, "Nucleotide Change"] += ";" + copy_NC_str

df['Nucleotide Change'] = df['Nucleotide Change'].str.split(';')
df = df.explode('Nucleotide Change')
df = df.drop_duplicates()
df['Phenotype'] = df['Phenotype'].ffill()
df['Allele Name'] = df['Allele Name'].ffill()

df.to_excel("Corrected ISBT Tables/001_ABO.xlsx", index=False)
df.to_csv("Corrected ISBT Tables/001_ABO.tsv", sep='\t', index=False)

print("Done")