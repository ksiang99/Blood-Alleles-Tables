# Add nucleotide changes not stated in ABO B alleles and RHCE alleles

import pandas as pd

file_ABO = "Extracted Tables/001_ABO.tsv"
file_RHCE = "Extracted Tables/004_RH (RHCE).tsv"
df_ABO = pd.read_csv(file_ABO, sep='\t')
df_RHCE = pd.read_csv(file_RHCE, sep='\t')

# For ABO B alleles
copy_NC_lst = df_ABO.loc[df_ABO["Allele Name"] == "ABO*B.01", "Nucleotide Change"].tolist()
copy_NC_str = ";".join(copy_NC_lst)
start_ind = df_ABO.loc[df_ABO["Allele Name"] == "ABO*B.02"].index[0]
end_ind = df_ABO.loc[df_ABO["Allele Name"] == "ABO*BEL.05"].index[0]

for ind in range(start_ind, end_ind + 1):
    df_ABO.loc[ind, "Nucleotide Change"] += ";" + copy_NC_str

df_ABO['Nucleotide Change'] = df_ABO['Nucleotide Change'].str.split(';')
df_ABO = df_ABO.explode('Nucleotide Change')
df_ABO = df_ABO.drop_duplicates()
df_ABO['Phenotype'] = df_ABO['Phenotype'].ffill()
df_ABO['Allele Name'] = df_ABO['Allele Name'].ffill()

df_ABO.to_excel("Corrected ISBT Tables/001_ABO.xlsx", index=False)
df_ABO.to_csv("Corrected ISBT Tables/001_ABO.tsv", sep='\t', index=False)

# For RHCE*02 alleles
copy_NC_lst = df_RHCE.loc[df_RHCE["Allele Name"] == "RHCE*02", "Nucleotide Change"].tolist()
copy_NC_str = ";".join(copy_NC_lst)
start_ind = df_RHCE.loc[df_RHCE["Allele Name"] == "RHCE*02.01 RHCE*Ce.01 RHCE*CeMA RHCE*CeJAL"].index[0]
end_ind = df_RHCE.loc[df_RHCE["Allele Name"] == "RHCE*02N.14 RHCE*CeN.14"].index[0]

for ind in range(start_ind, end_ind + 1):
    df_RHCE.loc[ind, "Nucleotide Change"] += ";" + copy_NC_str

# For RHCE*03 alleles
copy_NC_lst = df_RHCE.loc[df_RHCE["Allele Name"] == "RHCE*03", "Nucleotide Change"].tolist()
copy_NC_str = ";".join(copy_NC_lst)
start_ind = df_RHCE.loc[df_RHCE["Allele Name"] == "RHCE*03.01 RHCE*cE.01 RHCE*cEEW"].index[0]
end_ind = df_RHCE.loc[df_RHCE["Allele Name"] == "RHCE*03N.07 RHCE*cEN.07"].index[0]

for ind in range(start_ind, end_ind + 1):
    df_RHCE.loc[ind, "Nucleotide Change"] += ";" + copy_NC_str

# For RHCE*04 alleles
copy_NC_lst = df_RHCE.loc[df_RHCE["Allele Name"] == "RHCE*CE", "Nucleotide Change"].tolist()
copy_NC_str = ";".join(copy_NC_lst)
start_ind = df_RHCE.loc[df_RHCE["Allele Name"] == "RHCE*04.01 RHCE*CE.01"].index[0]
end_ind = df_RHCE.loc[df_RHCE["Allele Name"] == "RHCE*04.02 RHCE*CE.02"].index[0]

for ind in range(start_ind, end_ind + 1):
    df_RHCE.loc[ind, "Nucleotide Change"] += ";" + copy_NC_str

df_RHCE['Nucleotide Change'] = df_RHCE['Nucleotide Change'].str.split(';')
df_RHCE = df_RHCE.explode('Nucleotide Change')
df_RHCE = df_RHCE.drop_duplicates()
df_RHCE['Phenotype'] = df_RHCE['Phenotype'].ffill()
df_RHCE['Allele Name'] = df_RHCE['Allele Name'].ffill()

df_RHCE = df_RHCE[~df_RHCE["Nucleotide Change"].isin(["NM_020485.8:c.307C", "NM_020485.8:c.676G"]) | df_RHCE["Allele Name"].str.contains("RHCE\\*01", regex=True)]

df_RHCE.to_excel("Corrected ISBT Tables/004_RH (RHCE).xlsx", index=False)
df_RHCE.to_csv("Corrected ISBT Tables/004_RH (RHCE).tsv", sep='\t', index=False)

print("Done")