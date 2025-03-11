import pandas as pd

file_erythro = "Tables with coords/Erythrogene_Table_coords.tsv"
file_ABO = "Tables with coords/001_ABO_coords.tsv"
file_JMH = "Tables with coords/026_JMH_coords.tsv"
file_IN = "Tables with coords/023_IN_coords.tsv"
file_CROM = "Tables with coords/021_CROM_coords.tsv"
file_LU = "Tables with coords/005_LU_coords.tsv"
file_KLF1 = "Tables with coords/046_KLF1_coords.tsv"
file_H = "Tables with coords/018_H_coords.tsv"


df_e = pd.read_csv(file_erythro, sep='\t')
df_ABO = pd.read_csv(file_ABO, sep='\t')
df_JMH = pd.read_csv(file_JMH, sep='\t')
df_IN = pd.read_csv(file_IN, sep='\t')
df_CROM = pd.read_csv(file_CROM, sep='\t')
df_LU = pd.read_csv(file_LU, sep='\t')
df_KLF1 = pd.read_csv(file_KLF1, sep='\t')
df_H = pd.read_csv(file_H, sep='\t')

# Change genotype value of c.261delG in ABO (ISBT)
indices = df_ABO[df_ABO['Nucleotide Change'] == 'NM_020469.2:c.261delG'].index
df_ABO.loc[indices, 'GRCh38 Coordinates'] = 'NC_000009.12:g.133257521delC'
df_ABO.loc[indices, 'GRCh37 Coordinates'] = 'NC_000009.11:g.136132908delC'
df_ABO.loc[indices, 'GRCh38 VCF Position'] = '133257521'
df_ABO.loc[indices, 'GRCh37 VCF Position'] = '136132908'
df_ABO.loc[indices, 'GRCh38 Alt Allele'] = '0'
df_ABO.loc[indices, 'GRCh37 Alt Allele'] = '0'
df_ABO = df_ABO.astype(str)
df_ABO.loc[0, 'GRCh38 Coordinates'] = ""
df_ABO.loc[0, 'GRCh37 Coordinates'] = ""
df_ABO.loc[0, 'GRCh38 VCF Position'] = ""
df_ABO.loc[0, 'GRCh37 VCF Position'] = ""
df_ABO.to_excel("Tables with coords/001_ABO_coords.xlsx", index=False)
df_ABO.to_csv("Tables with coords/001_ABO_coords.tsv", sep='\t', index=False)

# Change genotype value of c.1545A>G in JMH
indices = df_JMH[df_JMH['Nucleotide Change'] == 'NM_003612.3:c.1545A>G'].index
df_JMH.loc[indices, 'GRCh38 Alt Allele'] = '2'
df_JMH.loc[indices, 'GRCh37 Alt Allele'] = '2'
df_JMH = df_JMH.astype(str)
df_JMH.loc[0, 'GRCh38 Coordinates'] = ""
df_JMH.loc[0, 'GRCh37 Coordinates'] = ""
df_JMH.loc[0, 'GRCh38 VCF Position'] = ""
df_JMH.loc[0, 'GRCh37 VCF Position'] = ""
df_JMH.to_excel("Tables with coords/026_JMH_coords.xlsx", index=False)
df_JMH.to_csv("Tables with coords/026_JMH_coords.tsv", sep='\t', index=False)

# Change genotype value of c.255C>G in IN
indices = df_IN[df_IN['Nucleotide Change'] == 'NM_001001391.1:c.255C>G'].index
df_IN.loc[indices, 'GRCh38 Alt Allele'] = '2'
df_IN.loc[indices, 'GRCh37 Alt Allele'] = '2'
df_IN = df_IN.astype(str)
df_IN.loc[1, 'GRCh38 Coordinates'] = ""
df_IN.loc[1, 'GRCh37 Coordinates'] = ""
df_IN.loc[1, 'GRCh38 VCF Position'] = ""
df_IN.loc[1, 'GRCh37 VCF Position'] = ""
df_IN.to_excel("Tables with coords/023_IN_coords.xlsx", index=False)
df_IN.to_csv("Tables with coords/023_IN_coords.tsv", sep='\t', index=False)

# Change genotype value of c.155G>T in CROM
indices = df_CROM[df_CROM['Nucleotide Change'] == 'NM_000574.5:c.155G>T'].index
df_CROM.loc[indices, 'GRCh38 Alt Allele'] = '2'
df_CROM.loc[indices, 'GRCh37 Alt Allele'] = '2'
df_CROM = df_CROM.astype(str)
df_CROM.loc[0, 'GRCh38 Coordinates'] = ""
df_CROM.loc[0, 'GRCh37 Coordinates'] = ""
df_CROM.loc[0, 'GRCh38 VCF Position'] = ""
df_CROM.loc[0, 'GRCh37 VCF Position'] = ""
df_CROM.to_excel("Tables with coords/021_CROM_coords.xlsx", index=False)
df_CROM.to_csv("Tables with coords/021_CROM_coords.tsv", sep='\t', index=False)

# Change genotype value of c.711C>A in LU
indices = df_LU[df_LU['Nucleotide Change'] == 'NM_005581.4:c.711C>A'].index
df_LU.loc[indices, 'GRCh38 Alt Allele'] = '2'
df_LU.loc[indices, 'GRCh37 Alt Allele'] = '2'
df_LU = df_LU.astype(str)
df_LU.loc[0, 'GRCh38 Coordinates'] = ""
df_LU.loc[0, 'GRCh37 Coordinates'] = ""
df_LU.loc[0, 'GRCh38 VCF Position'] = ""
df_LU.loc[0, 'GRCh37 VCF Position'] = ""
df_LU.to_excel("Tables with coords/005_LU_coords.xlsx", index=False)
df_LU.to_csv("Tables with coords/005_LU_coords.tsv", sep='\t', index=False)

# Change genotype value of in KLF1
indices = df_KLF1[df_KLF1['Nucleotide Change'] == 'NM_006563.3:c.809C>A'].index
df_KLF1.loc[indices, 'GRCh38 Alt Allele'] = '2'
df_KLF1.loc[indices, 'GRCh37 Alt Allele'] = '2'
indices = df_KLF1[df_KLF1['Nucleotide Change'] == 'NM_006563.3:c.114delC'].index
df_KLF1.loc[indices, 'GRCh38 Alt Allele'] = '2'
df_KLF1.loc[indices, 'GRCh37 Alt Allele'] = '2'
indices = df_KLF1[df_KLF1['Nucleotide Change'] == 'NM_006563.3:c.517_519delCCC'].index
df_KLF1.loc[indices, 'GRCh38 Alt Allele'] = '2'
df_KLF1.loc[indices, 'GRCh37 Alt Allele'] = '2'
indices = df_KLF1[df_KLF1['Nucleotide Change'] == 'NM_006563.3:c.519_520insC'].index
df_KLF1.loc[indices, 'GRCh38 Alt Allele'] = '3'
df_KLF1.loc[indices, 'GRCh37 Alt Allele'] = '3'
indices = df_KLF1[df_KLF1['Nucleotide Change'] == 'NM_006563.3:c.310_311insG'].index
df_KLF1.loc[indices, 'GRCh38 Alt Allele'] = '2'
df_KLF1.loc[indices, 'GRCh37 Alt Allele'] = '2'
df_KLF1 = df_KLF1.astype(str)
df_KLF1.loc[0, 'GRCh38 Coordinates'] = ""
df_KLF1.loc[0, 'GRCh37 Coordinates'] = ""
df_KLF1.loc[0, 'GRCh38 VCF Position'] = ""
df_KLF1.loc[0, 'GRCh37 VCF Position'] = ""
df_KLF1.to_excel("Tables with coords/046_KLF1_coords.xlsx", index=False)
df_KLF1.to_csv("Tables with coords/046_KLF1_coords.tsv", sep='\t', index=False)

# Change genotype value of c.13_19dup in LU
indices = df_H[df_H['Nucleotide Change'] == 'NM_000148.3:c.13_19dup'].index
df_H.loc[indices, 'GRCh38 Alt Allele'] = '2'
df_H.loc[indices, 'GRCh37 Alt Allele'] = '2'
df_H = df_H.astype(str)
df_H.loc[0, 'GRCh38 Coordinates'] = ""
df_H.loc[0, 'GRCh37 Coordinates'] = ""
df_H.loc[0, 'GRCh38 VCF Position'] = ""
df_H.loc[0, 'GRCh37 VCF Position'] = ""
df_H.to_excel("Tables with coords/018_H_coords.xlsx", index=False)
df_H.to_csv("Tables with coords/018_H_coords.tsv", sep='\t', index=False)

# Change genotype value in Erythro
indices = df_e[df_e['Nucleotide Change'] == 'NM_020469.2:c.261delG'].index
df_e.loc[indices, 'GRCh38 Coordinates'] = 'NC_000009.12:g.133257521delC'
df_e.loc[indices, 'GRCh37 Coordinates'] = 'NC_000009.11:g.136132908delC'
df_e.loc[indices, 'GRCh38 VCF Position'] = '133257521'
df_e.loc[indices, 'GRCh37 VCF Position'] = '136132908'
df_e.loc[indices, 'GRCh38 Alt Allele'] = '0'
df_e.loc[indices, 'GRCh37 Alt Allele'] = '0'

indices = df_e[df_e['Nucleotide Change'] == 'NM_003612.3:c.1545A>G'].index
df_e.loc[indices, 'GRCh38 Alt Allele'] = '2'
df_e.loc[indices, 'GRCh37 Alt Allele'] = '2'

indices = df_e[df_e['Nucleotide Change'] == 'NM_002049.3:c.1067G>T'].index
df_e.loc[indices, 'GRCh38 Alt Allele'] = '2'
df_e.loc[indices, 'GRCh37 Alt Allele'] = '2'

df_e = df_e.astype(str)

df_e.to_excel("Tables with coords/Erythrogene_Table_coords.xlsx", index=False)
df_e.to_csv("Tables with coords/Erythrogene_Table_coords.tsv", sep='\t', index=False)

print("Done")