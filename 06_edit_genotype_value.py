import pandas as pd

file_erythro = "Tables with coords/Erythrogene_Table_coords.tsv"
file_ABO = "Tables with coords/001_ABO_coords.tsv"
file_JMH = "Tables with coords/026_JMH_coords.tsv"
file_IN = "Tables with coords/023_IN_coords.tsv"
file_CROM = "Tables with coords/021_CROM_coords.tsv"
file_LU = "Tables with coords/005_LU_coords.tsv"
file_KLF1 = "Tables with coords/046_KLF1_coords.tsv"
file_H = "Tables with coords/018_H_coords.tsv"
file_CO = "Tables with coords/015_CO_coords.tsv"
file_KEL = "Tables with coords/006_KEL_coords.tsv"
file_JR = "Tables with coords/032_JR_coords.tsv"
file_LAN = "Tables with coords/033_LAN_coords.tsv"
file_RHD = "Tables with coords/004_RH (RHD)_coords.tsv"
file_RHCE = "Tables with coords/004_RH (RHCE)_coords.tsv"

df_e = pd.read_csv(file_erythro, sep='\t')
df_ABO = pd.read_csv(file_ABO, sep='\t')
df_JMH = pd.read_csv(file_JMH, sep='\t')
df_IN = pd.read_csv(file_IN, sep='\t')
df_CROM = pd.read_csv(file_CROM, sep='\t')
df_LU = pd.read_csv(file_LU, sep='\t')
df_KLF1 = pd.read_csv(file_KLF1, sep='\t')
df_H = pd.read_csv(file_H, sep='\t')
df_CO = pd.read_csv(file_CO, sep='\t')
df_KEL = pd.read_csv(file_KEL, sep='\t')
df_JR = pd.read_csv(file_JR, sep='\t')
df_LAN = pd.read_csv(file_LAN, sep='\t')
df_RHD = pd.read_csv(file_RHD, sep='\t')
df_RHCE = pd.read_csv(file_RHCE, sep='\t')

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
df_LU.loc[5, 'GRCh38 Coordinates'] = ""
df_LU.loc[5, 'GRCh37 Coordinates'] = ""
df_LU.loc[5, 'GRCh38 VCF Position'] = ""
df_LU.loc[5, 'GRCh37 VCF Position'] = ""
df_LU.to_excel("Tables with coords/005_LU_coords.xlsx", index=False)
df_LU.to_csv("Tables with coords/005_LU_coords.tsv", sep='\t', index=False)

# Change genotype value of in KLF1
indices = df_KLF1[df_KLF1['Nucleotide Change'].isin(['NM_006563.3:c.809C>A', 'NM_006563.3:c.114delC',
                                                     'NM_006563.3:c.517_519delCCC', 'NM_006563.3:c.310_311insG',
                                                     'NM_006563.3:c.519_525dupCGGCGC C'])].index
df_KLF1.loc[indices, 'GRCh38 Alt Allele'] = '2'
df_KLF1.loc[indices, 'GRCh37 Alt Allele'] = '2'
indices = df_KLF1[df_KLF1['Nucleotide Change'] == 'NM_006563.3:c.519_520insC'].index
df_KLF1.loc[indices, 'GRCh38 Alt Allele'] = '3'
df_KLF1.loc[indices, 'GRCh37 Alt Allele'] = '3'
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

# Change genotype value of c.601delG in CO
indices = df_CO[df_CO['Nucleotide Change'] == 'NM_198098.2:c.601delG'].index
df_CO.loc[indices, 'GRCh38 Alt Allele'] = '2'
df_CO.loc[indices, 'GRCh37 Alt Allele'] = '2'
df_CO = df_CO.astype(str)
df_CO.loc[0, 'GRCh38 Coordinates'] = ""
df_CO.loc[0, 'GRCh37 Coordinates'] = ""
df_CO.loc[0, 'GRCh38 VCF Position'] = ""
df_CO.loc[0, 'GRCh37 VCF Position'] = ""
df_CO.to_excel("Tables with coords/015_CO_coords.xlsx", index=False)
df_CO.to_csv("Tables with coords/015_CO_coords.tsv", sep='\t', index=False)

# Change genotype value in KEL
indices = df_KEL[df_KEL['Nucleotide Change'].isin(['NM_000420.3:c.904delG', 'NM_000420.3:c.924+1G>A',
                                                   'NM_000420.3:c.1477C>T', 'NM_000420.3:c.1546C>T',
                                                   'NM_000420.3:c.2175delC', 'NM_000420.3:c.2107G>A',
                                                   'NM_000420.3:c.578C>G'])].index
df_KEL.loc[indices, 'GRCh38 Alt Allele'] = '2'
df_KEL.loc[indices, 'GRCh37 Alt Allele'] = '2'
df_KEL = df_KEL.astype(str)
df_KEL.loc[4, 'GRCh38 Coordinates'] = ""
df_KEL.loc[4, 'GRCh37 Coordinates'] = ""
df_KEL.loc[4, 'GRCh38 VCF Position'] = ""
df_KEL.loc[4, 'GRCh37 VCF Position'] = ""
df_KEL.to_excel("Tables with coords/006_KEL_coords.xlsx", index=False)
df_KEL.to_csv("Tables with coords/006_KEL_coords.tsv", sep='\t', index=False)

# Change genotype value of c.420dupA in JR
indices = df_JR[df_JR['Nucleotide Change'] == 'NM_004827.3:c.420dupA'].index
df_JR.loc[indices, 'GRCh38 Alt Allele'] = '2'
df_JR.loc[indices, 'GRCh37 Alt Allele'] = '2'
df_JR = df_JR.astype(str)
df_JR.loc[0, 'GRCh38 Coordinates'] = ""
df_JR.loc[0, 'GRCh37 Coordinates'] = ""
df_JR.loc[0, 'GRCh38 VCF Position'] = ""
df_JR.loc[0, 'GRCh37 VCF Position'] = ""
df_JR.to_excel("Tables with coords/032_JR_coords.xlsx", index=False)
df_JR.to_csv("Tables with coords/032_JR_coords.tsv", sep='\t', index=False)

# Change genotype value of c.301dupG in JR
indices = df_LAN[df_LAN['Nucleotide Change'] == 'NM_005689.4:c.301dupG'].index
df_LAN.loc[indices, 'GRCh38 Alt Allele'] = '2'
df_LAN.loc[indices, 'GRCh37 Alt Allele'] = '2'
df_LAN = df_LAN.astype(str)
df_LAN.loc[0, 'GRCh38 Coordinates'] = ""
df_LAN.loc[0, 'GRCh37 Coordinates'] = ""
df_LAN.loc[0, 'GRCh38 VCF Position'] = ""
df_LAN.loc[0, 'GRCh37 VCF Position'] = ""
df_LAN.to_excel("Tables with coords/033_LAN_coords.xlsx", index=False)
df_LAN.to_csv("Tables with coords/033_LAN_coords.tsv", sep='\t', index=False)

# Change genotype value in RHD
indices = df_RHD[df_RHD['Nucleotide Change'].isin(["NM_016124.6:c.697G>A", "NM_016124.6:c.410C>A",
                                                   "NM_016124.6:c.809T>A", "NM_016124.6:c.745_757del13",
                                                   "NM_016124.6:c.510_511insG", "NM_016124.6:c.203G>C"])].index
df_RHD.loc[indices, 'GRCh38 Alt Allele'] = '2'
df_RHD.loc[indices, 'GRCh37 Alt Allele'] = '2'
df_RHD = df_RHD.astype(str)
df_RHD.loc[0, 'GRCh38 Coordinates'] = ""
df_RHD.loc[0, 'GRCh37 Coordinates'] = ""
df_RHD.loc[0, 'GRCh38 VCF Position'] = ""
df_RHD.loc[0, 'GRCh37 VCF Position'] = ""
df_RHD.to_excel("Tables with coords/004_RH (RHD)_coords.xlsx", index=False)
df_RHD.to_csv("Tables with coords/004_RH (RHD)_coords.tsv", sep='\t', index=False)

# Change genotype value in RHCE
indices = df_RHCE[df_RHCE['Nucleotide Change'].isin(["NM_020485.8:c.818C>A", "NM_020485.8:c.500T>A",
                                                     "NM_020485.8:c.380C>A"])].index
df_RHCE.loc[indices, 'GRCh38 Alt Allele'] = '2'
df_RHCE.loc[indices, 'GRCh37 Alt Allele'] = '2'
df_RHCE = df_RHCE.astype(str)
df_RHCE.loc[0, 'GRCh38 Coordinates'] = ""
df_RHCE.loc[0, 'GRCh37 Coordinates'] = ""
df_RHCE.loc[0, 'GRCh38 VCF Position'] = ""
df_RHCE.loc[0, 'GRCh37 VCF Position'] = ""
df_RHCE.to_excel("Tables with coords/004_RH (RHCE)_coords.xlsx", index=False)
df_RHCE.to_csv("Tables with coords/004_RH (RHCE)_coords.tsv", sep='\t', index=False)

# Change genotype value in Erythro
indices = df_e[df_e['Nucleotide Change'] == 'NM_020469.2:c.261delG'].index
df_e.loc[indices, 'GRCh38 Coordinates'] = 'NC_000009.12:g.133257521delC'
df_e.loc[indices, 'GRCh37 Coordinates'] = 'NC_000009.11:g.136132908delC'
df_e.loc[indices, 'GRCh38 VCF Position'] = '133257521'
df_e.loc[indices, 'GRCh37 VCF Position'] = '136132908'
df_e.loc[indices, 'GRCh38 Alt Allele'] = '0'
df_e.loc[indices, 'GRCh37 Alt Allele'] = '0'

indices = df_e[df_e['Nucleotide Change'].isin(['NM_003612.3:c.1545A>G', 'NM_002049.3:c.1067G>T'])].index
df_e.loc[indices, 'GRCh38 Alt Allele'] = '2'
df_e.loc[indices, 'GRCh37 Alt Allele'] = '2'

df_e = df_e.astype(str)

df_e.to_excel("Tables with coords/Erythrogene_Table_coords.xlsx", index=False)
df_e.to_csv("Tables with coords/Erythrogene_Table_coords.tsv", sep='\t', index=False)

print("Done")