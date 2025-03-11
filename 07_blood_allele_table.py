import pandas as pd
import numpy as np
import os

# Read files
folder = "Tables with coords"
remove_ISBT_variant_file = 'ISBT_variants_to_remove.tsv'
remove_Erythro_variant_file = 'Erythrogene_variants_to_remove.tsv'
overwrite_ISBT_variant_file = 'ISBT_variants_to_overwrite.tsv'
files = [f for f in os.listdir(folder) if f.endswith('.tsv')]

file_name_lst = [pd.read_csv(os.path.join(folder, file), sep="\t") for file in files]
combined_df = pd.concat(file_name_lst, ignore_index=True)
df1 = pd.read_csv(remove_ISBT_variant_file, sep="\t")
df2 = pd.read_csv(remove_Erythro_variant_file, sep="\t")
to_remove = pd.concat([df1, df2], ignore_index=True)
to_overwrite = pd.read_csv(overwrite_ISBT_variant_file, sep="\t")

removed_NC = combined_df[combined_df["Nucleotide Change"].isin(to_remove["Nucleotide Change"])]
removed_allele = removed_NC.loc[~removed_NC["Allele Name"].isin(["RHCE*01", "KN*01"]), "Allele Name"].unique()
df_filtered = combined_df[~combined_df["Nucleotide Change"].isin(to_remove["Nucleotide Change"])]
df_filtered = df_filtered[~df_filtered["Allele Name"].isin(removed_allele)]

overwrite_NC = df_filtered[df_filtered["Nucleotide Change"].isin(to_overwrite["Old Nucleotide Change"])]
overwrite_allele = overwrite_NC["Allele Name"].unique()

df_merged = pd.merge(df_filtered, to_overwrite[['Old Nucleotide Change', 'New Nucleotide Change', 
                    "GRCh38 Coordinates", "GRCh37 Coordinates", "GRCh38 VCF Position", 
                    "GRCh37 VCF Position", "GRCh38 Alt Allele", "GRCh37 Alt Allele"]],
                    how='left', left_on='Nucleotide Change', right_on='Old Nucleotide Change')

df_merged['Nucleotide Change'] = df_merged['New Nucleotide Change'].fillna(df_merged['Nucleotide Change'])
df_merged['GRCh38 Coordinates_x'] = df_merged['GRCh38 Coordinates_x'].fillna(df_merged['GRCh38 Coordinates_y'])
df_merged['GRCh37 Coordinates_x'] = df_merged['GRCh37 Coordinates_x'].fillna(df_merged['GRCh37 Coordinates_y'])
df_merged['GRCh38 VCF Position_x'] = df_merged['GRCh38 VCF Position_x'].fillna(df_merged['GRCh38 VCF Position_y'])
df_merged['GRCh37 VCF Position_x'] = df_merged['GRCh37 VCF Position_x'].fillna(df_merged['GRCh37 VCF Position_y'])
df_merged['GRCh38 Alt Allele_x'] = np.where(df_merged['GRCh38 Alt Allele_y'].notna(), df_merged['GRCh38 Alt Allele_y'], df_merged['GRCh38 Alt Allele_x'])
df_merged['GRCh37 Alt Allele_x'] = np.where(df_merged['GRCh37 Alt Allele_y'].notna(), df_merged['GRCh37 Alt Allele_y'], df_merged['GRCh37 Alt Allele_x'])

df_merged = df_merged.drop(columns=['Old Nucleotide Change', 'New Nucleotide Change', 'GRCh38 Coordinates_y',
                                    'GRCh37 Coordinates_y', 'GRCh38 VCF Position_y', 'GRCh37 VCF Position_y',
                                    'GRCh38 Alt Allele_y', 'GRCh37 Alt Allele_y'])

df_merged.rename(columns={'GRCh38 Coordinates_x': 'GRCh38 Coordinates', 
                          'GRCh37 Coordinates_x': 'GRCh37 Coordinates',
                          'GRCh38 VCF Position_x': 'GRCh38 VCF Position',
                          'GRCh37 VCF Position_x': 'GRCh37 VCF Position',
                          'GRCh38 Alt Allele_x': 'GRCh38 Alt Allele',
                          'GRCh37 Alt Allele_x': 'GRCh37 Alt Allele'}, inplace=True)

KN_remove = ['NM_000573.3:c.4768A', 'NM_000573.3:c.4801A', 'NM_000573.3:c.4223C',
             'NM_000573.3:c.4828T', 'NM_000573.3:c.4843A', 'NM_000573.3:c.3623A']

df_merged = df_merged[~((df_merged['Nucleotide Change'].isin(KN_remove)) & 
                            (df_merged['Allele Name'] == 'KN*01'))]

df_merged = df_merged.fillna("-")

df_merged['GRCh38 VCF Position'] = df_merged['GRCh38 VCF Position'].apply(
    lambda x: int(x) if isinstance(x, float) else x)
df_merged['GRCh37 VCF Position'] = df_merged['GRCh37 VCF Position'].apply(
    lambda x: int(x) if isinstance(x, float) else x)
df_merged['GRCh38 Alt Allele'] = df_merged['GRCh38 Alt Allele'].astype(int)
df_merged['GRCh37 Alt Allele'] = df_merged['GRCh37 Alt Allele'].astype(int)

df_collapsed = df_merged.groupby('Allele Name', as_index=False, sort=False).agg({
    'Blood System': 'first',
    'Phenotype': 'first',
    'Gene': 'first',
    'Chromosome': 'first',
    'Nucleotide Change': lambda x: ';'.join(x), 
    'GRCh38 Coordinates': lambda x: ';'.join(x.astype(str)),  
    'GRCh37 Coordinates': lambda x: ';'.join(x.astype(str)),  
    'GRCh38 VCF Position': lambda x: ';'.join(x.astype(str)),  
    'GRCh37 VCF Position': lambda x: ';'.join(x.astype(str)),  
    'GRCh38 Alt Allele': lambda x: ';'.join(x.astype(str)),  
    'GRCh37 Alt Allele': lambda x: ';'.join(x.astype(str))  
})

df_final = df_collapsed.drop_duplicates(subset=['Gene', 'Nucleotide Change'], keep='first')

df_final= df_final[['Blood System', 'Phenotype', 'Allele Name', 'Gene', 'Chromosome', 
                             'Nucleotide Change', 'GRCh38 Coordinates', 'GRCh37 Coordinates', 
                             'GRCh38 VCF Position', 'GRCh37 VCF Position', 'GRCh38 Alt Allele', 
                             'GRCh37 Alt Allele']]

df_final.to_excel('Blood_Allele_Table.xlsx', index=False)
df_final.to_csv('Blood_Allele_Table.tsv', sep='\t', index=False)

columns_to_explode = ['Nucleotide Change', 'GRCh38 Coordinates', 'GRCh37 Coordinates', 
                      'GRCh38 VCF Position', 'GRCh37 VCF Position', 'GRCh38 Alt Allele', 'GRCh37 Alt Allele']
df_final[columns_to_explode] = df_final[columns_to_explode].apply(lambda col: col.str.split(';'))
df_final = df_final.explode(columns_to_explode, ignore_index=True)

df_final.to_excel('Blood_Allele_Table_(Separated).xlsx', index=False)
df_final.to_csv('Blood_Allele_Table_(Separated).tsv', sep='\t', index=False)

print("Done")