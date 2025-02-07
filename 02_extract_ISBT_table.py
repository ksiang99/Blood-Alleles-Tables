# Files that require manual copy and paste
# All tables: MNS, H
# Some tables: RHD, KEL, LE, XG, KLF1
# Misc cleaning: RHD, LE, DO, RHAG, JR, EMM

import tabula
import pandas as pd
import os
import re

def clean_df(df):
    try:
        # Remove potential empty rows and columns during extraction
        df = df.dropna(how='all')
        desired_num_columns = len(extracted_df.columns)
        df = df[df.notnull().sum(axis=1) >= desired_num_columns - 2]
        df = df.dropna(axis=1,how='all')

        # Clean up 'Nucleotide Change' column
        # Each phenotype may be associated with > 1 SNP. Split the string such that one row has only 1 nucleotide change.
        df['Nucleotide Change'] = df['Nucleotide Change'].fillna('-')
        df['Nucleotide Change'] = df['Nucleotide Change'].str.split(';')
        df = df.explode('Nucleotide Change')
        df['Nucleotide Change'] = df['Nucleotide Change'].str.split(r'(?:\s|^)(?=\(?c\.)')
        df = df.explode('Nucleotide Change')
        df['Nucleotide Change'] = df['Nucleotide Change'].str.split(r'(exons 2-3 D)')
        df = df.explode('Nucleotide Change')
        df['Nucleotide Change'] = df['Nucleotide Change'].str.split(r'(RHD exon 6-10)')
        df = df.explode('Nucleotide Change')
        df['Nucleotide Change'] = df['Nucleotide Change'].replace(r'[\r\n]+', '', regex=True).str.strip()
        df['Nucleotide Change'] = df['Nucleotide Change'].replace(r'\*+\s*see.*', '', regex=True).str.strip()
        df['Nucleotide Change'] = df['Nucleotide Change'].replace(r'\(.*', '', regex=True).str.strip()
        df = df[~df['Nucleotide Change'].isin(['', 'change', 'Reference nucleotides', 'Obsolete'])]
        
        # Add "c." to nucleotide change that start with a number if "c." is not present
        df['Nucleotide Change'] = df['Nucleotide Change'].str.strip()
        df['Nucleotide Change'] = df['Nucleotide Change'].apply(
        lambda x: 'c.' + x if x[0].isdigit() and not x.startswith('c.') else x)

        # Clean up 'Phenotype' and 'Allele Name' columns
        # Fill the empty rows of "Phenotype" and "Allele Name" columns with the previous row values
        df = df[~df['Phenotype'].str.contains('Null phenotypes|Null alleles|Del defined|RHD weak D|Weak phenotypes|Weak FY|Weak JK|Altered|Mod|Activating|unconfirmed|Null Phenotypes|Null phenotype|\*Obsolete\*', na=False)]
        df = df[~(df['Phenotype'].str.contains('partial D', na=False) & ~df['Phenotype'].str.contains(r'(weak partial D|Weak partial D)', na=False))]
        df["Phenotype"] = df["Phenotype"].apply(clean_phenotype)
        df["Allele Name"] = df["Allele Name"].replace(r'\s*or\s*.*|†.*|‡.*|\(.*', '', regex=True).str.strip()
        df['Allele Name'] = df['Allele Name'].replace(r'[\r\n]+', ' ', regex=True).str.strip()
        df = df[~df['Allele Name'].str.contains('see', case=False, na=False)]
        df = df[~df['Allele Name'].str.contains('Allele Name', case=False, na=False)]
        df = df[~df['Allele Name'].str.contains('obsolete', case=False, na=False)]
        df['Allele Name'] = df['Allele Name'].replace(r'\ssee RHCE.*', '', regex=True)
        df['Phenotype'] = df['Phenotype'].ffill()
        df['Allele Name'] = df['Allele Name'].ffill()

        # Get Gene from Allele Name
        df['Gene'] = df['Allele Name'].str.split('*').str[0]
        df = df.reset_index(drop=True)
        df.loc[df['Gene'] == df['Allele Name'], 'Gene'] = df['Allele Name'].str.split('.').str[0]

        # Add HGVS notation in front of "c." in "Nucleotide Change" column
        df['Transcript'] = df['Gene'].map(HGVS_df.set_index('Allele Name')['Transcript'])
        df['Nucleotide Change'] = df.apply(
            lambda row: row['Transcript'] + ":" + row['Nucleotide Change'] if pd.notna(row['Transcript']) and row['Nucleotide Change'].startswith('c') else row['Nucleotide Change'],
            axis=1)
        df.drop(columns='Transcript', inplace=True)

        # Fill blood group, chromosome column, add HGVS genome notation columns 
        df['Blood System'] = df["Gene"].map(HGVS_df.set_index('Allele Name')['Blood Group'])
        df['Chromosome'] = df['Gene'].map(HGVS_df.set_index('Allele Name')['Chromosome'])
        df['Gene'] = df["Gene"].map(HGVS_df.set_index('Allele Name')['Gene'])
        
        return df
    
    except Exception as e:
        print(f"Skipping file due to error: {e}")
        pass

def clean_phenotype(value):
    if isinstance(value, str):
        value = value.strip()
        value = re.sub(r' {2,}(?=\s*or)', ' ', value)
        value = re.sub(r'(?<=or)(?=\S)', ' ', value)
        value = re.sub(r'(?<=or)\s+', ' ', value) 
        value = re.sub(r'comment.*', '', value)
        value = re.sub(r'(?<=Active).*', '', value)
        value = re.sub(r'erythroid.*', '', value)
        value = re.sub(r'\(i.*', '', value)
        value = re.sub(r'Variant.*', '', value)
        value = re.sub(r',\s*also.*', '', value)
        value = re.sub(r'[\r\n]+', ' ', value)
        value = re.sub(r'\(likely the.*', '', value)
        value = value.strip()

    return value

## Main Script##
group_name = ["ABO", "MNS", "P1PK", "RH (RHCE)", "RH (RHD)", "LU", "KEL", "LE", "FY", "JK", "DI",
               "YT", "XG", "SC", "DO", "CO", "LW", "CHRG", "H", "KX", "GE", "CROM", "KN",
                "IN", "OK", "RAPH", "JMH", "I", "GLOB", "GIL", "RHAG", "FORS", "JR", "LAN",
                "VEL", "CD59", "AUG", "KANNO", "SID", "CTL2", "PEL", "MAM", "EMM", "ABCC1",
                "ER", "GATA1", "KLF1"]

# Read files
save_folder = "Extracted Tables"
HGVS_df = pd.read_csv('HGVS Notation.txt', sep='\t')
pdf_folder = "ISBT Blood Alleles Table"
pdf_files = [f for f in os.listdir(pdf_folder) if f.endswith('.pdf')]

# Extract tables in ISBT pdf
for i, pdf_file in enumerate(pdf_files):
    
    # Skip ISBT pdf if extraction has been done
    if (i == 3 or i == 4) and os.path.exists(f"{save_folder}/004_{group_name[i]}.tsv"):
        continue
    elif i > 9 and os.path.exists(f"{save_folder}/0{i}_{group_name[i]}.tsv"):
        continue
    elif i > 4 and os.path.exists(f"{save_folder}/00{i}_{group_name[i]}.tsv"):
        continue
    elif os.path.exists(f"{save_folder}/00{i+1}_{group_name[i]}.tsv"):
        continue

    pdf_path = os.path.join(pdf_folder, pdf_file)
    tables = tabula.read_pdf(pdf_path, pages='all', multiple_tables=True)
    if tables:
        concatenated_df = pd.concat(tables, ignore_index=True)

    # Extract Phenotype, Allele Name and Nucleotide Change columns from table
    try:
        pheno_col = concatenated_df.iloc[:, 0]
        allele_col = concatenated_df.iloc[:, 1]
        nc_col = concatenated_df.iloc[:, 2]

    except IndexError:
        try:
            pheno_col = concatenated_df["Phenotype"]
            allele_col = concatenated_df["Allele name"]
            nc_col = concatenated_df["Nucleotide change"]

        except KeyError:
            pheno_col = pd.Series(dtype='object')
            allele_col = pd.Series(dtype='object')
            nc_col = pd.Series(dtype='object')

    # Order the columns in extracted_df
    extracted_df = pd.DataFrame({
        'Blood System' : "",
        'Phenotype': pheno_col,
        'Allele Name': allele_col,
        'Gene' : "",
        'Chromosome' : "",
        'Nucleotide Change': nc_col,
        'GRCh38 Coordinates' : "",
        'GRCh37 Coordinates' : "",
        'GRCh38 VCF Position' : "",
        'GRCh37 VCF Position' : ""
    })
    
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)

    # Save dataframe in excel and tsv format
    try:
        cleaned_df = clean_df(extracted_df)
        if i == 3 or i == 4:
            cleaned_df.to_excel(f"{save_folder}/004_{group_name[i]}.xlsx", index=False)
            cleaned_df.to_csv(f"{save_folder}/004_{group_name[i]}.tsv", sep='\t', index=False)
        
        elif i > 9:
            cleaned_df.to_excel(f"{save_folder}/0{i}_{group_name[i]}.xlsx", index=False)
            cleaned_df.to_csv(f"{save_folder}/0{i}_{group_name[i]}.tsv", sep='\t', index=False)        
        elif i > 4:
            cleaned_df.to_excel(f"{save_folder}/00{i}_{group_name[i]}.xlsx", index=False)
            cleaned_df.to_csv(f"{save_folder}/00{i}_{group_name[i]}.tsv", sep='\t', index=False)
        else:
            cleaned_df.to_excel(f"{save_folder}/00{i+1}_{group_name[i]}.xlsx", index=False)
            cleaned_df.to_csv(f"{save_folder}/00{i+1}_{group_name[i]}.tsv", sep='\t', index=False)

        print(f'Processed {pdf_file}')

    except Exception as e:
        print(f"Error processing {pdf_file}: {e}")
        continue

print("Done")