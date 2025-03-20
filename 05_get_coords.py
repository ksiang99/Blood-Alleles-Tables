from requests.exceptions import RequestException, ConnectionError
import pandas as pd
import numpy as np
import requests
import time
import os
import re

base_url = "https://rest.variantvalidator.org/VariantValidator/variantvalidator"

# Function to fetch variant data from VariantValidator API
def fetch_variant_data(genome_build, variant_description, select_transcripts, retries=3, delay=2, timeout=5):
    url = f"{base_url}/{genome_build}/{variant_description}/{select_transcripts}"
    attempt = 0
    
    # Provide 3 tries to request API in case of failure
    while attempt < retries:
        try:
            response = requests.get(url, headers={"content-type": "application/json"}, timeout=timeout)
            if response.status_code == 200:
                return response.json()
            else:
                return {"error": response.status_code, "message": response.text}
        except (RequestException, ConnectionError) as e:
            attempt += 1
            print(f"Request failed (Attempt {attempt}/{retries}). Error: {e}")
            if attempt < retries:
                print(f"Retrying in {delay} seconds...")
                time.sleep(delay)
            else:
                return {"error": "Request failed after multiple attempts", "message": str(e)}

# Function to extract position based on hgvs_genomic_description
def extract_position(data, NC_descript, GRCh38_HGVS, GRCh37_HGVS):
    ref_build = ['grch37', 'grch38']
    
    # Clean up nucleotide changes text
    NC_descript = NC_descript.replace(" ", "")
    NC_descript = NC_descript.replace("−", "-")
    NC_descript = re.sub(r"([a-z])>([a-z])", lambda match: match.group(1).upper() + ">" + match.group(2).upper(), NC_descript)
    print(NC_descript)
    
    try:
        # Handle case where nucleotide change given is in genome instead of transcript
        if "NC_" in NC_descript:
            primary_assembly_loci = data["intergenic_variant_1"]["primary_assembly_loci"]
            grch37_coord = primary_assembly_loci[ref_build[0]]["hgvs_genomic_description"]
            grch38_coord = primary_assembly_loci[ref_build[1]]["hgvs_genomic_description"]
            if GRCh37_HGVS in grch37_coord:
                grch37_pos = str(int(primary_assembly_loci[ref_build[0]]["vcf"]["pos"]))
                grch37_coord = (
                    grch37_coord.split('g.')[0] + "g." + 
                    str(int(grch37_pos) + 1) + "_" + 
                    str(int(grch37_pos) + len(primary_assembly_loci[ref_build[0]]["vcf"]["ref"]) - 1) + 
                    "del" + primary_assembly_loci[ref_build[0]]["vcf"]["ref"][1:10]
                )
            if GRCh38_HGVS in grch38_coord:
                grch38_pos = str(int(primary_assembly_loci[ref_build[1]]["vcf"]["pos"]))
                grch38_coord = (
                grch38_coord.split('g.')[0] + "g." + 
                str(int(grch38_pos) + 1) + "_" + 
                str(int(grch38_pos) + len(primary_assembly_loci[ref_build[1]]["vcf"]["ref"]) - 1) + 
                "del" + primary_assembly_loci[ref_build[1]]["vcf"]["ref"][1:10]
                )
            return grch38_coord, grch38_pos, grch37_coord, grch37_pos

        # Handle delins case
        elif "delins" in NC_descript:
            pass
        
        # Handle deletion case
        elif "del" in NC_descript:
            # Check if element before del is a digit. If not, skip the nucleotide change
            check_digit = re.search(r"c\.(\d+)(?:[-_+]\d+)*del", NC_descript)
            if not check_digit:
                print(1)
                return None, None, None, None
           
            # Determine if deletion is single base or multiple bases and provide HGVS compliant
            # nucleotide change to VariantValidator.
            # Eg. NM_016124.6:c.684delGAG will be NM_016124.6:c.684_686delGAG
            match = re.search(r"c\.(\d+)(?:[-_+]?\d*)del([A-Za-z]+)", NC_descript)
            if match:
                NC_sequence = match.group(2)
                if len(NC_sequence) > 1:
                    # Handle cases where '-' is used to separate deletion position instead of '_'
                    NC_descript = NC_descript.replace("-", "_")
                    start = int(match.group(1))
                    end = start + len(NC_sequence) - 1
                    NC_descript = NC_descript.split('c.')[0] + "c." + str(start) + "_" + str(end) + "del"
           
            NC_descript = NC_descript.split("del")[0] + "del"
            NC_descript.strip()
            print(NC_descript)
            primary_assembly_loci = data[NC_descript]["primary_assembly_loci"]
            grch37_coord = primary_assembly_loci[ref_build[0]]["hgvs_genomic_description"]
            grch38_coord = primary_assembly_loci[ref_build[1]]["hgvs_genomic_description"]
            us_count_grch37 = grch37_coord.count("_")
            us_count_grch38 = grch38_coord.count("_")
            
            if GRCh37_HGVS in grch37_coord:
                grch37_pos = str(int(primary_assembly_loci[ref_build[0]]["vcf"]["pos"]))
                if us_count_grch37 == 2:
                    grch37_coord = (
                    grch37_coord.split('g.')[0] + "g." + 
                    str(int(grch37_pos) + 1) + "_" + 
                    str(int(grch37_pos) + len(primary_assembly_loci[ref_build[0]]["vcf"]["ref"]) - 1) +
                    "del" + primary_assembly_loci[ref_build[0]]["vcf"]["ref"][1:10]
                    )
                else:
                    grch37_coord = (
                    grch37_coord.split('g.')[0] + "g." + 
                    str(int(grch37_pos) + 1) + "del" +
                    primary_assembly_loci[ref_build[0]]["vcf"]["ref"][1]
                    )
            
            if GRCh38_HGVS in grch38_coord:
                grch38_pos = str(int(primary_assembly_loci[ref_build[1]]["vcf"]["pos"]))
                if us_count_grch38 == 2:
                    grch38_coord = (
                    grch38_coord.split('g.')[0] + "g." + 
                    str(int(grch38_pos) + 1) + "_" + 
                    str(int(grch38_pos) + len(primary_assembly_loci[ref_build[1]]["vcf"]["ref"]) - 1) + 
                    "del" + primary_assembly_loci[ref_build[1]]["vcf"]["ref"][1:10]
                    )
                else:
                    grch38_coord = (
                    grch38_coord.split('g.')[0] + "g." + 
                    str(int(grch38_pos) + 1) + "del" +
                    primary_assembly_loci[ref_build[1]]["vcf"]["ref"][1]
                    )
            # if "=" in grch37_coord:
            #     grch37_coord = primary_assembly_loci[ref_build[0]]["hgvs_genomic_description"] + primary_assembly_loci[ref_build[0]]["vcf"]["ref"]
            # if "=" in grch38_coord:
            #     grch38_coord = primary_assembly_loci[ref_build[1]]["hgvs_genomic_description"] + primary_assembly_loci[ref_build[1]]["vcf"]["ref"]
            return grch38_coord, grch38_pos, grch37_coord, grch37_pos

        elif "ins" in NC_descript:
            # Check if element before ins is a digit. If not, skip the nucleotide change
            check_digit = re.search(r"c\.(\d+)[-_]?\d*ins", NC_descript)
            if not check_digit:
                return None, None, None, None

            # Check if nucleotide change is in HGVS compliant. If not, provide HGVS compliant nucleotide change
            # Eg. NM_020485.8:c.93insT will be NM_020485.8:c.93_94insT
            if NC_descript.count("_") != 2:
                match = re.search(r"c\.(\d+)ins([A-Za-z]+)", NC_descript)
                if match:
                    start = int(match.group(1))
                    NC_sequence = match.group(2)  
                    end = start + 1
                    NC_descript = NC_descript.split('ins')[0] + "_" + str(end) + "ins" + NC_sequence
                    NC_descript.strip()
            
            # Handle cases where API returns "dup" instead of "ins"
            # Eg. NM_020485.8:c.94dup instead of NM_020485.8:c.93_94ins
            try:
                primary_assembly_loci = data[NC_descript]["primary_assembly_loci"]
            
            except KeyError:
                match = re.search(r"(NM_\d+\.\d+:c\.\d+)_\d+ins(\w+)", NC_descript)
                if match:
                    first_pos = match.group(1)
                    second_pos = NC_descript.split('_')[-1].split('ins')[0]
                    NC_descript = first_pos.split(':')[0] + ":" + 'c.' + second_pos + 'dup'
                    NC_descript.strip()
            
            # Handle cases where wrong position is provided 
            # (Sometimes it takes the start position instead of end position)
            # Eg. NM_020485.8:c.93dup instead of NM_020485.8:c.94dup
            try:
                primary_assembly_loci = data[NC_descript]["primary_assembly_loci"]
           
            except KeyError:
                NC_descript = first_pos + 'dup'
                NC_descript.strip()
                primary_assembly_loci = data[NC_descript]["primary_assembly_loci"]
            
            # If nucleotide change is not found, return None
            if primary_assembly_loci == None:
                return None, None, None, None
            
            grch37_coord = primary_assembly_loci[ref_build[0]]["hgvs_genomic_description"]
            grch38_coord = primary_assembly_loci[ref_build[1]]["hgvs_genomic_description"]
            
            if GRCh37_HGVS in grch37_coord:
                if "dup" in grch37_coord:
                    grch37_coord = grch37_coord + primary_assembly_loci[ref_build[0]]["vcf"]["alt"][1:10]
                    grch37_coord = grch37_coord.replace("dup", "ins")
                grch37_pos = str(int(primary_assembly_loci[ref_build[0]]["vcf"]["pos"]))
            
            if GRCh38_HGVS in grch38_coord:
                if "dup" in grch38_coord:
                    grch38_coord = grch38_coord + primary_assembly_loci[ref_build[1]]["vcf"]["alt"][1:10]
                    grch38_coord = grch38_coord.replace("dup", "ins")
                grch38_pos = str(int(primary_assembly_loci[ref_build[1]]["vcf"]["pos"]))
            
            return grch38_coord, grch38_pos, grch37_coord, grch37_pos

        # Handle dup case
        elif "dup" in NC_descript:
            NC_descript = NC_descript.split("dup")[0] + "dup"
            primary_assembly_loci = data[NC_descript]["primary_assembly_loci"]
            grch37_coord = primary_assembly_loci[ref_build[0]]["hgvs_genomic_description"]
            grch38_coord = primary_assembly_loci[ref_build[1]]["hgvs_genomic_description"]
            
            if GRCh37_HGVS in grch37_coord:
                grch37_coord = grch37_coord + primary_assembly_loci[ref_build[0]]["vcf"]["alt"][1:]
                grch37_pos = str(int(primary_assembly_loci[ref_build[0]]["vcf"]["pos"]))
            
            if GRCh38_HGVS in grch38_coord:
                grch38_coord = grch38_coord + primary_assembly_loci[ref_build[1]]["vcf"]["alt"][1:]
                grch38_pos = str(int(primary_assembly_loci[ref_build[1]]["vcf"]["pos"]))
            
            return grch38_coord, grch38_pos, grch37_coord, grch37_pos
        
        # Handle delin and SNPs case
        primary_assembly_loci = data[NC_descript]["primary_assembly_loci"]
        grch37_coord = primary_assembly_loci[ref_build[0]]["hgvs_genomic_description"]
        grch38_coord = primary_assembly_loci[ref_build[1]]["hgvs_genomic_description"]
        
        if GRCh37_HGVS in grch37_coord:
            grch37_pos = str(int(primary_assembly_loci[ref_build[0]]["vcf"]["pos"]))
        
        if GRCh38_HGVS in grch38_coord:
            grch38_pos = str(int(primary_assembly_loci[ref_build[1]]["vcf"]["pos"]))
        
        if "=" in grch37_coord:
                grch37_coord = grch37_coord + primary_assembly_loci[ref_build[0]]["vcf"]["alt"]
        
        if "=" in grch38_coord:
                grch38_coord = grch38_coord + primary_assembly_loci[ref_build[1]]["vcf"]["alt"]
        
        return grch38_coord, grch38_pos, grch37_coord, grch37_pos
        
    except KeyError:
        return None, None, None, None

# Function to fetch GRCh38 and GRCh37 data for each variant
def fetch_variant_for_row(row, retries=3, delay=2):
    attempt = 0
    while attempt < retries:
        data = fetch_variant_data(
            genome_build="GRCh37",
            variant_description=row['Nucleotide Change'],
            select_transcripts="refseq_select"
        )

        # Extract positions for both GRCh38 and GRCh37
        grch38_coord, grch38_pos, grch37_coord, grch37_pos = extract_position(data, row["Nucleotide Change"], row["GRCh38 HGVS"], row["GRCh37 HGVS"])
        
        if None not in (grch38_coord, grch38_pos, grch37_coord, grch37_pos):
            print(grch38_coord, grch38_pos, grch37_coord, grch37_pos)
            return grch38_coord, grch37_coord, grch38_pos, grch37_pos

        else:
            attempt += 1
            time.sleep(delay)

    return None, None, None, None

## Main Script ##
HGVS_df = pd.read_csv("HGVS_Notation.txt", delimiter="\t")
save_folder = "Tables with coords"
folder = "Extracted Tables"
files = [f for f in os.listdir(folder) if f.endswith('.tsv')]

if not os.path.exists(save_folder):
    os.makedirs(save_folder)

for i, file in enumerate(files):
    base_name = os.path.splitext(file)[0]
    
    if os.path.exists(f"{save_folder}/{base_name}_coords.tsv"):
        continue

    path = os.path.join(folder, file)
    df = pd.read_csv(path, sep='\t')
    print(f"Processing {file}")

    # Mapping genes to their corresponding HGVS notation
    df['GRCh38 HGVS'] = df['Gene'].map(HGVS_df.set_index('Gene')['GRCh38_HGVS'])
    df['GRCh37 HGVS'] = df['Gene'].map(HGVS_df.set_index('Gene')['GRCh37_HGVS'])

    # Drop duplicate nucleotide changes to reduce run time
    temp_df = df.drop_duplicates(subset=["Nucleotide Change", "Gene"])

    # Retrieve coordinates using VariantValidator API
    temp_df['GRCh38 Coordinates'], temp_df['GRCh37 Coordinates'], temp_df['GRCh38 VCF Position'], temp_df['GRCh37 VCF Position'] = zip(*temp_df.apply(fetch_variant_for_row, axis=1))

    # Merge results from VariantValidator API with dataframe
    merged_df = pd.merge(
        df, temp_df[["Nucleotide Change", "GRCh38 Coordinates", "GRCh37 Coordinates", "GRCh38 VCF Position", "GRCh37 VCF Position"]],
        on="Nucleotide Change",
        how="left")
    merged_df = merged_df.drop_duplicates().reset_index(drop=True)
    df["GRCh38 Coordinates"] = merged_df["GRCh38 Coordinates_y"]
    df["GRCh37 Coordinates"] = merged_df["GRCh37 Coordinates_y"]
    df["GRCh38 VCF Position"] = merged_df["GRCh38 VCF Position_y"]
    df["GRCh37 VCF Position"] = merged_df["GRCh37 VCF Position_y"]
    df.drop(['GRCh38 HGVS', 'GRCh37 HGVS'], axis=1, inplace=True)

    # Create 2 columns "GRCh38 Alt Allele" and "GRch37 Alt Allele" to indicate changes in base.
    # 0 refers to no change in base (Same as ref genome), 1 refers to change in base
    df['GRCh38 Alt Allele'] = np.where(df['GRCh38 Coordinates'].isna() | (df['GRCh38 Coordinates'] == '') | df['GRCh38 Coordinates'].str.contains('='), 0, 1)
    df['GRCh37 Alt Allele'] = np.where(df['GRCh37 Coordinates'].isna() | (df['GRCh37 Coordinates'] == '') | df['GRCh37 Coordinates'].str.contains('='), 0, 1)
    
    # Save dataframe onto excel and text file
    df.to_excel(f'{save_folder}/{base_name}_coords.xlsx', index=False)
    df.to_csv(f'{save_folder}/{base_name}_coords.tsv', sep='\t', index=False)
    print(f"Processed {file}")

print("Done")