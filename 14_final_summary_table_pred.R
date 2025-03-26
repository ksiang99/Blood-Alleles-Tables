rm(list = ls())

setwd("/home/svu/e0726996/Blood-type-GWAS/")

# Set random seed
set.seed(0)

library(dplyr)
library(stringr)

# Set paths
PATH_ALLELE_TABLE <- "./R/Blood_Allele_Table.tsv"
df_allele_table <- read.delim(PATH_ALLELE_TABLE, header = TRUE)
df_allele_table$Pheno_Allele <- paste0(df_allele_table$Phenotype, "/", df_allele_table$Allele.Name)
PREDICT_PATH <- "./R/results/predict"
pred_files <- list.files(PREDICT_PATH, pattern = "summary_table\\.Rdata$", full.names = TRUE)
pred_files <- pred_files[!grepl("final_summary_table", pred_files)]

# Function to sort the alleles based on number of nucleotide changes (Largest to smallest)
process_hier <- function(idx) {
    name <- df_allele_table$Pheno_Allele[idx]
    if (!any(grepl("KN\\:", name))) {
        name <- sub(".* or ([^|]+)$", "\\1", name)
    }
    nc <- df_allele_table$Nucleotide.Change[idx]
    
    # Count number of nucleotide changes
    num_nc <- str_count(nc, ";") + 1
    num_nc[nc == "-"] <- 0
    
    # numbering <- as.numeric(str_extract(name, "\\d+"))
    
    # Sorting based on criteria.
    # 1. Reference allele is placed last.
    # 2. Alleles sorted by descending number of nucleotide changes.
    # 3. If number of nucleotide changes are the same, Erythrogene alleles are placed after ISBT alleles
    # 4. Keep ISBT alleles in original order
    sorted_indices <- order(nc == "-", -num_nc, str_detect(name, "-/\\*1000G only\\*.*"), seq_along(name))
    hier <- name[sorted_indices]
    
    return(hier)
}

# Function to remove reference if there are more than 1 potential phenotypes
remove_ref <- function(row, ref) {
    substrings <- strsplit(row, " \\| ")[[1]]

    if (length(substrings) > 1) {
        # Handle H blood group cases as there are 2 reference allele FUT1*01 and FUT2*01
        if (any(grepl("FUT1\\*", ref))) {
            if (sum(grepl("FUT1\\*", substrings)) == 1) {
                return(substrings)
            }
        }
        if (any(grepl("FUT2\\*", ref))) {
            if (sum(grepl("FUT2\\*", substrings)) == 1) {
                return(substrings)
            }
        }
        # Handle CHRG blood group cases as there are 2 reference allele C4A*03 and C4B*03
        if (any(grepl("C4A\\*", ref))) {
            if (sum(grepl("C4A\\*", substrings)) == 1) {
                return(substrings)
            }
        }
        if (any(grepl("C4B\\*", ref))) {
            if (sum(grepl("C4B\\*", substrings)) == 1) {
                return(substrings)
            }
        }

        # Handle MNS blood group cases as there are 2 reference allele GYPA*01 and GYPB*04.05
        if (any(grepl("GYPA\\*", ref))) {
            if (sum(grepl("GYPA\\*", substrings)) == 1) {
                return(substrings)
            }
        }
        if (any(grepl("GYPB\\*", ref))) {
            if (sum(grepl("GYPB\\*", substrings)) == 1) {
                return(substrings)
            }
        }

        # Handle RH blood group cases as there are 2 reference allele RHCE*01 and RHD*01
        if (any(grepl("RHCE\\*", ref))) {
            if (sum(grepl("RHCE\\*", substrings)) == 1) {
                return(substrings)
            }
        }
        if (any(grepl("RHD\\*", ref))) {
            if (sum(grepl("RHD\\*", substrings)) == 1) {
                return(substrings)
            }
        }

    }

    else if (length(substrings) == 1) {
        return(substrings)
    }
    
    substrings <- substrings[substrings != ref]
    substrings <- substrings[substrings != ""]
    return(substrings)
}


# Function to find generalised phenotype name
find_key <- function(dict, value) {
  keys <- names(dict)
  for (key in keys) {
    if (key != "hier" && value %in% dict[[key]]) {
      return(key)
    }
  }
}

# Function to keep 1 phenotype based on hierarchy
final_pheno <- function(row, dict) {
    substrings <- unlist(strsplit(row, " \\| "))
    
    for (pheno in dict[["hier"]]) {
        if (any(grepl(pheno, substrings, fixed = TRUE))) {
            pattern <- paste0(dict[["allele_name"]], collapse = "|")
            pheno <- sub(pattern, "", pheno)
            key <- find_key(dict, pheno)
            result <- paste0(key, ": ", pheno)
            return(result)
        }
    }
    return(row)
}

# Function that uses the functions above
process_group <- function(row, dict) {
    process_single_row <- function(single_row) {
        single_row <- remove_ref(single_row, dict[["ref"]])
        if (!any(grepl("KN\\:", single_row))) {
            single_row <- sub(".* or ([^|]+)$", "\\1", single_row)
        }
        return(final_pheno(single_row, dict))
    }
    return(sapply(row, process_single_row, USE.NAMES = FALSE))
}

# Create named list for each blood group
# ref: reference allele of blood group
# allele_name: pattern used to clean up string at the end
# others: Provide a key to ISBT phenotype name
# TGP: To represent 1000G alleles

ABO <- within(list(
    ref = c("A1/ABO*A1.01"), allele_name = c("/ABO.*", "/\\*1000G.*"),
    "1" = c("A1"), "2" = c("A2"), "3" = c("Aweak", "Ax/Aweak", "Ael"), "4" = c("B"), "5" = c("B3", "Bel"),
    "6" = c("B(A)"), "7" = c("cisAB"), "8" = c("O"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "ABO")
    })

MNS_GYPA <- within(list(
    ref = c("MNS:1 or M+/GYPA*01"), allele_name = c("/GYPA.*", "/\\*1000G.*"),
    "1" = c("M+"), "2" = c("N+"), "3" = c("Mc+"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "GYPA")
    })

MNS_GYPB <- within(list(
    ref = c("MNS:4 or s+/GYPB*04"), allele_name = c("/GYPB.*", "/\\*1000G.*"),
    "1" = c("s+"), "2" = c("S+"), "3" = c("sD+"), "4" = c("Mit+"), "5" = c("S–U+var"),
    "6" = c("Mv+"), "7" = c("S-s-U-"), "8" = c("s+w, U+w"), "9" = c("MNS:4,5 (s+, U+) altered GPB"),
    "10" = c("S+w"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "GYPB")
    })

P1PK <- within(list(
    ref = c("P1+ Pk+/A4GALT*01"), allele_name = c("/A4GALT.*", "/\\*1000G.*"),
    "1" = c("P1+ Pk+"),"2" = c("p"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "A4GALT")
    })

RHCE <- within(list(
    ref = c("RH:4 or c RH:5 or e RH:6 or f (ce)/RHCE*01"),
    allele_name = c("/RHCE.*", "/\\*1000G.*"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "RHCE")
    }
)

RHD <- within(list(
    ref = c("D RH:1/RHD*01"),
    allele_name = c("/RHD.*", "/\\*1000G.*"),
    "1" = c("D RH:1"), "2" = c("DAU0"), "3" = c("DAU3"), "4" = c("Del"), "5" = c("DAU0.01"), 
    "6" = c("DAR3.1 (weak partial D 4.0)"), "7" = c("DIVa RH30+ (Goa+) DIV type 1.0"), "8" = c("DIII type 4"), "9" = c("DV Type 4 RH:23 (Dw+)"),"10" = c("DAU5"),
    "11" = c("DIIIa RH:54 (DAK+)"), "12" = c("RHD(F175L)"), "13" = c("Type 45"), "14" = c("DVII RH:40 (Tar+)"), "15" = c("Type 1"),
    "16" = c("Type 3"),"17" = c("DAU6"), "18" = c("DNB"), "19" = c("Type 28"), "20" = c("Type 66"),
    "21" = c("Type 25"), "22" = c("Type 33"), "23" = c("DIII type 6"), "24" = c("DFV"), "25" = c("DAU0.02"),
    "26" = c("DAU14"), "27" = c("Weak partial Type 15"), "28" = c("DUC2"),
    TGP = c("-")
    ), {
        hier <- process_hier(df_allele_table$Gene == "RHD")
    })


LU <- within(list(
    ref = c("LU:2 or Lu(b+)/LU*02"), allele_name = c("/LU.*", "/\\*1000G.*"),
    "1" = c("Lu(b+)"), "2" = c("Au(a−b+)"), "3" = c("LU:-8,14"), "4" = c("LU:-26, LUBI−"),
    "5" = c("LU:-6,9"), "6" = c("LU:-22, LURC−"), "7" = c("Lu(a+)"), "8" = c("LU:1,19"),
    "9" = c("LU:-16"), "10" = c("LU:-13"), "11" = c("LU:-29, LURA−"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "BCAM")
    })

KEL <- within(list(
    ref = c("KEL:2 or k+/KEL*02"), allele_name = c("/KEL.*", "/\\*1000G.*"),
    "1" = c("k+"), "2" = c("KEL:1weak", "Kmod", "K0 phenotype"), "3" = c("Js(a+b–)"), "4" = c("K+k-"),
    "5" = c("Kp(a+b–c–)"), "6" = c("Ul(a+)"), "7" = c("KEL:-19"), "8" = c("KEL:-18"), "9" = c("KTIM-"), 
    "10" = c("KYO+, KYOR-"), "11" = c("TOU-"), "12" = c("KHUL-,KEAL+"), "13" = c("K0"), 
    "14" = c("KHIZ-, KHOZ+"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "KEL")
    })

LE <- within(list(
    ref = c("FUT3 Active/FUT3*01.01"), allele_name = c("/FUT3.*", "/\\*1000G.*"),
    "1" = c("FUT3 Active"), "2" = c("FUT3 Inactive - Le(a-b-)"), "3" = c("FUT3 Active (Weak)"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "FUT3")
    })

FY <- within(list(
    ref = c("FY:1 or Fy(a+)/FY*01"), allele_name = c("/FY.*", "/\\*1000G.*"),
    "1" = c("Fy(a+)"), "2" = c("Fy(a+w)"), "3" = c("Fy(b+)"), "4" = c("Fy(b+w), Fyx"),
    "5" = c("Fy(a−b−)"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "ACKR1")
    })

JK <- within(list(
    ref = c("JK:1 or Jk(a+)/JK*01"), allele_name = c("/JK.*", "/\\*1000G.*"),
    "1" = c("Jk(a+)"), "2" = c("Jk(a+W)", "Jk(a+w)"), "3" = c("Jk(b+)"),
    "4" = c("Jk(b+W)", "Jk(b+w)"), "5" = c("Jk(a–b–)"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "SLC14A1")
    })   

DI <- within(list(
    ref = c("DI:–1,2 or Di(a–b+)/DI*02"), allele_name = c("/DI.*", "/\\*1000G.*"),
    "1" = c("Di(a–b+)"), "2" = c("Di(a+b–)"), "3" = c("Wr(a+b–)"), "4" = c("Wd(a+)"),
    "5" = c("DI:23"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "SLC4A1")
    })

YT <- within(list(
    ref = c("YT:1,-2 or Yt(a+b-)/YT* 01"), allele_name = c("/YT.*", "/\\*1000G.*"),
    "1" = c("Yt(a+b-)"), "2" = c("Yt(a-b+)"), "3" = c("YTEG-"), "4" = c("YTLI-"),
    "5" = c("YTOT-"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "ACHE")
    })

XG <- within(list(
    ref = c("Xga/XG*01"), allele_name = c("/XG.*", "/\\*1000G.*"),
    "1" = c("Xga"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "XG")
    })

SC <- within(list(
    ref = c("SC:1 or Sc1+/SC*01"), allele_name = c("/SC.*", "/\\*1000G.*"),
    "1" = c("Sc1+"), "2" = c("Sc2+"), "3" = c("STAR–"), "4" = c("SCAC–"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "ERMAP")
    })

DO <- within(list(
    ref = c("DO:2 or Do(b+)/DO*02"), allele_name = c("/DO.*", "/\\*1000G.*"),
    "1" = c("Do(b+)"), "2" = c("Do(a+)"), "3" = c("Hy–"), "4" = c("Jo(a–)"),
    "5" = c("DODE–"), "6" = c("DOLG–"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "ART4")
    })

CO <- within(list(
    ref = c("CO:1 or Co(a+)/CO*01.01"), allele_name = c("/CO.*", "/\\*1000G.*"),
    "1" = c("Co(a+)"), "2" = c("Co(b+)"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "AQP1")
    })

LW <- within(list(
    ref = c("LW:5 or LW(a+)/LW*05"), allele_name = c("/LW.*", "/\\*1000G.*"),
    "1" = c("LW(a+)"), "2" = c("LW(b+)"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "ICAM4")
    })

CHRG_A <- within(list(
    ref = c("Ch–Rg+ or CH:–1,–2,–3,–4,–5,–6 RG:1,2/C4A*03"), allele_name = c("/C4A.*", "/\\*1000G.*"),
    "1" = c("CH:–1,–2,–3,–4,–5,–6 RG:1,2"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "C4A")
    })

CHRG_B <- within(list(
    ref = c("Ch+Rg– or CH:1,2,3,4,5,6 RG:–1,–2/C4B*03"), allele_name = c("/C4B.*", "/\\*1000G.*"),
    "1" = c("CH:1,2,3,4,5,6 RG:–1,–2"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "C4B")
    })

H_FUT1 <- within(list(
    ref = c("H+/FUT1*01"), allele_name = c("/FUT1.*", "/\\*1000G.*"),
    "1" = c("H+"), "2" = c("H+weak"), "3" = c("H–"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "FUT1")
    })

H_FUT2 <- within(list(
    ref = c("H+/FUT2*01"),allele_name = c("/FUT2.*", "/\\*1000G.*"),
    "1" = c("H+"), "2" = c("H+w"), "3" = c("H–"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "FUT2")
    })

KX <- within(list(
    ref = c("XK:1 or Kx+/XK*01"), allele_name = c("/XK.*", "/\\*1000G.*"),
    "1" = c("Kx+"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "XK")
    })

GE <- within(list(
    ref = c("GE:2,3,4/GE*01"), allele_name = c("/GE.*", "/\\*1000G.*"),
    "1" = c("GE:2,3,4"), "2" = c("GEPL–"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "GYPC")
    })

CROM <- within(list(
    ref = c("CROM:1 or Cra+)/CROM*01"), allele_name = c("/CROM.*", "/\\*1000G.*"),
    "1" = c("Cra+)"), "2" = c("Cr(a–)"), "3" = c("Tc(b+)"), "4" = c("WES(a+)"),
    "5" = c("SERF–"), "6" = c("Tc(c+)"), "7" = c("UMC–"), "8" = c("CRUE−"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "CD55")
    })

KN <- within(list(
    ref = c("KN:1 or Kn(a+) KN:3 or McC(a+) KN:4 or Sl1+ KN:5 or Yk(a+) KN:8 or Sl3+ KN:9 or KCAM+ KN:11 or DACY+/KN*01"),
    allele_name = c("/KN.*", "/\\*1000G.*"),
    "1" = c("KN:1 or Kn(a+) KN:3 or McC(a+) KN:4 or Sl1+ KN:5 or Yk(a+) KN:8 or Sl3+ KN:9 or KCAM+ KN:11 or DACY+"),
    "2" = c("KN:–5 or Yk(a–)"), "3" = c("KN:-9 or KCAM- KN:10 or KDAS+"), "4" = c("KN:2 or Kn(a-b+) KN:-9 or KCAM- KN:10 or KDAS+"),
    "5" = c("KN:-4 or Sl1- KN:-3,6 or McC(a-b+) KN:7 or Vil+ KN:-9 or KCAM- KN:10 or KDAS+"), "6" = c("KN:7 or Vil+ KN:-9 or KCAM- KN:10 or KDAS+"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "CR1")
    })

IN <- within(list(
    ref = c("In(a–b+)/IN*02"), allele_name = c("/IN.*", "/\\*1000G.*"),
    "1" = c("In(a–b+)"), "2" = c("In(a+b–)"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "CD44")
    })

OK <- within(list(
    ref = c("OK:1 or Ok(a+)/OK*01.01"), allele_name = c("/OK.*", "/\\*1000G.*"),
    "1" = c("Ok(a+)"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "BSG")
    })

RAPH <- within(list(
    ref = c("RAPH:1 or MER2+/RAPH*01"), allele_name = c("/RAPH.*", "/\\*1000G.*"),
    "1" = c("MER2+"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "CD151")
    })

JMH <- within(list(
    ref = c("JMH:1 or JMH+/JMH*01"), allele_name = c("/JMH.*", "/\\*1000G.*"),
    "1" = c("JMH+"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "SEMA7A")
    })

I <- within(list(
    ref = c("I:1 or I+/GCNT2*01"), allele_name = c("/GCNT2.*", "/\\*1000G.*"),
    "1" = c("I+"), "2" = c("I+W"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "GCNT2")
    })

GLOB <- within(list(
    ref = c("GLOB:1 (P+)/GLOB*01"), allele_name = c("/GLOB.*", "/\\*1000G.*"),
    "1" = c("GLOB:1 (P+)"), "2" = c("GLOB:–1 (P–)"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "B3GALNT1")
    })

GIL <- within(list(
    ref = c("GIL:1 or GIL+/GIL*01", "/\\*1000G.*"), allele_name = c("/GIL.*", "/\\*1000G.*"),
    "1" = c("GIL+"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "AQP3")
    })

RHAG <- within(list(
    ref = c("RHAG:1 or Duclos+/RHAG*01"), allele_name = c("/RHAG.*", "/\\*1000G.*"),
    "1" = c("Duclos+"), "2" = c("DSLK−, Kg+"), "3" = c("Rhmod"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "RHAG")
    })

FORS <- within(list(
    ref = c("FORS:-1 (FORS-)/GBGT1*01N.01"), allele_name = c("/GBGT1.*", "/\\*1000G.*"),
    "1" = c("FORS:-1 (FORS-)"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "GBGT1")
    })

JR <- within(list(
    ref = c("Jr(a+)/ABCG2*01"), allele_name = c("/ABCG2.*", "/N/A.*", "/\\*1000G.*"),
    "1" = c("Jr(a+)"), "2" = c("Jr(a+w)"), "3" = c("Unclear Jra phnotype"), "4" = c("Jr(a−)"), TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "ABCG2")
    })

LAN <- within(list(
    ref = c("Lan+/ABCB6*01"), allele_name = c("/ABCB6.*", "/\\*1000G.*"),
    "1" = c("Lan+"), "2" = c("Lan(+wk)"), "3" = c("Lan-"), TGP = c("-")
    ), {
        hier <- process_hier(df_allele_table$Gene == "ABCB6")
    })

VEL <- within(list(
    ref = c("VEL:1 (Vel+)/VEL*01"), allele_name = c("/VEL.*", "/\\*1000G.*"),
    "1" = c("VEL:1 (Vel+)"), TGP = c("-")
    ), {
        hier <- process_hier(df_allele_table$Gene == "SMIM1")
    })

CD59 <- within(list(
    ref = c("CD59:+1 or CD59.1+/CD59*01"), allele_name = c("/CD59.*", "/\\*1000G.*"),
    "1" = c("CD59.1+"), TGP = c("-")
    ), {
        hier <- process_hier(df_allele_table$Gene == "CD59")
    })

AUG <- within(list(
    ref = c("AUG1+, At(a+), ATML−, ATAM+ AUG:1,2,−3,4/AUG*01"), allele_name = c("/AUG.*", "/\\*1000G.*"),
    "1" = c("AUG1+, At(a+), ATML−, ATAM+ AUG:1,2,−3,4"), TGP = c("-")
    ), {
        hier <- process_hier(df_allele_table$Gene == "SLC29A1")
    })

GATA1 <- within(list(
    ref = c("Common/GATA1*01"), allele_name = c("/GATA1.*", "/\\*1000G.*"),
    "1" = c("Common"), TGP = c("-")
    ), {
        hier <- process_hier(df_allele_table$Gene == "GATA1")
    })

KLF1 <- within(list(
    ref = c("Common/KLF1*01"), allele_name = c("/KLF1.*", "/\\*1000G.*"),
    "1" = c("Common"), "2" = c("In(Lu)"), TGP = c("-")
    ), {
        hier <- process_hier(df_allele_table$Gene == "KLF1")
    })


my_dict <- list(
    "ABO" = ABO, "P1PK" = P1PK, "LU" = LU, "KEL" = KEL, "LE" = LE, "FY" = FY, "JK" = JK, "DI" = DI, "YT" = YT, "XG" = XG,
    "SC" = SC, "DO" = DO, "CO" = CO, "LW" = LW, "CHRG_A" = CHRG_A, "CHRG_B" = CHRG_B, "H_FUT1" = H_FUT1, "H_FUT2" = H_FUT2, "KX" = KX, "GE" = GE, 
    "CROM" = CROM, "KN" = KN, "IN" = IN, "OK" = OK, "RAPH" = RAPH, "JMH" = JMH, "I" = I, "GLOB" = GLOB, "GIL" = GIL, "RHAG" = RHAG, "FORS" = FORS, 
    "LAN" = LAN, "VEL" = VEL, "CD59" = CD59, "AUG" = AUG, "GATA1" = GATA1, "KLF1" = KLF1, "JR" = JR, "MNS_GYPA" = MNS_GYPA, "MNS_GYPB" = MNS_GYPB,
    "RHCE" = RHCE, "RHD" = RHD
)

for (file in pred_files) {
    
    file_name <- sub("(.*)_summary.*", "\\1", file)

    if (file.exists(paste0(file_name, "_final_summary_table.Rdata"))) {
        next
    }

    print(paste0("Processing", file))
    load(file)
    col <- colnames(df_sum_table)
    
    if (col == "RH") {
        col <- sub(".*/(.*)\\*.*", "\\1", df_sum_table[[1]][1])
        colnames(df_sum_table)[1] <- col
    }

    else if (col %in% c("H", "MNS")) {
        col <- paste0(col, "_", sub(".*/(.*)\\*.*", "\\1", df_sum_table[[1]][1]))
        colnames(df_sum_table)[1] <- col
    }

    else if (col == "CHRG") {
        col <- paste0(col, "_", sub(".*C4([A-Za-z])\\*.*", "\\1", df_sum_table[[1]][1]))
        colnames(df_sum_table)[1] <- col
    }

    df_sum_table <- df_sum_table %>%
        mutate({{ col }} := process_group(.data[[col]], my_dict[[col]]))

    
    save(df_sum_table, file = paste0(file_name, "_final_summary_table.Rdata"))
    write.table(df_sum_table, file = paste0(file_name, "_final_summary_table.tsv"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

}

print("Done")