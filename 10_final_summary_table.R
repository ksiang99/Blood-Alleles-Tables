# R script to generate summary table of inferred phenotype for each sample in each blood group
rm(list = ls())

setwd("/home/svu/e0726996/Blood-type-GWAS/")

# Set random seed
set.seed(0)

library(dplyr)
library(stringr)

# Set paths
PATH_SUM_TABLE <- "./R/results/summary_table.tsv"
PATH_ALLELE_TABLE <- "./R/Blood Allele Table.tsv"
df_sum_table <- read.delim(PATH_SUM_TABLE, header = TRUE)
df_allele_table <- read.delim(PATH_ALLELE_TABLE, header = TRUE)
df_allele_table$Pheno_Allele <- paste0(df_allele_table$Phenotype, "/", df_allele_table$Allele.Name)

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
# others: Provide a general phenotype name to ISBT phenotype name
# TGP: To represent 1000G alleles
# Currently, the general phenotype names are long/weird

ABO <- within(list(
    ref = c("A1/ABO*A1.01"),
    allele_name = c("/ABO.*", "/\\*1000G.*"),
    cisAB = c("cisAB"),
    B_A = c("B(A)"),
    A = c("A1", "A2", "A3"),
    A_weak = c("Aweak", "Ax/Aweak", "Ael"),
    B = c("B"),
    B_weak = c("B3", "Bel"),
    O = c("O"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "ABO")
    }
)

MNS_GYPA <- within(list(
    ref = c("MNS:1 or M+/GYPA*01"),
    allele_name = c("/GYPA.*", "/\\*1000G.*"),
    M = c("M+"),
    N = c("N+"),
    Mc = c("Mc+"),
    TGP_GYPA = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "GYPA")
    }
)

MNS_GYPB <- within(list(
    ref = c("MNS:4 or s+/GYPB*04"),
    allele_name = c("/GYPB.*", "/\\*1000G.*"),
    s_small = c("s+"),
    S = c("S+"),
    sD = c("sD+"),
    Mit = c("Mit+"),
    U_var = c("S–U+var"),
    Mv = c("Mv+"),
    s_null = c("S-s-U-"),
    S_weak = c("S+w"),
    s_weak = c("s+w, U+w"),
    S_altered = c("MNS:4,5 (s+, U+) altered GPB"),
    TGP_GYPB = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "GYPB")
    }
)

P1PK <- within(list(
    ref = c("P1+ Pk+/A4GALT*01"),
    allele_name = c("/A4GALT.*", "/\\*1000G.*"),
    P1 = c("P1+ Pk+"),
    P2 = c("P1– Pk+ (P2)"),
    p = c("p"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "A4GALT")
    }
)

RHCE <- list()

RHD <- list()

LU <- within(list(
    ref = c("LU:2 or Lu(b+)/LU*02"),
    allele_name = c("/LU.*", "/\\*1000G.*"),
    Lu_a = c("Lu(a+)"),
    Lu_sixteen = c("LU:-16"),
    Lu_nineteen = c("LU:1,19"),
    Au_b = c("Au(a−b+)"),
    Lu_nine = c("LU:-6,9"),
    Lu_fourteen = c("LU:-8,14"),
    LURC_null = c("LU:-22, LURC−"),
    Lu_b = c("Lu(b+)"),
    Lu_null = c("Lunull"),
    LUBI_null = c("LU:-26, LUBI−"),
    LURA_null = c("LU:-29, LURA−"),
    Lu_thirteen = c("LU:-13"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "BCAM")
    }
)

KEL <- within(list(
    ref = c("KEL:2 or k+/KEL*02"),
    allele_name = c("/KEL.*", "/\\*1000G.*"),
    K_big = c("K+k-"),
    K_small = c("k+"),
    Kel_null = c("K0"), 
    Kp_a = c("Kp(a+b–c–)"),
    Js_a = c("Js(a+b–)"), 
    Ul_a = c("Ul(a+)"),
    KHOZ = c("KHIZ-, KHOZ+"),
    K_eighteen_null = c("KEL:-18"), 
    K_nineteen_null = c("KEL:-19"), 
    TOU_null = c("TOU-"), 
    KTIM_null = c("KTIM-"), 
    KYO = c("KYO+, KYOR-"), 
    KEAL = c("KHUL-,KEAL+"), 
    Kel_weak = c("KEL:1weak", "Kmod", "K0 phenotype"), # Affects both K and k
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "KEL")
    }
)

LE <- within(list(
    ref = c("FUT3 Active/FUT3*01.01"),
    allele_name = c("/FUT3.*", "/\\*1000G.*"),
    FUT3_active = c("FUT3 Active"),
    FUT3_weak = c("FUT3 Active (Weak)"),
    Le_null = c("FUT3 Inactive - Le(a-b-)"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "FUT3")
    }
)

FY <- within(list(
    ref = c("FY:1 or Fy(a+)/FY*01"),
    allele_name = c("/FY.*", "/\\*1000G.*"),
    Fy_a = c("Fy(a+)"),
    Fy_a_weak = c("Fy(a+w)"),
    Fy_b = c("Fy(b+)"),
    Fy_b_weak = c("Fy(b+w), Fyx"),
    Fy_null = c("Fy(a−b−)"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "ACKR1")
    }
)

JK <- within(list(
    ref = c("JK:1 or Jk(a+)/JK*01"),
    allele_name = c("/JK.*", "/\\*1000G.*"),
    Jk_a = c("Jk(a+)"),
    Jk_a_weak = c("Jk(a+W)", "Jk(a+w)"),
    Jk_b = c("Jk(b+)"),
    Jk_b_weak = c("Jk(b+W)", "Jk(b+w)"),
    Jk_null = c("Jk(a–b–)"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "SLC14A1")
    }
)   


DI <- within(list(
    ref = c("DI:–1,2 or Di(a–b+)/DI*02"),
    allele_name = c("/DI.*", "/\\*1000G.*"),
    Di_a = c("Di(a+b–)"),
    Di_b = c("Di(a–b+)"), # "DI:23",
    Wr_a = c("Wr(a+b–)"),
    Wd_a = c("Wd(a+)"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "SLC4A1")
    }
)

YT <- within(list(
    ref = c("YT:1,-2 or Yt(a+b-)/YT* 01"),
    allele_name = c("/YT.*", "/\\*1000G.*"),
    Yt_a = c("Yt(a+b-)"),
    Yt_b = c("Yt(a-b+)"),
    YTEG_null = c("YTEG-"),
    YTLI_null = c("YTLI-"),
    YTOT_null = c("YTOT-"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "ACHE")
    }
)

XG <- within(list(
    ref = c("Xga/XG*01"),
    allele_name = c("/XG.*", "/\\*1000G.*"),
    Xg_a = c("Xga"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "XG")
    }
)

SC <- within(list(
    ref = c("SC:1 or Sc1+/SC*01"),
    allele_name = c("/SC.*", "/\\*1000G.*"),
    Sc_1 = c("Sc1+"),
    Sc_2 = c("Sc2+"),
    STAR_null = c("STAR–"),
    SCAC_null = c("SCAC–"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "ERMAP")
    }
)

DO <- within(list(
    ref = c("DO:2 or Do(b+)/DO*02"),
    allele_name = c("/DO.*", "/\\*1000G.*"),
    Do_a = c("Do(a+)"),
    Do_b = c("Do(b+)"),
    Jo_a_null = c("Jo(a–)"),
    DODE_null = c("DODE–"),
    DOLG_null = c("DOLG–"),
    Hy_null = c("Hy–"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "ART4")
    }
)

CO <- within(list(
    ref = c("CO:1 or Co(a+)/CO*01.01"),
    allele_name = c("/CO.*", "/\\*1000G.*"),
    Co_a = c("Co(a+)"),
    Co_b = c("Co(b+)"),
    Co_null = c("Co(a–b–)"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "AQP1")
    }
)

LW <- within(list(
    ref = c("LW:5 or LW(a+)/LW*05"),
    allele_name = c("/LW.*", "/\\*1000G.*"),
    Lw_a = c("LW(a+)"),
    Lw_b = c("LW(b+)"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "ICAM4")
    }
)

CHRG_A <- within(list(
    ref = c("Ch–Rg+ or CH:–1,–2,–3,–4,–5,–6 RG:1,2/C4A*03"),
    allele_name = c("/C4A.*", "/\\*1000G.*"),
    Chrg_a = c("CH:–1,–2,–3,–4,–5,–6 RG:1,2"),
    TGP_a = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "C4A")
    }
)

CHRG_B <- within(list(
    ref = c("Ch+Rg– or CH:1,2,3,4,5,6 RG:–1,–2/C4B*03"),
    allele_name = c("/C4B.*", "/\\*1000G.*"),
    Chrg_b = c("CH:1,2,3,4,5,6 RG:–1,–2"),
    TGP_b = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "C4B")
    }
)

H1 <- within(list(
    ref = c("H+/FUT1*01"),
    allele_name = c("/FUT1.*", "/\\*1000G.*"),
    H1 = c("H+"),
    H1_weak = c("H+weak"),
    H1_null = c("H–"),
    H1_TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "FUT1")
    }
)

H2 <- within(list(
    ref = c("H+/FUT2*01"),
    allele_name = c("/FUT2.*", "/\\*1000G.*"),
    H2 = c("H+"),
    H2_weak = c("H+w"),
    H2_null = c("H–"),
    H2_TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "FUT2")
    }
)

KX <- within(list(
    ref = c("XK:1 or Kx+/XK*01"),
    allele_name = c("/XK.*", "/\\*1000G.*"),
    Xk = c("Kx+"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "XK")
    }
)

GE <- within(list(
    ref = c("GE:2,3,4/GE*01"),
    allele_name = c("/GE.*", "/\\*1000G.*"),
    Ge = c("GE:2,3,4"),
    GEPL_null = c("GEPL–"),
    GECT_null = c("GECT–"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "GYPC")
    }
)

CROM <- within(list(
    ref = c("CROM:1 or Cra+)/CROM*01"),
    allele_name = c("/CROM.*", "/\\*1000G.*"),
    Crom_a = c("Cra+)"),
    Crom_a_null = c("Cr(a–)"),
    Tc_b = c("Tc(b+)"),
    Wes_a = c("WES(a+)"),
    SERF_null = c("SERF–"),
    UMC_null = c("UMC–"),
    CRUE_null = c("CRUE−"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "CD55")
    }
)

KN <- within(list(
    ref = c("KN:1 or Kn(a+) KN:3 or McC(a+) KN:4 or Sl1+ KN:5 or Yk(a+) KN:8 or Sl3+ KN:9 or KCAM+ KN:11 or DACY+/KN*01"),
    allele_name = c("/KN.*", "/\\*1000G.*"),
    Kn_a = c("KN:1 or Kn(a+) KN:3 or McC(a+) KN:4 or Sl1+ KN:5 or Yk(a+) KN:8 or Sl3+ KN:9 or KCAM+ KN:11 or DACY+"),
    Kn_two_ten = c("KN:2 or Kn(a-b+) KN:-9 or KCAM- KN:10 or KDAS+"),
    Yk_a_null = c("KN:–5 or Yk(a–)"),
    Kn_one_six = c("KN:-4 or Sl1- KN:-3,6 or McC(a-b+) KN:7 or Vil+ KN:-9 or KCAM- KN:10 or KDAS+"),
    Kn_one_seven = c("KN:7 or Vil+ KN:-9 or KCAM- KN:10 or KDAS+"),
    Kn_one_ten = c("KN:-9 or KCAM- KN:10 or KDAS+"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "CR1")
    }
)

IN <- within(list(
    ref = c("In(a–b+)/IN*02"),
    allele_name = c("/IN.*", "/\\*1000G.*"),
    In_a = c("In(a+b–)"),
    In_b = c("In(a–b+)"),
    INFI_null = c("INFI–"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "CD44")
    }
)

OK <- within(list(
    ref = c("OK:1 or Ok(a+)/OK*01.01"),
    allele_name = c("/OK.*", "/\\*1000G.*"),
    Ok_a = c("Ok(a+)"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "BSG")
    }
)

RAPH <- within(list(
    ref = c("RAPH:1 or MER2+/RAPH*01"),
    allele_name = c("/RAPH.*", "/\\*1000G.*"),
    RAPH_1 = c("MER2+"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "CD151")
    }
)

JMH <- within(list(
    ref = c("JMH:1 or JMH+/JMH*01"),
    allele_name = c("/JMH.*", "/\\*1000G.*"),
    Jmh = c("JMH+"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "SEMA7A")
    }
)

I <- within(list(
    ref = c("I:1 or I+/GCNT2*01"),
    allele_name = c("/GCNT2.*", "/\\*1000G.*"),
    I = c("I+"),
    I_weak = c("I+W"),
    I_null = c("I–"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "GCNT2")
    }
)

GLOB <- within(list(
    ref = c("GLOB:1 (P+)/GLOB*01"),
    allele_name = c("/GLOB.*", "/\\*1000G.*"),
    P = c("GLOB:1 (P+)"),
    P_null = c("GLOB:–1 (P–)"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "B3GALNT1")
    }
)

GIL <- within(list(
    ref = c("GIL:1 or GIL+/GIL*01", "/\\*1000G.*"),
    allele_name = c("/GIL.*", "/\\*1000G.*"),
    Gil = c("GIL+"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "AQP3")
    }
)

RHAG <- within(list(
    ref = c("RHAG:1 or Duclos+/RHAG*01"),
    allele_name = c("/RHAG.*", "/\\*1000G.*"),
    Rhag = c("Duclos+"),
    DSLK_null = c("DSLK−, Kg+"),
    Rhmod = c("Rhmod"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "RHAG")
    }
)

FORS <- within(list(
    ref = c("FORS:-1 (FORS-)/GBGT1*01N.01"),
    allele_name = c("/GBGT1.*", "/\\*1000G.*"),
    Fors_null = c("FORS:-1 (FORS-)"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "GBGT1")
    }
)

JR <- within(list(
    ref = c("Jr(a+)/ABCG2*01"),
    allele_name = c("/ABCG2.*", "/N/A.*", "/\\*1000G.*"),
    Jr_a = c("Jr(a+)"),
    Jr_weak = c("Jr(a+w)"),
    Jr_null = c("Jr(a−)"),
    Jr_unclear = c("Unclear Jra phnotype"),
    TGP = c("-")
    ), {
    hier <- process_hier(df_allele_table$Gene == "ABCG2")
    }

)

LAN <- within(list(
    ref = c("Lan+/ABCB6*01"),
    allele_name = c("/ABCB6.*", "/\\*1000G.*"),
    Lan = c("Lan+"),
    Lan_weak = c("Lan(+wk)"),
    Lan_null = c("Lan-"),
    TGP = c("-")
    ), {
        hier <- process_hier(df_allele_table$Gene == "ABCB6")
    }
)

VEL <- within(list(
    ref = c("VEL:1 (Vel+)/VEL*01"),
    allele_name = c("/VEL.*", "/\\*1000G.*"),
    hier = c("VEL:1 (Vel+)/VEL*01"),
    Vel = c("VEL:1 (Vel+)"),
    TGP = c("-")
    ), {
        hier <- process_hier(df_allele_table$Gene == "SMIM1")
    }
)

CD59 <- within(list(
    ref = c("CD59:+1 or CD59.1+/CD59*01"),
    allele_name = c("/CD59.*", "/\\*1000G.*"),
    Cd59 = c("CD59.1+"),
    TGP = c("-")
    ), {
        hier <- process_hier(df_allele_table$Gene == "CD59")
    }

)

AUG <- within(list(
    ref = c("AUG1+, At(a+), ATML−, ATAM+ AUG:1,2,−3,4/AUG*01"),
    allele_name = c("/AUG.*", "/\\*1000G.*"),
    Aug = c("AUG1+, At(a+), ATML−, ATAM+ AUG:1,2,−3,4"),
    TGP = c("-")
    ), {
        hier <- process_hier(df_allele_table$Gene == "SLC29A1")
    }
)

GATA1 <- within(list(
    ref = c("Common/GATA1*01"),
    allele_name = c("/GATA1.*", "/\\*1000G.*"),
    Common = c("Common"),
    TGP = c("-")
    ), {
        hier <- process_hier(df_allele_table$Gene == "GATA1")
    }
)

KLF1 <- within(list(
    ref = c("Common/KLF1*01"),
    allele_name = c("/KLF1.*", "/\\*1000G.*"),
    Common = c("Common"),
    In_lu = c("In(Lu)"),
    TGP = c("-")
    ), {
        hier <- process_hier(df_allele_table$Gene == "KLF1")
})


my_dict <- list(
    "ABO" = ABO, "P1PK" = P1PK, "LU" = LU, "KEL" = KEL, "LE" = LE, "FY" = FY, "JK" = JK, "DI" = DI, "YT" = YT, "XG" = XG,
    "SC" = SC, "DO" = DO, "CO" = CO, "LW" = LW, "CHRG_A" = CHRG_A, "CHRG_B" = CHRG_B, "H1" = H1, "H2" = H2, "KX" = KX, "GE" = GE, 
    "CROM" = CROM, "KN" = KN, "IN" = IN, "OK" = OK, "RAPH" = RAPH, "JMH" = JMH, "I" = I, "GLOB" = GLOB, "GIL" = GIL, "RHAG" = RHAG, "FORS" = FORS, 
    "LAN" = LAN, "VEL" = VEL, "CD59" = CD59, "AUG" = AUG, "GATA1" = GATA1, "KLF1" = KLF1, "JR" = JR, "MNS_GYPA" = MNS_GYPA, "MNS_GYPB" = MNS_GYPB
)

# Iterate through each blood system to infer the correct phenotype for each sample
df_sum_table <- df_sum_table %>%
    mutate(across(
        c("ABO", "P1PK", "LU", "KEL", "LE", "FY", "JK", "DI", "YT", "SC", "DO", "CO", "LW", 
          "KX", "GE", "CROM", "IN", "OK", "RAPH", "JMH", "I", "GLOB", "GIL", "RHAG", "FORS", 
          "LAN", "VEL", "CD59", "AUG", "GATA1", "KLF1", "JR", "KN", "XG"), 
        ~ process_group(.x, my_dict[[cur_column()]])
    )) %>%
    mutate(H = paste0(
        sub(" \\| .*", "", process_group(H, my_dict[["H1"]])), 
        " | ", 
        process_group(H, my_dict[["H2"]])
    )) %>%
    mutate(CHRG = paste0(
        sub(" \\| .*", "", process_group(CHRG, my_dict[["CHRG_A"]])), 
        " | ", 
        process_group(CHRG, my_dict[["CHRG_B"]])
    )) %>%
    mutate(MNS = paste0(
        sub(" \\| .*", "", process_group(MNS, my_dict[["MNS_GYPA"]])), 
        " | ", 
        process_group(MNS, my_dict[["MNS_GYPB"]])
    ))

save(df_sum_table, file = "./R/results/final_summary_table.Rdata")

print("Done")

# df_sum_table <- df_sum_table %>%
#     mutate(
#         LE = case_when(
#             grepl("Le_null", LE) ~ "Le_null: Le(a-b-)",
#             grepl("FUT3_active", LE) & grepl("H2_null", H) ~ "Le_a: Le(a+b-)",
#             grepl("FUT3_active", LE) & grepl("H2_weak", H) ~ "Le_a_b: Le(a+b+)",
#             TRUE ~ "Le_b: Le(a-b+)"
#         ),
        
#         ABO = ifelse(grepl("H1_null", H), "O: O (Due to H1 null)", ABO),
        
#         GLOB = ifelse(grepl("p: p", P1PK), "P_null: P– (Due to p)", GLOB),
        
#         LU = ifelse(grepl("In_lu", KLF1), "Lu_null: Lunull (Due to KLF1 Mutation)", LU),
        
#         IN = ifelse(grepl("In_lu", KLF1), "In_null: In(a–b–) (Due to KLF1 Mutation)", IN)
#     )