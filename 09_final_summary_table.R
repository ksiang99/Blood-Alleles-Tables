rm(list = ls())

setwd("/home/svu/e0726996/Blood-type-GWAS/")

# Set random seed
set.seed(0)

library(dplyr)

# Create "dictionary" for each blood group
# ref: reference allele of blood group
# allele_name: pattern used to clean up string at the end
# hier: Alleles found in df_sum_table, ranked based on # nucleotide changes associated to allele
# Null alleles largest # NCs -> smallest, Remaining alleles largest -> smallest)
# Need to add in remaining alleles!!!!!
# Others: Provide a generalised name for phenotypes For example: "AB: cis(AB)"

ABO <- list(
    ref = c("A1/ABO*A1.01"),
    allele_name = c("/ABO.*"),
    hier = c("O/ABO*O.02.03", "O/ABO*O.02.01", "O/ABO*O.06", "B/ABO*B.02", "B3/ABO*B3.02", 
    "B/ABO*B.01", "A2/ABO*A2.04", "Aweak/ABO*AW.09", "cisAB/ABO*cisAB.05", "cisAB/ABO*cisAB.06",
    "B(A)/ABO*BA.03", "Ax/Aweak/ABO*AW.31.01", "B(A)/ABO*BA.01", "Ax/Aweak/ABO*AW31.02- 05", "Ael/ABO*AEL.07", 
    "cisAB/ABO*cisAB.02", "A2/ABO*A2.09", "Aweak/ABO*AW.25", "Ael/ABO*AEL.02", "A2/ABO*A2.01",
    "A2/ABO*A2.20", "A3/ABO*A3.02", "Aweak/ABO*AW.27", "Ax/Aweak/ABO*AW.30.02", "A1/ABO*A1.02", 
    "A2/ABO*A2.02", "A2/ABO*A2.03", "A2/ABO*A2.06", "Ax/Aweak/ABO*AW.30.01", "A1/ABO*A1.01"),
    AB = c("cisAB", "B(A)"),
    A = c("A1", "A2", "A3", "Aweak", "Ax/Aweak", "Ael"),
    B = c("B", "B3", "Bel"),
    O = c("O")
)

P1PK <- list(
    ref = c("P1+ Pk+/A4GALT*01"),
    allele_name = c("/A4GALT.*"),
    hier = c("p/A4GALT*01N.16", "P1+ Pk+/A4GALT*01.02", "P1+ Pk+/A4GALT*01"),
    P1 = c("P1+ Pk+", "P1+ Pk+ NOR+"),
    P2 = c("P1– Pk+ (P2)"),
    p = c("p")
)

RHCE <- list()

RHD <- list()

LU <- list(
    ref = c("LU:2 or Lu(b+)/LU*02"),
    allele_name = c("/LU.*"),
    hier = c("Lunull/LU*02N.03", "LU:-13/LU*02.-13", "LU:-16/LU*01.-16", "LU:1,19/LU*01.19", "Lu(a+)/LU*01",  
    "Au(a−b+)/LU*02.19", "LU:-5/LU*02.-05", "LU:-6,9/LU*02.09", "LU:-8,14/LU*02.14", "LU:-22, LURC−/LU*02.–22", 
    "LU:-26, LUBI−/LU*02.-26", "LU:-29, LURA−/LU*02.-29", "Lu(b+)/LU*02"),
    Lu_a = c("Lu(a+)", "LU:-16", "LU:1,19"),
    Lu_b = c("Lu(b+)", "LU:-5", "LU:-6,9", "LU:-8,14", "Au(a−b+)", "LU:-22, LURC−", "LU:-26, LUBI−", "LU:-29, LURA−"),
    Lu_null = c("Lunull")
)

# 02N.17 has the same NC position as 02.02 but different base change.
# Base change for 02.02 occurs more frequent. Hence, 02N.17 is listed at the back
KEL <- list(
    ref = c("KEL:2 or k+/KEL*02"),
    allele_name = c("/KEL.*"),
    hier = c("K0/KEL*02N.11", "K0/KEL*02N.12", "K0/KEL*02N.13", "K0/KEL*02N.16", "K0/KEL*02N.39", 
    "K0/KEL*02N.42", "K0/KEL*02N.44", "Kmod/KEL*02M.04", "K0 phenotype/KEL*02M.05", "K+k-/KEL*01.01", 
    "k+/KEL*02.02", "Kp(a+b–c–)/KEL*02.03", "Js(a+b–)/KEL*02.06", "Ul(a+)/KEL*02.10", "KEL: -11,17/KEL*02.17", 
    "KEL:-18/KEL*02.-18.1", "KEL:-18/KEL*02.-18.2", "KEL:-19/KEL*02.-19", "TOU-/KEL*02.-26", "KTIM-/KEL*02.-30",
    "KYO+, KYOR-/ KEL*02.31", "KETI-/KEL*02.-36", "KHUL-,KEAL+/KEL*02.39", "KHIZ-, KHOZ+/KEL*02.41", "K0/KEL*02N.17",
    "KEL:1weak/KEL*01M.01", "k+/KEL*02"),
    K = c("K+k-"),
    k = c("k+", "Kp(a+b–c–)", "Js(a+b–)", "Ul(a+)", "KEL: -11,17", "KEL:-18", "KEL:-19", "TOU-", "KTIM-", "KYO+, KYOR-", "KETI-", 
    "KHUL-,KEAL+", "KHIZ-, KHOZ+"),
    Kel_null = c("K0"),
    Kel_weak = c("KEL:1weak", "Kmod", "K0 phenotype") # Affects both K and k
)

LE <- list(
    ref = c("FUT3 Active/FUT3*01.01"),
    allele_name = c("/FUT3.*"),
    hier = c("FUT3 Inactive - Le(a-b-)/FUT3*01N.03.12", "FUT3 Inactive - Le(a-b-)/FUT3*01N.01.08", "FUT3 Inactive - Le(a-b-)/FUT3*01N.01.09", "FUT3 Inactive - Le(a-b-)/FUT3*01N.01.10", "FUT3 Inactive - Le(a-b-)/FUT3*01N.01.11",
    "FUT3 Inactive - Le(a-b-)/FUT3*01N.17.16", "FUT3 Inactive - Le(a-b-)/FUT3*01N.01.02", "FUT3 Inactive - Le(a-b-)/FUT3*01N.01.05", "FUT3 Inactive - Le(a-b-)/FUT3*01N.01.12", "FUT3 Inactive - Le(a-b-)/FUT3*01N.03.05", 
    "FUT3 Inactive - Le(a-b-)/FUT3*01N.03.08", "FUT3 Inactive - Le(a-b-)/FUT3*01N.03.11", "FUT3 Inactive - Le(a-b-)/FUT3*01N.17.07", "FUT3 Inactive - Le(a-b-)/FUT3*01N.17.12", "FUT3 Inactive - Le(a-b-)/FUT3*01N.17.15",
    "FUT3 Inactive - Le(a-b-)/FUT3*01N.17.17", "FUT3 Inactive - Le(a-b-)/FUT3*01N.17.19", "FUT3 Inactive - Le(a-b-)/FUT3*01N.01.01", "FUT3 Inactive - Le(a-b-)/FUT3*01N.03.02", "FUT3 Inactive - Le(a-b-)/FUT3*01N.03.04",
    "FUT3 Inactive - Le(a-b-)/FUT3*01N.17.02", "FUT3 Inactive - Le(a-b-)/FUT3*01N.17.03", "FUT3 Inactive - Le(a-b-)/FUT3*01N.17.04", "FUT3 Inactive - Le(a-b-)/FUT3*01N.17.05", "FUT3 Inactive - Le(a-b-)/FUT3*01N.03.01",
    "FUT3 Active/FUT3*01.16", "FUT3 Active/FUT3*01.15", "FUT3 Active/FUT3*01.03", "FUT3 Active/FUT3*01.04", "FUT3 Active/FUT3*01.05",
    "FUT3 Active/FUT3*01.06", "FUT3 Active/FUT3*01.10", "FUT3 Active/FUT3*01.11", "FUT3 Active/FUT3*01.13", "FUT3 Active/FUT3*01.17.01",
    "FUT3 Active/FUT3*01.01"),
    FUT3_active = c("FUT3 Active"),
    Le_null = c("FUT3 Inactive - Le(a-b-)")
)

FY <- list(
    ref = c("FY:1 or Fy(a+)/FY*01"),
    allele_name = c("/FY.*"),
    hier = c("Fy(a−b−)/FY*02N.01", "Fy(a−b−)/FY*01N.01", "Fy(a−b−)/FY*01N.06", "Fy(b+w), Fyx/FY*02W.01", "Fy(a+w)/FY*01W.02", 
    "Fy(b+)/FY*02", "Fy(a+w)/FY*01W.01", "Fy(a+)/FY*01"),
    Fy_a = c("Fy(a+)", "Fy(a+w)"),
    Fy_b = c("Fy(b+)", "Fy(b+w), Fyx"),
    Fy_null = c("Fy(a−b−)")
)

JK <- list(
    ref = c("JK:1 or Jk(a+)/JK*01"),
    allele_name = c("/JK.*"),
    hier = c("Jk(a–b–)/JK*01N.20", "Jk(a–b–)/JK*02N.13", "Jk(a–b–)/JK*01N.17", "Jk(a–b–)/JK*02N.01", "Jk(a–b–)/JK*02N.02",
    "Jk(a–b–)/JK*02N.06", "Jk(a–b–)/JK*02N.07", "Jk(a–b–)/JK*02N.08", "Jk(a–b–)/JK*02N.09", "Jk(a–b–)/JK*02N.17",
    "Jk(a–b–)/JK*02N.21", "Jk(a–b–)/JK*01N.03", "Jk(a–b–)/JK*01N.04", "Jk(a–b–)/JK*01N.06", "Jk(a–b–)/JK*01N.18",
    "Jk(a–b–)/JK*01N.19", "Jk(b+w)/JK*02W.03", "Jk(b+W)/JK*02W.04", "Jk(a+w)/JK*01W.06", "Jk(a+W)/JK*01W.11",
    "Jk(b+)/JK*02", "Jk(a+W)/JK*01W.01", "Jk(a+W)/JK*01W.02", "Jk(a+W)/JK*01W.03", "Jk(a+W)/JK*01W.04",
    "Jk(a+W)/JK*01W.09", "Jk(a+)/JK*01"),
    Jk_a = c("Jk(a+)", "Jk(a+W)", "Jk(a+w)"),
    Jk_b = c("Jk(b+)", "Jk(b+W)", "Jk(b+w)"),
    Jk_null = c("Jk(a–b–)")
)   


DI <- list(
    ref = c("DI:–1,2 or Di(a–b+)/DI*02"),
    allele_name = c("/DI.*"),
    hier = c("Di(a+b–)/DI*01", "Wr(a+b–)/DI*02.03", "Wr(a-b+)/DI*02.04", "Wd(a+)/DI*02.05", "DI:23/DI*02.23",
    "Di(a–b+)/DI*02"),
    Di_a = c("Di(a+b–)"),
    Di_b = c("Wr(a+b–)", "Wr(a-b+)", "Wd(a+)", "DI:23", "Di(a–b+)")
)

YT <- list(
    ref = c("YT:1,-2 or Yt(a+b-)/YT* 01"),
    allele_name = c("/YT.*"),
    hier = c("Yt(a-b+)/YT*02", "YTEG-/YT*01.-03", "YTLI-/YT*01.-04", "YTOT-/YT*01.-05", "Yt(a+b-)/YT* 01"),
    Yt_a = c("YTEG-", "YTLI-", "YTOT-", "Yt(a+b-)"),
    Yt_b = c("Yt(a-b+)")
)

XG <- list(
    ref = c("Xga/XG*01"),
    allele_name = c("/XG.*"),
    hier = c("Xga/XG*01"),
    Xg_a = c("Xga")
)

SC <- list(
    ref = c("SC:1 or Sc1+/SC*01"),
    allele_name = c("/SC.*"),
    hier = c("Sc2+/SC*02", "SCAC–/SC*01.–09", "STAR–/SC*01.–05", "SCAN–/SC*01.–07", "Sc1+/SC*01"),
    Sc_1 = c("SCAC–", "STAR–", "SCAN–", "Sc1+"),
    Sc_2 = c("Sc2+")
)

# Too many ors
DO <- list(
    ref = c("DO:2 or Do(b+)/DO*02 or DO*B"),
    allele_name = c("/DO.*"),
    hier = c("Jo(a–)/DO*01.–05", "DOLG–/DO*01.–08", "DODE–/DO*01.–10", "Do(a+)/DO*01", "Do(b+)/DO*02"),
    Do_a = c("Jo(a–)", "DOLG–", "DODE–", "Do(a+)"),
    Do_b = c("Do(b+)")
)

CO <- list(
    ref = c("CO:1 or Co(a+)/CO*01.01"),
    allele_name = c("/CO.*"),
    hier = c("Co(a–b–)/CO*01N.06", "Co(b+)/CO*02.01", "Co(b+)/CO*02.02", "Co(a+)/CO*01.01"),
    Co_a = c("Co(a+)"),
    Co_b = c("Co(b+)"),
    Co_null = c("Co(a–b–)")
)

LW <- list(
    ref = c("LW:5 or LW(a+)/LW*05"),
    allele_name = c("/LW.*"),
    hier = c("LW(b+)/LW*07", "LW(a+)/LW*05"),
    Lw_a = c("LW(a+)"),
    Lw_b = c("LW(b+)")
)

CHRG_A <- list(
    ref = c("Ch–Rg+ or CH:–1,–2,–3,–4,–5,–6 RG:1,2/C4A*03"),
    allele_name = c("/C4A.*"),
    hier = c("CH:–1,–2,–3,–4,–5,–6 RG:1,2/C4A*03"),
    Chrg_a = c("CH:–1,–2,–3,–4,–5,–6 RG:1,2")
)

CHRG_B <- list(
    ref = c("Ch+Rg– or CH:1,2,3,4,5,6 RG:–1,–2/C4B*03"),
    allele_name = c("/C4B.*"),
    hier = c("CH:1,2,3,4,5,6 RG:–1,–2/C4B*03"),
    Chrg_b = c("CH:1,2,3,4,5,6 RG:–1,–2")
)

H1 <- list(
    ref = c("H+/FUT1*01"),
    allele_name = c("/FUT1.*"),
    hier = c("H–/FUT1*01N.09", "H–/FUT1*01N.21", "H+weak/FUT1*01W.32", "H+weak/FUT1*01W.02", "H+weak/FUT1*01W.09",
    "H+weak/FUT1*01W.24", "H+/FUT1*01.02", "H+/FUT1*01"),
    H1 = c("H+"),
    H1_weak = c("H+weak"),
    H1_null = c("H–")
)

H2 <- list(
    ref = c("H+/FUT2*01"),
    allele_name = c("/FUT2.*"),
    hier = c("H–/FUT2*01N.02", "H–/FUT2*01N.03", "H–/FUT2*01N.04", "H–/FUT2*01N.05", "H–/FUT2*01N.06",
    "H–/FUT2*01N.12", "H–/FUT2*01N.13", "H–/FUT2*01N.14", "H–/FUT2*01N.15", "H–/FUT2*01N.16",
    "H–/FUT2*01N.17","H+/FUT2*01.03.02", "H+/FUT2*01.02", "H+/FUT2*01.03.01", "H+/FUT2*01.04",
    "H+/FUT2*01.06", "H+/FUT2*01.08", "H+w/FUT2*01W.01", "H+w/FUT2*01W.02.01",  "H+/FUT2*01"),
    H2 = c("H+", "H+w"),
    H2_weak = c("H+w")
    H2_null = c("H–")
)

KX <- list(
    ref = c("XK:1 or Kx+/XK*01"),
    allele_name = c("/XK.*"),
    hier = c("Kx+/XK*01"),
    Xk = c("Kx+")
)

GE <- list(
    ref = c("GE:2,3,4/GE*01"),
    allele_name = c("/GE.*"),
    hier = c("GEPL–/GE*01.–10", "GECT–/GE*01.–13", "GE:2,3,4/GE*01"),
    Ge = c("GEPL–", "GECT–", "GE:2,3,4")
)

CROM <- list(
    ref = c("CROM:1 or Cra+)/CROM*01"),
    allele_name = c("/CROM.*"),
    hier = c("Cr(a–)/CROM*–01", "Tc(b+)/CROM*01.03", "Dr(a–)/CROM*01.–05", "WES(a+)/CROM*01.08", "UMC–/CROM*01.–10",
    "SERF–/CROM*01.–12", "CRUE−/CROM*01.–17", "Cra+)/CROM*01"),
    Crom_a = c("Tc(b+)", "Dr(a–)", "WES(a+)", "UMC–", "SERF–", "CRUE−", "Cra+)"),
    Crom_a_null = c("Cr(a–)")
)

# Too many "or", function will not work properly.
KN <- list(
    ref = c("KN:1 or Kn(a+) KN:3 or McC(a+) KN:4 or Sl1+ KN:5 or Yk(a+) KN:8 or Sl3+ KN:9 or KCAM+ KN:11 or DACY+/KN*01"),
    allele_name = c("/KN.*"),
    hier = c("KN:-4 or Sl1- KN:-3,6 or McC(a-b+) KN:7 or Vil+ KN:-9 or KCAM- KN:10 or KDAS+/KN*01.06", "KN:7 or Vil+ KN:-9 or KCAM- KN:10 or KDAS+/KN*01.07",
    "KN:2 or Kn(a-b+) KN:-9 or KCAM- KN:10 or KDAS+/KN*02.10", "KN:2 or Kn(a-b+)/KN*02", "KN:–5 or Yk(a–)/KN*01.-05", "KN:-9 or KCAM- KN:10 or KDAS+/KN*01.10",
    "KN:-11 or DACY- KN:12 or YCAD+/KN*01.12", "KN:1 or Kn(a+) KN:3 or McC(a+) KN:4 or Sl1+ KN:5 or Yk(a+) KN:8 or Sl3+ KN:9 or KCAM+ KN:11 or DACY+/KN*01"),
    Kn_a = c(""),
    Kn_b = c("")
)

IN <- list(
    ref = c("In(a–b+)/IN*02"),
    allele_name = c("/IN.*"),
    hier = c("In(a+b–)/IN*01", "INFI–/IN*02.–03", "In(a–b+)/IN*02"),
    In_a = c("In(a+b–)"),
    In_b = c("In(a–b+)", "INFI–")
)

OK <- list(
    ref = c("OK:1 or Ok(a+)/OK*01.01"),
    allele_name = c("/OK.*"),
    hier = c("Ok(a+)/OK*01.01"),
    Ok_a = c("Ok(a+)")
)

RAPH <- list(
    ref = c("RAPH:1 or MER2+/RAPH*01"),
    allele_name = c("/RAPH.*"),
    hier = c("MER2+/RAPH*01"),
    RAPH_1 = c("MER2+")
)

JMH <- list(
    ref = c("JMH:1 or JMH+/JMH*01"),
    allele_name = c("/JMH.*"),
    hier = c("JMH+/JMH*01"),
    Jmh = c("JMH")
)

I <- list(
    ref = c("I:1 or I+/GCNT2*01"),
    allele_name = c("/GCNT2.*"),
    hier = c("I+W/GCNT2*01W.01", "I+W/GCNT2*01W.02", "I+/GCNT2*01"),
    i = c("I+", "I+W"),
    i_null = c()
)

GLOB <- list(
    ref = c("GLOB:1 (P+)/GLOB*01"),
    allele_name = c("/GLOB.*"),
    hier = c("GLOB:–1 (P–)/GLOB*01N.01", "GLOB:1 (P+)/GLOB*01.02", "GLOB:1 (P+)/GLOB*01"),
    P = c("GLOB:1 (P+)"),
    P_null = c("GLOB:–1 (P–)")
)

GIL <- list(
    ref = c("GIL:1 or GIL+/GIL*01"),
    allele_name = c("/GIL.*"),
    hier = c("GIL+/GIL*01"),
    Gil = c("GIL+")
)

RHAG <- list(
    ref = c("RHAG:1 or Duclos+/RHAG*01"),
    allele_name = c("/RHAG.*"),
    hier = c("DSLK−, Kg+/RHAG*01.−03", "Rhmod/RHAG*01M.10", "Duclos+/RHAG*01"),
    Rhag = c("DSLK−, Kg+", "Rhmod", "Duclos+")
)

FORS <- list(
    ref = c("FORS:-1 (FORS-)/GBGT1*01N.01"),
    allele_name = c("/GBGT1.*"),
    hier = c("FORS:-1 (FORS-)/GBGT1*01N.04", "FORS:-1 (FORS-)/GBGT1*01N.02", "FORS:-1 (FORS-)/GBGT1*01N.03", "FORS:-1 (FORS-)/GBGT1*01N.01"),
    Fors_null = c("FORS:-1 (FORS-)")
)

# "Unclear phenotype" does not have proper allele_name.
JR <- list(
    ref = c("Jr(a+)/ABCG2*01"),
    allele_name = c("/ABCG2.*"),
    hier = c("Jr(a−)/ABCG2*01N.02.02", "Jr(a−)/ABCG2*01N.01", "Jr(a−)/ABCG2*01N.02.01", "Jr(a−)/ABCG2*01N.15", "Jr(a+w)/ABCG2*01W.01",
    "Jr(a+w)/ABCG2*01W.02", "Unclear Jra phnotype/N/A (1)", "Unclear Jra phnotype/N/A (6)", "Unclear Jra phnotype/N/A (7)", "Jr(a+)/ABCG2*01"),
    Jr_a = c("Jr(a+)", "Jr(a+w)"),
    Jr_null = c("Jr(a−)"),
    Jr_unclear = c("Unclear Jra phnotype")

)

LAN <- list(
    ref = c("Lan+/ABCB6*01"),
    allele_name = c("/ABCB6.*"),
    hier = c("Lan-/ABCB6*01N.02", "Lan-/ABCB6*01N.18", "Lan-/ABCB6*01N.26", "Lan-/ABCB6*01N.13", "Lan(+wk)/ABCB6*01W.01",
    "Lan(+wk)/ABCB6*01W.02", "Lan(+wk)/ABCB6*01W.03", "Lan(+wk)/ABCB6*01W.04", "Lan+/ABCB6*01"),
    Lan = c("Lan+", "Lan(+wk)"),
    Lan_null = c("Lan-")

)

VEL <- list(
    ref = c("VEL:1 (Vel+)/VEL*01"),
    allele_name = c("/VEL.*"),
    hier = c("VEL:1 (Vel+)/VEL*01"),
    Vel = c("VEL:1 (Vel+)")

)

CD59 <- list(
    ref = c("CD59:+1 or CD59.1+/CD59*01"),
    allele_name = c("/CD59.*"),
    hier = c("CD59.1+/CD59*01 "),
    Cd59 = c("CD59.1+")

)

AUG <- list(
    ref = c("AUG1+, At(a+), ATML−, ATAM+ AUG:1,2,−3,4/AUG*01"),
    allele_name = c("/AUG.*"),
    hier = c("AUG1+, At(a+), ATML−, ATAM+ AUG:1,2,−3,4"),
    Aug = c("AUG1+, At(a+), ATML−, ATAM+ AUG:1,2,−3,4")

)

GATA1 <- list(
    ref = c("Common/GATA1*01"),
    allele_name = c("/GATA1.*"),
    hier = c("Common/GATA1*01 "),
    Common = c("Common")

)

KLF1 <- list(
    ref = c("Common/KLF1*01"),
    allele_name = c("/KLF1.*"),
    hier = c("In(Lu)/KLF1*BGM13", "In(Lu)/KLF1*BGM20", "In(Lu)/KLF1*BGM23", "In(Lu)/KLF1*BGM32", "In(Lu)/KLF1*BGM33",
    "In(Lu)/KLF1*BGM34", "Common/KLF1*01"),
    Common = c("Common"),
    In_lu = c("In(Lu)")
)
# KN between CROM and IN
# JR between FORS and LAN
my_dict <- list(
    "ABO" = ABO, "P1PK" = P1PK, "LU" = LU, "KEL" = KEL, "LE" = LE, "FY" = FY, "JK" = JK, "DI" = DI, "YT" = YT, "XG" = XG,
    "SC" = SC, "DO" = DO, "CO" = CO, "LW" = LW, "CHRG_A" = CHRG_A, "CHRG_B" = CHRG_B, "H1" = H1, "H2" = H2, "KX" = KX, "GE" = GE, 
    "CROM" = CROM, "IN" = IN, "OK" = OK, "RAPH" = RAPH, "JMH" = JMH, "I" = I, "GLOB" = GLOB, "GIL" = GIL, "RHAG" = RHAG, "FORS" = FORS, 
    "LAN" = LAN, "VEL" = VEL, "CD59" = CD59, "AUG" = AUG, "GATA1" = GATA1, "KLF1" = KLF1, "JR" = JR
)

# Function to remove reference if there is more than 1 potential phenotypes
remove_ref <- function(row, ref) {
    substrings <- strsplit(row, " \\| ")[[1]]

    # Ensure that reference is kept when there are only erythrogene alleles
    if (length(substrings) > 1) {
        if (any(grepl("1000G only", substrings))) {
            count_1000G <- sum(grepl("1000G only", substrings))
            count_substrings <- length(substrings)
            if (count_substrings - count_1000G == 1) {
                return(substrings)
            }
        }
        # Handle H blood group cases as there are 2 reference allele FUT1.01 and FUT2.01
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
final_pheno <- function(row, dict, erythro) {
    substrings <- unlist(strsplit(row, " \\| "))
    
    for (pheno in dict[["hier"]]) {
        if (any(grepl(pheno, substrings, fixed = TRUE))) {
            pheno <- sub(dict[["allele_name"]], "", pheno)
            key <- find_key(dict, pheno)
            result <- paste0(key, ": ", pheno)
            if (length(erythro) > 0) {
                result <- paste0(result, " | ", paste(erythro, collapse = " | "))
            }
            return(result)
        }
    }
    return(row)
}

# Function that uses the functions above
process_group <- function(row, dict) {
    process_single_row <- function(single_row) {
        single_row <- remove_ref(single_row, dict[["ref"]])
        single_row <- sub(".* or ([^|]+).*", "\\1", single_row)
        erythro <- grep("1000G only", single_row, value = TRUE)
        return(final_pheno(single_row, dict, erythro))
    }
    return(sapply(row, process_single_row, USE.NAMES = FALSE))
}

# Set paths
PATH_SUM_TABLE <- "./R/results/summary_table.tsv"
df_sum_table <- read.delim(PATH_SUM_TABLE, header = TRUE)

df_sum_table <- df_sum_table %>%
    mutate(across(
        c("ABO", "P1PK", "LU", "KEL", "LE", "FY", "JK", "DI", "YT", "SC", "DO", "CO", "LW", 
          "KX", "GE", "CROM", "IN", "OK", "RAPH", "JMH", "I", "GLOB", "GIL", "RHAG", "FORS", 
          "LAN", "VEL", "CD59", "AUG", "GATA1", "KLF1", "JR"), 
        ~ process_group(.x, my_dict[[cur_column()]])
    )) %>%
    mutate(H = paste0(
        sub(" \\| .*", "", process_group(H, my_dict[["H1"]])), 
        " | ", 
        process_group(H, my_dict[["H2"]])
    ))

df_sum_table <- df_sum_table %>%
    mutate(
        LE = case_when(
            grepl("Le_null", LE) ~ "Le_null: Le(a-b-)",
            grepl("FUT3_active", LE) & grepl("H2_null", H) ~ "Le_a: Le(a+b-)",
            grepl("FUT3_active", LE) & grepl("H2_weak", H) ~ "Le_a_b: Le(a+b+)",
            TRUE ~ "Le_b: Le(a-b+)"
        ),
        
        ABO = ifelse(grepl("H1_null", H), "O: O (Due to H1 null)", ABO),
        
        GLOB = ifelse(grepl("p: p", P1PK), "P_null: P– (Due to p)", GLOB),
        
        LU = ifelse(grepl("In_lu", KLF1), "Lu_null: Lunull (Due to KLF1 Mutation)", LU),
        
        IN = ifelse(grepl("In_lu", KLF1), "In_null: In(a–b–) (Due to KLF1 Mutation)", IN)
    )

print("Done")