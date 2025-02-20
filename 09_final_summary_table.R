# R script to infer phenotype based on variant positions
rm(list = ls())

setwd("/home/svu/e0726996/Blood-type-GWAS/")

# Set random seed
set.seed(0)

# Create "dictionary" for each blood group
# ref: reference allele of blood group
# allele_name: pattern used to clean up string at the end
# hier: Alleles found in df_sum_table, ranked based on # nucleotide changes associated to allele (Largest -> Smallest)
# Need to add in remaining alleles!!!!!
# Others: Provide a generalised name for phenotypes For example: "AB: cis(AB)"
ABO <- list(
    ref = c("A1/ABO*A1.01"),
    allele_name = c("/ABO.*"),
    hier = c("Aweak/ABO*AW.09", "A2/ABO*A2.09", "Aweak/ABO*AW.25", "A2/ABO*A2.01", "Aweak/ABO*AW.27",
    "A3/ABO*A3.02", "A2/ABO*A2.06", "O/ABO*O.10", "B/ABO*B.02", "B3/ABO*B3.02", "B/ABO*B.01",
    "A2/ABO*A2.04", "cisAB/ABO*cisAB.05", "cisAB/ABO*cisAB.06", "B(A)/ABO*BA.03", "O/ABO*O.02.03",
    "B(A)/ABO*BA.01", "O/ABO*O.02.01", "Ax/Aweak/ABO*AW.31.01", "Ax/Aweak/ABO*AW31.02- 05",
    "Ael/ABO*AEL.07", "cisAB/ABO*cisAB.02", "Ael/ABO*AEL.02", "A2/ABO*A2.20", "Ax/Aweak/ABO*AW.30.02",
    "A1/ABO*A1.02", "A2/ABO*A2.02", "A2/ABO*A2.03", "Ax/Aweak/ABO*AW.30.01", "O/ABO*O.06", "O/ABO*O.13", "A1/ABO*A1.01"),
    AB = c("cisAB", "B(A)"),
    A = c("A1", "A2"),
    B = c("B"),
    ABO_weak = c("A3", "Aweak", "Ax/Aweak", "Ael", "B3", "Bel"),
    O = c("O")
)

P1PK <- list(
    ref = c("P1+ Pk+/A4GALT*01"),
    allele_name = c("/A4GALT.*"),
    hier = c("P1+ Pk+/A4GALT*01.02", "p/A4GALT*01N.16", "P1+ Pk+/A4GALT*01"),
    P1 = c("P1+ Pk+", "P1+ Pk+ NOR+"),
    P2 = c("P1– Pk+ (P2)"),
    p = c("p")
)

RHCE <- list()

RHD <- list()

LU <- list(
    ref = c("LU:2 or Lu(b+)/LU*02"),
    allele_name = c("/LU.*"),
    hier = c("LU:-13/LU*02.-13", "LU:-16/LU*01.-16", "LU:1,19/LU*01.19", "Au(a−b+)/LU*02.19", "LU:1 or Lu(a+)/LU*01",
    "LU:-5/LU*02.-05", "LU:-6,9/LU*02.09", "LU:-8,14/LU*02.14", "LU:-22, LURC−/LU*02.–22", "LU:-26, LUBI−/LU*02.-26",
    "LU:-29, LURA−/LU*02.-29", "Lunull/LU*02N.03", "LU:2 or Lu(b+)/LU*02"),
    Lu_a = c("Lu(a+)", "LU:-16", "LU:1,19"),
    Lu_b = c("Lu(b+)", "LU:-5", "LU:-6,9", "LU:-8,14", "Au(a−b+)", "LU:-22, LURC−", "LU:-26, LUBI−", "LU:-29, LURA−"),
    Lu_null = c("Lunull")
)

KEL <- list(
    ref = c("KEL:2 or k+/KEL*02"),
    allele_name = c("/KEL.*"),
    hier = c("K0/KEL*02N.11", "K0/KEL*02N.44", "Kmod/KEL*02M.04", "K0 phenotype/KEL*02M.05", "K+k-/KEL*01.01", "k+/KEL*02.02",
    "Kp(a+b–c–)/KEL*02.03", "Js(a+b–)/KEL*02.06", "Ul(a+)/KEL*02.10", "KEL: -11,17/KEL*02.17", "KEL:-18/KEL*02.-18.1",
    "KEL:-18/KEL*02.-18.2", "KEL:-19/KEL*02.-19", "TOU-/KEL*02.-26", "KTIM-/KEL*02.-30", "KYO+, KYOR-/ KEL*02.31",
    "KETI-/KEL*02.-36", "KHUL-,KEAL+/KEL*02.39", "KHIZ-, KHOZ+/KEL*02.41", "K0/KEL*02N.12", "K0/KEL*02N.13", "K0/KEL*02N.16",
    "K0/KEL*02N.17", "K0/KEL*02N.39", "K0/KEL*02N.42", "KEL:1weak/KEL*01M.01", "k+/KEL*02"),
    K = c("K+k-"),
    k = c("k+", "Kp(a+b–c–)", "Js(a+b–)", "Ul(a+)", "KEL: -11,17", "KEL:-18", "KEL:-19", "TOU-", "KTIM-", "KYO+, KYOR-", "KETI-", 
    "KHUL-,KEAL+", "KHIZ-, KHOZ+"),
    Kel_null = c("K0"),
    Kel_weak = c("KEL:1weak", "Kmod", "K0 phenotype") # Affects both K and k
)

LE <- list(
    ref = c("FUT3 Active/FUT3*01.01"),
    allele_name = c("/FUT3.*"),
    hier = c("FUT3 Inactive - Le(a-b-)/FUT3*01N.03.12", "FUT3 Inactive - Le(a-b-)/FUT3*01N.01.08", "FUT3 Inactive - Le(a-b-)/FUT3*01N.01.09",
    "FUT3 Inactive - Le(a-b-)/FUT3*01N.01.10", "FUT3 Inactive - Le(a-b-)/FUT3*01N.01.11", "FUT3 Inactive - Le(a-b-)/FUT3*01N.17.16",
    "FUT3 Active/FUT3*01.16", "FUT3 Inactive - Le(a-b-)/FUT3*01N.01.02", "FUT3 Inactive - Le(a-b-)/FUT3*01N.01.05", "FUT3 Inactive - Le(a-b-)/FUT3*01N.01.12",
    "FUT3 Inactive - Le(a-b-)/FUT3*01N.03.05", "FUT3 Inactive - Le(a-b-)/FUT3*01N.03.08", "FUT3 Inactive - Le(a-b-)/FUT3*01N.03.11",
    "FUT3 Inactive - Le(a-b-)/FUT3*01N.17.07", "FUT3 Inactive - Le(a-b-)/FUT3*01N.17.12", "FUT3 Inactive - Le(a-b-)/FUT3*01N.17.15",
    "FUT3 Inactive - Le(a-b-)/FUT3*01N.17.17", "FUT3 Inactive - Le(a-b-)/FUT3*01N.17.19", "FUT3 Active/FUT3*01.15", "FUT3 Inactive - Le(a-b-)/FUT3*01N.01.01",
    "FUT3 Inactive - Le(a-b-)/FUT3*01N.03.02", "FUT3 Inactive - Le(a-b-)/FUT3*01N.03.04", "FUT3 Inactive - Le(a-b-)/FUT3*01N.17.02",
    "FUT3 Inactive - Le(a-b-)/FUT3*01N.17.03", "FUT3 Inactive - Le(a-b-)/FUT3*01N.17.04", "FUT3 Inactive - Le(a-b-)/FUT3*01N.17.05",
    "FUT3 Active/FUT3*01.03", "FUT3 Active/FUT3*01.04", "FUT3 Active/FUT3*01.05", "FUT3 Active/FUT3*01.06", "FUT3 Active/FUT3*01.10",
    "FUT3 Active/FUT3*01.11", "FUT3 Active/FUT3*01.13", "FUT3 Active/FUT3*01.17.01", "FUT3 Inactive - Le(a-b-)/FUT3*01N.03.01",
    "FUT3 Active/FUT3*01.01"),
    FUT3_active = c("FUT3 Active"),
    Le_null = c("FUT3 Inactive - Le(a-b-)")
)

FY <- list(
    ref = c("FY:1 or Fy(a+)/FY*01"),
    allele_name = c("/FY.*"),
    hier = c("Fy(b+w), Fyx/FY*02W.01", "Fy(a−b−)/FY*02N.01", "Fy(a+w)/FY*01W.02", "Fy(b+)/FY*02", "Fy(a+w)/FY*01W.01", "Fy(a−b−)/FY*01N.01",
    "Fy(a−b−)/FY*01N.06", "Fy(a+)/FY*01"),
    Fy_a = c("Fy(a+)", "Fy(a+w)"),
    Fy_b = c("Fy(b+)", "Fy(b+w), Fyx"),
    Fy_null = c("Fy(a−b−)")
)

JK <- list(
    ref = c("JK:1 or Jk(a+)/JK*01"),
    allele_name = c("/JK.*"),
    hier = c("Jk(a–b–)/JK*01N.20", "Jk(b+w)/JK*02W.03", "Jk(a–b–)/JK*02N.13", "Jk(b+W)/JK*02W.04", "Jk(a+w)/JK*01W.06", "Jk(a+W)/JK*01W.11",
    "Jk(a–b–)/JK*01N.17", "Jk(a–b–)/JK*02N.01", "Jk(a–b–)/JK*02N.02", "Jk(a–b–)/JK*02N.06", "Jk(a–b–)/JK*02N.07", "Jk(a–b–)/JK*02N.08",
    "Jk(a–b–)/JK*02N.09", "Jk(a–b–)/JK*02N.17", "Jk(a–b–)/JK*02N.21", "Jk(b+)/JK*02", "Jk(a+W)/JK*01W.01", "Jk(a+W)/JK*01W.02",
    "Jk(a+W)/JK*01W.03", "Jk(a+W)/JK*01W.04", "Jk(a+W)/JK*01W.09", "Jk(a–b–)/JK*01N.03", "Jk(a–b–)/JK*01N.04", "Jk(a–b–)/JK*01N.06",
    "Jk(a–b–)/JK*01N.18", "Jk(a–b–)/JK*01N.19", "Jk(a+)/JK*01"),
    Jk_a = c("Jk(a+)", "Jk(a+W)"),
    Jk_b = c("Jk(b+)", "Jk(b+W)"),
    Jk_null = c("Jk(a–b–)")
)

DI <- list(
    ref = c("DI:–1,2 or Di(a–b+)/DI*02"),
    allele_name = c("/DI.*"),
    hier = c("Di(a+b–)", "Di(a–b+)"),
    Di_a = c("Di(a+b–)"),
    Di_b = c("Di(a–b+)")
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
    Xg_a = c("Xga/XG*01")
)

SC <- list(
    ref = c("SC:1 or Sc1+/SC*01"),
    allele_name = c("/SC.*"),
    hier = c("Sc2+/SC*02", "SCAC–/SC*01.–09", "STAR–/SC*01.–05", "SCAN–/SC*01.–07", "Sc1+/SC*01"),
    Sc_1 = c("SCAC–", "STAR–", "SCAN–", "Sc1+"),
    Sc_2 = c("Sc2+")
)

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
    hier = c("H–/FUT1*01N.21", "H+weak/FUT1*01W.32", "H+weak/FUT1*01W.02", "H+weak/FUT1*01W.09", "H+weak/FUT1*01W.24", "H+/FUT1*01.02",
    "H–/FUT1*01N.09", "H+/FUT1*01"),
    H1 = c("H+"),
    H1_weak = c("H+weak"),
    H1_null = c("H–")
)

H2 <- list(
    ref = c("H+/FUT2*01"),
    allele_name = c("/FUT2.*"),
    hier = c("H+/FUT2*01.03.02", "H+/FUT2*01.02", "H+/FUT2*01.03.01", "H+/FUT2*01.04", "H+/FUT2*01.06", "H+/FUT2*01.08", "H+w/FUT2*01W.01",
    "H+w/FUT2*01W.02.01", "H–/FUT2*01N.02", "H–/FUT2*01N.03", "H–/FUT2*01N.04", "H–/FUT2*01N.05", "H–/FUT2*01N.06", "H–/FUT2*01N.12",
    "H–/FUT2*01N.13", "H–/FUT2*01N.14", "H–/FUT2*01N.15", "H–/FUT2*01N.16", "H–/FUT2*01N.17", "H+/FUT2*01"),
    H2 = c("H+"),
    H2_weak = c("H+w"),
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
    hier = c("Cr(a–)/CROM*–01", "Tc(b+)/CROM*01.03", "Dr(a–)/CROM*01.–05", "WES(a+)/CROM*01.08", "UMC–/CROM*01.–10", "SERF–/CROM*01.–12",
    "CRUE−/CROM*01.–17", "Cra+)/CROM*01"),
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
    Kn = c("")
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
    i = c("I+"),
    i_weak = c("I+W")
)

GLOB <- list(
    ref = c("GLOB:1 (P+)/GLOB*01"),
    allele_name = c("/GLOB.*"),
    hier = c("GLOB:1 (P+)/GLOB*01.02", "GLOB:–1 (P–)/GLOB*01N.01", "GLOB:1 (P+)/GLOB*01"),
    Glob_P = c("GLOB:1 (P+)"),
    Glob_P_null = c("GLOB:–1 (P–)")
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
    Rhag = c("DSLK−", "Duclos+"),
    Rhag_weak = c("Rhmod")
)

FORS <- list(
    ref = c("FORS:-1 (FORS-)/GBGT1*01N.01"),
    allele_name = c("/GBGT1.*"),
    hier = c("FORS:-1 (FORS-)/GBGT1*01N.04", "FORS:-1 (FORS-)/GBGT1*01N.02", "FORS:-1 (FORS-)/GBGT1*01N.03", "FORS:-1 (FORS-)/GBGT1*01N.01"),
    Fors_null = c("FORS:-1 (FORS-)")
)

# "Unclear phenotype" does not have proper allele_name. Function will not work properly.
JR <- list(
    ref = c("Jr(a+)/ABCG2*01"),
    allele_name = c("/ABCG2.*"),
    hier = c("Jr(a−)/ABCG2*01N.15", "Jr(a−)/ABCG2*01N.02.02", "Jr(a−)/ABCG2*01N.01", "Jr(a−)/ABCG2*01N.02.01", "Jr(a+w)/ABCG2*01W.01",
    "Jr(a+w)/ABCG2*01W.02", "Unclear Jra phnotype/N/A (1)", "Unclear Jra phnotype/N/A (6)", "Unclear Jra phnotype/N/A (7)", "Jr(a+)/ABCG2*01"),
    Jr_a = c("Jr(a+)"),
    Jr_null = c("Jr(a−)"),
    Jr_unclear = c("Unclear Jra phnotype")

)

LAN <- list(
    ref = c("Lan+/ABCB6*01"),
    allele_name = c("/ABCB6.*"),
    hier = c("Lan-/ABCB6*01N.02", "Lan-/ABCB6*01N.18", "Lan-/ABCB6*01N.26", "Lan-/ABCB6*01N.13", "Lan(+wk)/ABCB6*01W.01",
    "Lan(+wk)/ABCB6*01W.02", "Lan(+wk)/ABCB6*01W.03", "Lan(+wk)/ABCB6*01W.04", "Lan+/ABCB6*01"),
    Lan = c("Lan+"),
    Lan_null = c("Lan-"),
    Lan_weak = c("Lan(+wk)")

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
    AUG_1 = c("AUG1+, At(a+), ATML−, ATAM+ AUG:1,2,−3,4")

)

GATA1 <- list(
    ref = c("Common/GATA1*01 "),
    allele_name = c("/GATA1.*"),
    hier = c("Common/GATA1*01 "),
    Common = c("Common")

)

KLF1 <- list(
    ref = c("Common/GATA1*01 "),
    allele_name = c("/GATA1.*"),
    hier = c("Common/GATA1*01 "),
    Common = c("Common")

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
            count_FUT1 <- sum(grepl("FUT1\\*", substrings))
            if (count_FUT1 == 1) {
                return(substrings)
            }
        }
        if (any(grepl("FUT2\\*", ref))) {
            count_FUT2 <- sum(grepl("FUT2\\*", substrings))
            if (count_FUT2 == 1) {
                return(substrings)
            }
        }
    }
    # Do nothing if there is only reference
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
    for (pheno in dict[["hier"]]) {
        if (pheno %in% row) {
            pheno <- sub(dict[["allele_name"]], "", pheno)
            row <- find_key(dict, pheno)
            row <- paste0(row, ": ", pheno)
            if (length(erythro) > 0) {
                row <- paste0(row, " | ", paste(erythro, collapse = " | "))
            }
            return(row)
        }
    }
    return(row)
}

# Function that uses the functions above
process_group <- function(row, dict) {
    row <- remove_ref(row, dict[["ref"]])
    simplified_substrings <- sub(".* or ", "", row)
    erythro <- grep("1000G only", simplified_substrings, value = TRUE)
    row <- final_pheno(unique(simplified_substrings), dict, erythro)
    return(row)
}

# Set paths
PATH_SUM_TABLE <- "./R/results/summary_table.tsv"
df_sum_table <- read.delim(PATH_SUM_TABLE, header = TRUE)

for (i in seq(1, nrow(df_sum_table))) {
    for (col in colnames(df_sum_table)) {
        row <- df_sum_table[[col]][i]
        
        df_sum_table[[col]][i] <- switch(col,
        "ABO" = process_group(row, ABO),
        #"MNS" = {},
        "P1PK" = process_group(row, P1PK),
        #"RH" = {},
        "LU" = process_group(row, LU),
        "KEL" = process_group(row, KEL),
        "LE" = process_group(row, LE),
        "FY" = process_group(row, FY),
        "JK" = process_group(row, JK),
        "DI" = process_group(row, DI),
        "YT" = process_group(row, YT),
        # "XG" = process_group(row, XG),
        "SC" = process_group(row, SC),
        "DO" = process_group(row, DO),
        "CO" = process_group(row, CO),
        "LW" = process_group(row, LW),
        # "CRHG" = {}",
        "H" = paste0(sub(" \\| .*", "", process_group(row, H1)), " | ", process_group(row, H2)),
        "KX" = process_group(row, KX),
        "GE" = process_group(row, GE),
        "CROM" = process_group(row, CROM),
        #"KN" = process_group(row, KN),
        "IN" = process_group(row, IN),
        "OK" = process_group(row, OK),
        "RAPH" = process_group(row, RAPH),
        "JMH" = process_group(row, JMH),
        "I" = process_group(row, I),
        "GLOB" = process_group(row, GLOB),
        "GIL" = process_group(row, GIL),
        "RHAG" = process_group(row, RHAG),
        "FORS" = process_group(row, FORS),
        #"JR" = process_group(row, JR),
        "LAN" = process_group(row, LAN),
        "VEL" = process_group(row, VEL),
        "CD59" = process_group(row, CD59),
        "AUG" = process_group(row, AUG),
        "GATA1" = process_group(row, GATA1),
        "KLF1" = process_group(row, KLF1),
        row
    )
    }
}

for (i in seq(1, nrow(df_sum_table))) {
    cell_ABO <- df_sum_table[["ABO"]][i]
    cell_LE <- df_sum_table[["LE"]][i]
    cell_H <- df_sum_table[["H"]][i]
            
    if (any(grepl("Le_null", cell_LE))) { # FUT3- -> Le(a-b-)
        df_sum_table[["LE"]][i] <- "Le_null: Le(a-b-)"
    }

    else if (any(grepl("FUT3_active", cell_LE)) && any(grepl("H2_null", cell_H))) { # FUT2-, FUT3+ -> Le(a+b-)
        df_sum_table[["LE"]][i] <- "Le_a: Le(a+b-)"
    }

    else if (any(grepl("FUT3_active", cell_LE)) && any(grepl("H2_weak", cell_H))) { # weak FUT2, FUT3+ -> Le(a+b+)
        df_sum_table[["LE"]][i] <- "Le_a_b: Le(a+b+)"
    }

    else { # FUT2+, FUT3+ -> Le(a-b+)
        df_sum_table[["LE"]][i] <- "Le_b: Le(a-b+)"
    }
}


print("Done")


# Ignore
# MNS_MN <- list(
#     ref = c("MNS:1 or M+/GYPA*01 or GYPA*M"),
#     allele_name = c("/GYPA.*"),
#     hier = c("Vr+", "Mc", "Mt(a+)", "M+", "M+w", "N"),
#     M = c("Vr+", "Mc", "Mt(a+)", "M+", "M+w"),
#     N = ("N+")
# )

# MNS_Ss <- list(
#     ref = c("MNS:4 or s+/GYPB*04 or GYPB*s"),
#     allele_name = c("/GYPB.*"),
#     hier = c("SD+", "S+", "s+", "S+w", "s+w", "MNS:4,5 (s+, U+) Altered GPB", "s+w, U+w", "Mv+", "Mit+", "s-U+var", "S-U+var", "S-s-U-"),
#     Ss = c("SD+"),
#     S = c("S+", "S+w"),
#     s = c("s+", "s+w", "MNS:4,5 (s+, U+) Altered GPB", "s+w, U+w"),
#     MNS_others = c("Mv+", "Mit+", "s-U+var", "S-U+var", "S-s-U-")
# )