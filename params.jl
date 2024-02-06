# Parameters for PAM of a cell growing on a linear glucan with only 1 linkage

using JLD

#Random.seed!(20)  # Set a random seed for reproducibility

# cell parameters and limits
# membrane parameters
sa_TBDT = 0.01*0.005 # μm2, specific surface area of TBDT molecules PDB 8AA2
sa_UT = 0.005*0.01    # μm2, specific surface area of membrane transporter molecules (BioNumbers ID: 114685)
sa_lm = 0.4e-6     # μm2, specific surface area of membrane lipid molecules

lm_mw = 860         # g/mol, molecular weight of membrane lipid molecules
Mmax = 0.0045        # unitless, maximum membrane transporter to lipid ratio (counts/counts)
Mmin = 3e-6         # unitless, minimum membrane transporter to lipid ratio (counts/counts)

Av = 6.0221409e+23; #Avogadro's number (numbers/mole)
Dmax = 500; # g/L to accommodate (BioNumbers ID: 115317), maximum cell density 
Dmin = 200; # g/L, made up minimum cell density
r_cy_max = 0.6 # μm, maximum cytoplasm radius 
DNA_mw = 650 # g/mol per basepair
DNA_length = 2.5e6 # basepairs

# average carbon content of proteins
Qcp = 0.5; # g carbon/g protein 

# diffusion coefficient of substrate polysaccharide laminarin (DC_eP), laminarin oligosaccharides (DC_O), and periplasmic glucose (DC_S)
DC_Plam = 150 # μm^2/s, based on Elyakova et al. 1994
DC_Olam = 300 # μm^2/s, made up
DC_Glc = 600 # μm^2/s, BioNumbers
DC_Palg = 0.01 # μm^2/s Minami et al. 2010
DC_LOalg = 1 # μm^2/s 
DC_SOalg = 100 # μm^2/s
DC_Uro = 900 # μm^2/s, Liu et al. 2004

# extracellular boundary layer thickness
delta = 0.005 # μm, made up

# molecular weights of metabolites
Glc_mw = 180.156 # g/mol, molecular weight of glucose
Pyr_mw = 88.06 # g/mol, molecular weight of pyruvate
cC_mw = 12.0107 # g/mol, molecular weight of carbon
Plam_mw = 30 * 180.156 - 29 * 18 # g/mol, molecular weight of laminarin polysaccharide
Olam_mw = 6*180.156-5*18 # g/mol, molecular weight of laminarin oligosaccharide
dbOlam_mw = 5*180.156-4*18 # g/mol, molecular weight of debranched laminarin oligosaccharide
Palg_mw = 90000 # g/mol, molecular weight of alginate polysaccharide
LOalg_mw = 18000 # g/mol, molecular weight of large alginate oligosaccharides
SOalg_mw = 900 # g/mol, molecular weight of small alginate oligosaccharides
Uro_mw = 196 # g/mol, molecular weight of unsaturated uronate
ATP_mw = 507.181 # g/mol, molecular weight of ATP

# minima and maxima for modelled metabolites
eOlam_min = 0 # numbers, minimum extracellular laminarin oligosaccharide concentration
eOlam_max = 2e6 # numbers, maximum extracellular laminarin oligosaccharide concentration
pOlam_min = 0 # numbers, minimum periplasmic laminarin oligosaccharide concentration
pOlam_max = 2e6 # numbers, maximum periplasmic laminarin oligosaccharide concentration
pdbOlam_min = 0 # numbers, minimum periplasmic debranched laminarin oligosaccharide concentration
pdbOlam_max = 2e6 # numbers, maximum periplasmic debranched laminarin oligosaccharide concentration
eLO_min = 0 # numbers, minimum extracellular large alginate oligosaccharide concentration
eLO_max = 1e6 # numbers, maximum extracellular large alginate oligosaccharide concentration
eSO_min = 0 # numbers, minimum extracellular small alginate oligosaccharide concentration
eSO_max = 1e7 # numbers, maximum extracellular small alginate oligosaccharide concentration
pSO_min = 0 # numbers, minimum periplasmic small alginate oligosaccharide concentration
pSO_max = 1e7 # numbers, maximum periplasmic small alginate oligosaccharide concentration
pU_min = 0 # numbers, minimum periplasmic debranched alginate oligosaccharide concentration
pU_max = 5e6 # numbers, maximum periplasmic debranched alginate oligosaccharide concentration

cU_min = 0 # numbers, minimum periplasmic uronate concentration
cU_max = 2e7 # numbers, maximum periplasmic uronate concentration
pGlc_min = 0 # numbers, minimum periplasmic glucose concentration
pGlc_max = 2e7 # numbers, maximum periplasmic glucose concentration
cGlc_min = 0 # numbers, minimum intracellular glucose concentration
cGlc_max = 2* 0.01 * 1e-15 * Av # numbers, maximum intracellular glucose concentration: ~10 mM (BioNumbers ID: 104695) * 1 μm^3 max cytoplasm volume # MAY BE TOO LOW
cPyr_min = 0 # numbers, minimum intracellular pyruvate concentration
cPyr_max = 0.05 * 1e-15 * Av # numbers, maximum intracellular pyruvate concentration: ~50 mM * 1 μm^3 max cytoplasm volume
cC_min = 0 # numbers, minimum intracellular carbon concentration
cC_max = 0.3 * 3 * 1e-15 * Av # numbers, maximum intracellular carbon concentration: ~300 mM (BioNumbers ID: 104678) * 1 μm^3 max cytoplasm volume, assuing 3 carbon per molecule
cP_min = 0 # numbers, minimum intracellular protein concentration
cP_max = 2*500 * Qcp / cC_mw * 1e-15 * Av # 500 mg / mL volume of cytoplasm (BioNumbers ID: 115317) MAY BE TOO LOW
lm_min = 0 # numbers, minimum membrane lipid concentration
lm_max = 2*(4 * pi * (1)^2) / sa_lm # numbers, maximum membrane lipid concentration: 2 monolayers of 1 μm radius composed of only lipids

# stoichiometry matrix

            #   ENDOlam SENTalg ENDOalg IMPRTlam DEBRlam EXOlam IMPRTalg EXOalg UGlc CGlc GGlc UUro ACB CPyr  GPyr  P    M    
stoich_matrix = [-0.2     0       0        0       0       0       0       0     0    0    0    0    0    0    0    0    0      # extracellular laminarin polysaccharide ... ePlam (DP30)
                  1       0       0       -1       0       0       0       0     0    0    0    0    0    0    0    0    0      # extracellular laminarin oligosaccharide ... eOlam (DP6)
                  0      -0.2     0        0       0       0       0       0     0    0    0    0    0    0    0    0    0      # extracellular alginate polysaccharide ... eP (DP500, ~90,000 Da)
                  0       0.5    -0.05     0       0       0       0       0     0    0    0    0    0    0    0    0    0      # extracellular M-rich large alginate oligosaccharide ... eLMO (DP100, ~18,000 Da)
                  0       0.5    -0.05     0       0       0       0       0     0    0    0    0    0    0    0    0    0      # extracellular G-rich large alginate oligosaccharide ... eLGO (DP100, ~18,000 Da)
                  0       0       1        0       0       0      -1       0     0    0    0    0    0    0    0    0    0      # extracellular M-rich small alginate oligosaccharide ... eSMO (DP5, ~900 Da)
                  0       0       1        0       0       0      -1       0     0    0    0    0    0    0    0    0    0      # extracellular G-rich small alginate oligosaccharide ... eSGO (DP5, ~900 Da)
                  0       0       0        1      -1       0       0       0     0    0    0    0    0    0    0    0    0      # periplasmic laminarin oligosaccharide ... pOlam
                  0       0       0        0       1      -0.2     0       0     0    0    0    0    0    0    0    0    0      # periplasmic debranched laminarin oligosaccharide ... pdbOlam (DP5)
                  0       0       0        0       1       1       0       0    -1    0    0    0    0    0    0    0    0      # periplasmic glucose ... pGlc
                  0       0       0        0       0       0       1      -0.25  0    0    0    0    0    0    0    0    0      # periplasmic M-rich alginate oligosaccharide ... pSMO
                  0       0       0        0       0       0       1      -0.25  0    0    0    0    0    0    0    0    0      # periplasmic G-rich alginate oligosaccharide ... pSGO
                  0       0       0        0       0       0       0       1     0    0    0   -1    0    0    0    0    0      # periplasmic unsaturated uronate ... pU
                  0       0       0        0       0       0       0       0     1   -1   -1    0    0    0    0    0    0      # intracellular glucose ... cGlc
                  0       0       0        0       0       0       0       0     0    0    0    1   -1    0    0    0    0      # intracellular unsaturated uronate ... cU
                  0       0       0        0       0       0       0       0     0    0    0    0    1   -1   -1    0    0      # intracellular pyruvate ... cPyr
                  0       0       0        0       0       0       0       0     0    4    0    0    0    2    0   -1   -1      # cellular carbon ... cC
                  0       0       0       -0.25    0       0      -0.25    0    -2    7.2 30   -2    0    3.6 15  -13   -4      # energy in ATP ... ATP (based on Su and refs therein, Rees et al. 2010 10.1038/nrm2646)
                  0       0       0        0       0       0       0       0     0    0    0    0    0    0    0    1    0      # carbon atoms in proteins ... cP
                  0       0       0        0       0       0       0       0     0    0    0    0    0    0    0    0    0.02 ] # membrane lipids ... lm


# reaction rates CAZymes
# based on normal distributions for endo and exo enzymes

mean_log_kM_endo = -6.48 # M, mean kM for endo enzymes
sd_log_kM_endo = 2       # M, standard deviation kM for endo enzymes, times unexplained variance
sd_log_kcat_endo_kMdpdnt = 1.19 # 1/s, standard deviation kcat for endo enzymes, kM dependent

mean_log_kcat_exo = 3.75 # 1/s, mean kcat for exo enzymes
sd_log_kcat_exo = 1.28 # 1/s, standard deviation kcat for exo enzymes
mean_log_kM_exo = -6.39 # M, mean kM for exo enzymes
sd_log_kM_exo = 1.66 # M, standard deviation kM for exo enzymes

# laminarin cascade
# extracellular endolaminarases GH16
#GH16A_km = exp(mean_log_kM_endo + sd_log_kM_endo * randn()) # M, devide by DP of substrate to relate to monomer
#GH16A_kcat = exp(log(GH16A_km) * 0.5 + 5.8 + sd_log_kcat_endo_kMdpdnt * randn())  # 1/s
#GH16A_mw = 30000 # g/mol
#GH16B_km = exp(mean_log_kM_endo + sd_log_kM_endo * randn()) # M, devide by DP of substrate to relate to monomer
#GH16B_kcat = exp(log(GH16B_km) * 0.5 + 5.8 + sd_log_kcat_endo_kMdpdnt * randn())  # 1/s
#GH16B_mw = 30000 # g/mol
#GH16C_km = exp(mean_log_kM_endo + sd_log_kM_endo * randn()) # M, devide by DP of substrate to relate to monomer
#GH16C_kcat = exp(log(GH16C_km) * 0.5 + 5.8 + sd_log_kcat_endo_kMdpdnt * randn())  # 1/s
#GH16C_mw = 30000 # g/mol

# periplasmic debranching laminarases GH30
GH30A_kcat = exp(mean_log_kcat_exo + sd_log_kcat_exo * randn())  # 1/s
GH30A_km = exp(mean_log_kM_exo + sd_log_kM_exo * randn()) # M, devide by DP of substrate to relate to monomer
GH30A_mw = 55000 # g/mol
GH30B_kcat = exp(mean_log_kcat_exo + sd_log_kcat_exo * randn())  # 1/s
GH30B_km = exp(mean_log_kM_exo + sd_log_kM_exo * randn()) # M, devide by DP of substrate to relate to monomer
GH30B_mw = 55000 # g/mol

# periplasmic exolaminarases GH17
GH17A_kcat = exp(mean_log_kcat_exo + sd_log_kcat_exo * randn())  # 1/s
GH17A_km = exp(mean_log_kM_exo + sd_log_kM_exo * randn()) # M, devide by DP of substrate to relate to monomer
GH17A_mw = 90000 # g/mol
GH17B_kcat = exp(mean_log_kcat_exo + sd_log_kcat_exo * randn())  # 1/s
GH17B_km = exp(mean_log_kM_exo + sd_log_kM_exo * randn()) # M, devide by DP of substrate to relate to monomer
GH17B_mw = 90000 # g/mol

# alginate cascade
# extracellular sentinel endo-alginate lyases PL7
#PL7A_km = exp(mean_log_kM_endo + sd_log_kM_endo * randn()) # M, devide by DP of substrate to relate to monomer
#PL7A_kcat = exp(log(PL7A_km) * 0.5 + 5.8 + sd_log_kcat_endo_kMdpdnt * randn())  # 1/s
#PL7A_mw = 40000 # g/mol
#PL7B_km = exp(mean_log_kM_endo + sd_log_kM_endo * randn()) # M, devide by DP of substrate to relate to monomer
#PL7B_kcat = exp(log(PL7B_km) * 0.5 + 5.8 + sd_log_kcat_endo_kMdpdnt * randn())  # 1/s
#PL7B_mw = 40000 # g/mol
#PL7C_km = exp(mean_log_kM_endo + sd_log_kM_endo * randn()) # M, devide by DP of substrate to relate to monomer
#PL7C_kcat = exp(log(PL7C_km) * 0.5 + 5.8 + sd_log_kcat_endo_kMdpdnt * randn())  # 1/s
#PL7C_mw = 40000 # g/mol

# extracellular M-specific endo-alginate lyases PL5
PL5A_km = exp(mean_log_kM_endo + sd_log_kM_endo * randn()) # M, devide by DP of substrate to relate to monomer
PL5A_kcat = exp(log(PL5A_km) * 0.5 + 5.8 + sd_log_kcat_endo_kMdpdnt * randn())  # 1/s
PL5A_mw = 35000 # g/mol
PL5B_km = exp(mean_log_kM_endo + sd_log_kM_endo * randn()) # M, devide by DP of substrate to relate to monomer
PL5B_kcat = exp(log(PL5B_km) * 0.5 + 5.8 + sd_log_kcat_endo_kMdpdnt * randn())  # 1/s
PL5B_mw = 35000 # g/mol
PL5C_km = exp(mean_log_kM_endo + sd_log_kM_endo * randn()) # M, devide by DP of substrate to relate to monomer
PL5C_kcat = exp(log(PL5C_km) * 0.5 + 5.8 + sd_log_kcat_endo_kMdpdnt * randn())  # 1/s
PL5C_mw = 35000 # g/mol

# extracellular G-specific endo-alginate lyases PL38
PL38A_km = exp(mean_log_kM_endo + sd_log_kM_endo * randn()) # M, devide by DP of substrate to relate to monomer
PL38A_kcat = exp(log(PL38A_km) * 0.5 + 5.8 + sd_log_kcat_endo_kMdpdnt * randn())  # 1/s
PL38A_mw = 46000 # g/mol
PL38B_km = exp(mean_log_kM_endo + sd_log_kM_endo * randn()) # M, devide by DP of substrate to relate to monomer
PL38B_kcat = exp(log(PL38B_km) * 0.5 + 5.8 + sd_log_kcat_endo_kMdpdnt * randn())  # 1/s
PL38B_mw = 46000 # g/mol
PL38C_km = exp(mean_log_kM_endo + sd_log_kM_endo * randn()) # M, devide by DP of substrate to relate to monomer
PL38C_kcat = exp(log(PL38C_km) * 0.5 + 5.8 + sd_log_kcat_endo_kMdpdnt * randn())  # 1/s
PL38C_mw = 46000 # g/mol

# periplasmic M-specific exo-alginate lyases PL17
PL17A_kcat = exp(mean_log_kcat_exo + sd_log_kcat_exo * randn())  # 1/s
PL17A_km = exp(mean_log_kM_exo + sd_log_kM_exo * randn()) # M, devide by DP of substrate to relate to monomer
PL17A_mw = 80500 # g/mol
PL17B_kcat = exp(mean_log_kcat_exo + sd_log_kcat_exo * randn())  # 1/s
PL17B_km = exp(mean_log_kM_exo + sd_log_kM_exo * randn()) # M, devide by DP of substrate to relate to monomer
PL17B_mw = 80500 # g/mol

# periplasmic G-specific exo-alginate lyase PL6
PL6A_kcat = exp(mean_log_kcat_exo + sd_log_kcat_exo * randn())  # 1/s
PL6A_km = exp(mean_log_kM_exo + sd_log_kM_exo * randn()) # M, devide by DP of substrate to relate to monomer
PL6A_mw = 85200 # g/mol
PL6B_kcat = exp(mean_log_kcat_exo + sd_log_kcat_exo * randn())  # 1/s
PL6B_km = exp(mean_log_kM_exo + sd_log_kM_exo * randn()) # M, devide by DP of substrate to relate to monomer
PL6B_mw = 85200 # g/mol

# reaction rates core enyzmes
# TBDT in outer membrane (oligosaccharide uptake)
TBDT_kcat = 0.1 # 1/s
TBDT_km = 0.0008  # M, Balhesteros et al. 2017 https://doi.org/10.1128/jb.00723-16
TBDT_mw = 190000 # g/mol, Silale and Van den Berg 2023 10.1146/annurev-micro-032421-111116 and Sverzhinsky et al. 2015 https://doi.org/10.1128%2FJB.00069-15 (Chimento et al. 2003 https://doi.org/10.1038/nsb914)

# monosaccharide uptake
UT_kcat = 5  # 1/s, questionable (from Garcia-Alles et al. 2002 doi.org/10.1021/bi025928d), updated from Walmsey et al. 1998 10.1016/S0968-0004(98)01326-7
UT_km = 0.0001  # M, Garcia-Alles et al. 2002 doi.org/10.1021/bi025928d, updated from Walmsey et al. 1998 10.1016/S0968-0004(98)01326-7
UT_mw = 176000 # g/mol, adapted from Norris et al (dismissed https://biocyc.org/ECOLI/NEW-IMAGE?type=NIL&object=CPLX-157)

# alginate-catabolism bridge
ACB_kcat = 100.0  # 1/s
ACB_km = 0.0001 # M
ACB_mw = 120000 # g/mol

# catabolism
C_kcat = 100.0  # 1/s, Cerisy et al. 2019 https://doi.org/10.1128/jb.00241-19 and Norris et al
C_km = 0.0001 # M, ?! adapted from Norris et al
C_mw = 4400000 # g/mol, adapted from Norris et al

# protein synthesis
P_kcat = 15.0  # 1/s, adapted from Norris et al
P_km = 0.0003 # M, ?! adapted from Norris et al
P_mw = 2200000 # g/mol

# membrane synthesis
M_kcat = 3.0  # 1/s, adapted from Norris et al 
M_km = 0.0002 # M, ?! adapted from Norris et al
M_mw = 2310000 # g/mol, adapted from Norris et al


