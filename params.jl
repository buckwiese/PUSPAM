# Parameters for PAM of a cell growing on a linear glucan with only 1 linkage

using JLD

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
DC_eP = 150 # μm^2/s, based on Elyakova et al. 1994
DC_O = 300 # μm^2/s, made up
DC_S = 600 # μm^2/s, BioNumbers

# extracellular boundary layer thickness
delta = 0.005 # μm, made up

# molecular weights of metabolites
S_c_mw = 180.156 # g/mol, molecular weight of intracellular substrate
C_c_mw = 12.0107 # g/mol, molecular weight of cellular carbon
eP_mw = 30 * 180.156 - 29 * 18 # g/mol, molecular weight of extracellular polysaccharide
pO_mw = 6*180.156-5*18 # g/mol, molecular weight of periplasmic oligosaccharide
pdbO_mw = 5*180.156-4*18 # g/mol, molecular weight of periplasmic debranched oligosaccharide
ATP_mw = 507.181 # g/mol, molecular weight of ATP

# minima and maxima for modelled metabolites
eO_min = 0 # numbers, minimum extracellular oligosaccharide concentration
eO_max = 2e6 # numbers, maximum extracellular oligosaccharide concentration
pO_min = 0 # numbers, minimum periplasmic oligosaccharide concentration
pO_max = 2e6 # numbers, maximum periplasmic oligosaccharide concentration
pdbO_min = 0 # numbers, minimum periplasmic debranched oligosaccharide concentration
pdbO_max = 2e6 # numbers, maximum periplasmic debranched oligosaccharide concentration

S_e_min = 0 # numbers, minimum periplasmic sugar concentration
S_e_max = 2e7 # numbers, maximum periplasmic sugar concentration
S_c_min = 0 # numbers, minimum intracellular sugar concentration
S_c_max = 0.01 * 1e-15 * Av # numbers, maximum intracellular sugar concentration: ~10 mM (BioNumbers ID: 104695) * 1 μm^3 max cytoplasm volume
C_c_min = 0 # numbers, minimum intracellular carbon concentration
C_c_max = 0.3 * 3 * 1e-15 * Av # numbers, maximum intracellular carbon concentration: ~300 mM (BioNumbers ID: 104678) * 1 μm^3 max cytoplasm volume, assuing 3 carbon per molecule
P_c_min = 0 # numbers, minimum intracellular protein concentration
P_c_max = 500 * Qcp / C_c_mw * 1e-15 * Av # 500 mg / mL volume of cytoplasm (BioNumbers ID: 115317) 
lm_min = 0 # numbers, minimum membrane lipid concentration
lm_max = 2*(4 * pi * (1)^2) / sa_lm # numbers, maximum membrane lipid concentration: 2 monolayers of 1 μm radius composed of only lipids

# stoichiometry matrix

            #   ENDO IMPRT DEBR EXO   U    C    G    P    M    
stoich_matrix = [-0.2  0    0    0    0    0    0    0    0       # extracellular polysaccharide ... eP
                  1   -1    0    0    0    0    0    0    0       # extracellular oligosaccharide ... eO
                  0    1   -1    0    0    0    0    0    0       # periplasmic oligosaccharide ... pO
                  0    0    1  -0.2   0    0    0    0    0       # periplasmic debranched oligosaccharide pdbO
                  0    0    1    1   -1    0    0    0    0       # periplasmic sugars ... S_e
                  0    0    0    0    1   -1   -1    0    0       # intracellular substrate... S_c
                  0    0    0    0    0    4    0   -1   -1       # cellular carbon... C_c
                  0   -0.25 0    0   -2    7.2 30  -13   -4       # energy in ATP ... ATP (based on Su and refs therein, Rees et al. 2010 10.1038/nrm2646)
                  0    0    0    0    0    0    0    1    0       # carbon atoms in proteins ... P_c
                  0    0    0    0    0    0    0    0    0.02 ]  # membrane lipids ... lm





# reaction rates CAZymes
# extracellular endolaminarase GH16
GH16_kcat = 30  # 1/s
GH16_km = 0.0008/30 # M
GH16_mw = 30000 # g/mol

# TBDT in outer membrane
TBDT_kcat = 0.1 # 1/s
TBDT_km = 0.0008  # M, Balhesteros et al. 2017 https://doi.org/10.1128/jb.00723-16
TBDT_mw = 190000 # g/mol, Silale and Van den Berg 2023 10.1146/annurev-micro-032421-111116 and Sverzhinsky et al. 2015 https://doi.org/10.1128%2FJB.00069-15 (Chimento et al. 2003 https://doi.org/10.1038/nsb914)

# periplasmic debranching laminarase GH30
GH30_kcat = 50  # 1/s
GH30_km = 0.05 # M
GH30_mw = 55000 # g/mol

# periplasmic exolaminarase GH17
GH17_kcat = 80 # 1/s
GH17_km = 0.003 # M
GH17_mw = 90000 # g/mol

# reaction rates core enyzmes
# uptake
UT_kcat = 5  # 1/s, questionable (from Garcia-Alles et al. 2002 doi.org/10.1021/bi025928d), updated from Walmsey et al. 1998 10.1016/S0968-0004(98)01326-7
UT_km = 0.0001  # M, Garcia-Alles et al. 2002 doi.org/10.1021/bi025928d, updated from Walmsey et al. 1998 10.1016/S0968-0004(98)01326-7
UT_mw = 176000 # g/mol, adapted from Norris et al (dismissed https://biocyc.org/ECOLI/NEW-IMAGE?type=NIL&object=CPLX-157)

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


