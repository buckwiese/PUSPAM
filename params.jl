# Parameters for PAM of a cell growing on a linear glucan with only 1 linkage

using JLD

# cell parameters and limits
# membrane parameters
sa_UT = 1.26e-5    # μm2, specific surface area of membrane transporter molecules
sa_lm = 0.4e-6     # μm2, specific surface area of membrane lipid molecules

lm_mw = 860         # g/mol, molecular weight of membrane lipid molecules
Mmax = 0.00045        # unitless, maximum membrane transporter to lipid ratio (counts/counts)
Mmin = 3e-6         # unitless, minimum membrane transporter to lipid ratio (counts/counts)

Av = 6.0221409e+23; #Avogadro's number (numbers/mole)
Dmax = 340; # g/L (BioNumbers ID: 109049), maximum cell density 
Dmin = 300; # g/L, made up minimum cell density

# average carbon content of proteins
Qcp = 0.5; # g carbon/g protein 


# sugar matrix

            #    ENDO IMPRT DEBR EXO 
sugar_matrix = [-0.2    0    0    0  # extracellular polysaccharide ... eP
                  1    -1    0    0  # extracellular oligosaccharide ... eO
                  0     1   -1    0  # periplasmic oligosaccharide ... pO
                  0     0    1  -0.2 # periplasmic debranched oligosaccharide pdbO
                  0     0    1    1] # periplasmic sugars ... S_e

# stoichiometry matrix

         # U  C  P  M
matrix = [ 1 -1  0  0   # intracellular substrate... S_c
           0  6 -1 -50   # cellular carbon... C_c
           0  0  1  0   # protein ... P_c
           0  0  0  1   # cell membrane ... lm
          -1  0  0  0 ] # periplasmic sugars... S_e

S_c_mw = 180.156 # g/mol, molecular weight of intracellular substrate
C_c_mw = 12.0107 # g/mol, molecular weight of cellular carbon
pO_mw = 6*180.156-5*18 # g/mol, molecular weight of periplasmic oligosaccharide
pdbO_mw = 5*180.156-4*18 # g/mol, molecular weight of periplasmic debranched oligosaccharide

# reaction rates CAZymes
# extracellular endolaminarase GH16
GH16_kcat = 10  # 1/s
GH16_km = 0.001 # M
GH16_mw = 50000 # g/mol

# TBDT in outer membrane
TBDT_kcat = 1000.0  # 1/s
TBDT_km = 0.0001  # M
TBDT_mw = 80500 # g/mol

# periplasmic debranching laminarase GH30
GH30_kcat = 40  # 1/s
GH30_km = 0.020 # M
GH30_mw = 65000 # g/mol

# periplasmic exolaminarase GH17
GH17_kcat = 200 # 1/s
GH17_km = 0.001 # M
GH17_mw = 38000 # g/mol

# reaction rates core enyzmes
# uptake
UT_kcat = 1000.0  # 1/s
UT_km = 0.005  # M
UT_mw = 115060 # g/mol

# catabolism
C_kcat = 500.0  # 1/s
C_km = 0.001 # M
C_mw = 3348400 # g/mol

# protein synthesis
P_kcat = 50.0  # 1/s
P_km = 0.01 # M
P_mw = 2133230 # g/mol

# membrane synthesis
M_kcat = 100.0  # 1/s
M_km = 0.01 # M
M_mw = 1511620 # g/mol


