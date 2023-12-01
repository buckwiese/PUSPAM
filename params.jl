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

# stoichiometry matrix

         # U  C  P  M
matrix = [ 1 -1  0  0   # intracellular substrate... S_c
           0  6 -1 -50   # cellular carbon... C_c
           0  0  1  0   # protein ... P_c
           0  0  0  1   # cell membrane ... lm
          -1  0  0  0 ] # extracellular substrate... S_e

S_c_mw = 180.156 # g/mol, molecular weight of intracellular substrate
C_c_mw = 12.0107 # g/mol, molecular weight of cellular carbon

# reaction rates
# uptake
UT_kcat = 100.0  # 1/s
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


