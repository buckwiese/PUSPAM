# Parameters for PAM of a cell growing on a linear glucan with only 1 linkage

using JLD

# cell parameters and limits
# membrane parameters
sa_UT = 1.26e-5    # μm2, specific surface area of membrane transporter molecules
sa_lm = 0.5e-6     # μm2, specific surface area of membrane lipid molecules

lm_mw = 220        # g/mol, molecular weight of membrane lipid molecules
Mmax = 0.00045      # unitless, maximum membrane transporter to lipid ratio (counts/counts)
Mmin = 3e-6         # unitless, minimum membrane transporter to lipid ratio (counts/counts)

Av = 6.0221409e+23; #Avogadro's number (numbers/mole)
Dmax = 18e9; # Da/μm3 (BioNumbers ID: 109049), maximum cell density 
Dmin = 1e9; # fg/μm3, made up minimum cell density

# stoichiometry matrix

         # U  C  P  M
matrix = [ 1 -1  0  0   # intracellular substrate... S_c
           0  6 -1 -1   # cellular carbon... C_c
           0  0  1  0   # protein ... P_c
           0  0  0  1   # cell membrane ... lm
          -1  0  0  0 ] # extracellular substrate... S_e

S_c_mw = 180.156 # g/mol, molecular weight of intracellular substrate
C_c_mw = 12.0107 # g/mol, molecular weight of cellular carbon

# reaction rates
# uptake
UT_kcat = 10000.0  # 1/s
UT_km = 100.0  # μM
UT_mw = 20000 # g/mol

# catabolism
C_kcat = 10000.0  # 1/s
C_km = 5000.0 # μM
C_mw = 100000 # g/mol

# protein synthesis
P_kcat = 20000.0  # 1/s
P_km = 5000.0 # μM
P_mw = 80000 # g/mol

# membrane synthesis
M_kcat = 10000.0  # 1/s
M_km = 5000.0 # μM
M_mw = 50000 # g/mol


