# Global cell parameters
using JLD

# Constants
Av = 6.0221409e+23; #Avogadro's number (numbers/mole)
# average carbon content of proteins
Qcp = 0.5; # g carbon/g protein
# Cells per liter
cell_count_L = 1e9; # cells/L, typical cell concentration in marine systems (10^6 cells/mL)

# Cell parameters and limits
# cell density and size limits
Dmax = 500; # g/L to accommodate (BioNumbers ID: 115317, Oldewurtel et al 2023 10.1073/pnas.2021416118), maximum cell density 
Dmin = 250; # g/L, (Oldewurtel et al 2023 10.1073/pnas.2021416118) minimum cell density
r_cy_min = 0.1 # μm, minimum cytoplasm radius, given filter poresizes
r_cy_max = 1.0 # μm, maximum cytoplasm radius 
r_pp_min = 0.01 # μm, minimum periplasm radius (BioNumbers ID: 114122)
r_pp_max = 0.025 # μm, maximum periplasm radius (BioNumbers ID: 114122)
# extracellular layer thickness
delta = 0.005 # μm, approximate diameter of extracellular enzymes (may need supp figure?)

# membrane parameters
sa_TBDT = 0.01*0.005 # μm2, specific surface area of TBDT molecules (PDB 8AA2)
sa_CAZY = 0.005*0.005  # μm2, specific surface area occupied by extracellular enzymes (conservative)
sa_UT = 0.005*0.01    # μm2, specific surface area of membrane transporter molecules (BioNumbers ID: 114685)
sa_lm = 0.4e-6     # μm2, specific surface area of membrane lipid molecules (BioNumbers ID: 114186)

lm_mw = 860         # g/mol, molecular weight of membrane lipid molecules
Mmax = 0.0045        # unitless, maximum membrane transporter to lipid ratio (counts/counts)
Mmin = 3e-6         # unitless, minimum membrane transporter to lipid ratio (counts/counts)

# molecular weights of metabolites
S_mw = 180.156 # g/mol, molecular weight of glucose
C_c_mw = 12.0107 # g/mol, molecular weight of cellular carbon
ATP_mw = 507.181 # g/mol, molecular weight of ATP

# limits of compound pools: upper limits should be maximum density * max volume of compartment / molecular weight * Avogadro's number
cy_vol_max = (4/3)*pi*(r_cy_max)^3 # μm^3, maximum cytoplasm volume
pp_vol_max = (4/3)*pi*((r_cy_max + r_pp_max)^3 - (r_cy_max)^3) # μm^3, maximum periplasm volume
S_e_min = 0 # numbers, minimum periplasmic sugar concentration
S_e_max = Dmax * pp_vol_max * 1e-15 / S_mw * Av # numbers, maximum periplasmic sugar concentration
S_c_min = 0 # numbers, minimum intracellular sugar concentration
S_c_max = Dmax * cy_vol_max * 1e-15 / S_mw * Av # numbers, maximum intracellular sugar concentration (a lot higher than ~10 mM (BioNumbers ID: 104695)!)
C_c_min = 0 # numbers, minimum intracellular carbon concentration
C_c_max = Dmax * cy_vol_max * 1e-15 / C_c_mw * Av # numbers, maximum intracellular carbon concentration (one order of mag. higher than ~300 mM (BioNumbers ID: 104678))
P_c_min = 0 # numbers, minimum intracellular protein concentration
P_c_max = Dmax * cy_vol_max * 1e-15 / (C_c_mw*Qcp) * Av # 500 mg / mL volume of cytoplasm (BioNumbers ID: 115317) 
lm_min = 0 # numbers, minimum membrane lipid concentration
lm_max = ((2 * 4*pi*(r_cy_max)^2) + (2 * 4*pi*(r_cy_max + r_pp_max)^2)) / sa_lm # numbers, maximum membrane lipid concentration: bilayers of 0.6 μm radius and 0.625 μm radius composed of only lipids

# uptake
UT_kcat = 10  # 1/s, questionable (from Garcia-Alles et al. 2002 doi.org/10.1021/bi025928d), updated from Walmsey et al. 1998 10.1016/S0968-0004(98)01326-7
UT_KM = 0.0001  # M, Garcia-Alles et al. 2002 doi.org/10.1021/bi025928d, updated from Walmsey et al. 1998 10.1016/S0968-0004(98)01326-7
UT_mw = 176000 # g/mol, adapted from Norris et al (dismissed https://biocyc.org/ECOLI/NEW-IMAGE?type=NIL&object=CPLX-157)

# catabolism
C_kcat = 50.0  # 1/s, Cerisy et al. 2019 https://doi.org/10.1128/jb.00241-19 and Norris et al
C_KM = 0.001 # M, ?! adapted from Norris et al
C_mw = 4400000 # g/mol, adapted from Norris et al

# protein synthesis
P_kcat = 25.0  # 1/s, adapted from Norris et al
P_KM = 0.003 # M, ?! adapted from Norris et al
P_mw = 2200000 # g/mol

# membrane synthesis
M_kcat = 10.0  # 1/s, adapted from Norris et al 
M_KM = 0.002 # M, ?! adapted from Norris et al
M_mw = 2310000 # g/mol, adapted from Norris et al
