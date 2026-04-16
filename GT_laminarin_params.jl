# Parameters for a cell growing on laminarin
using JLD

# diffusion coefficient of substrate polysaccharide laminarin (DC_eP), laminarin oligosaccharides (DC_O), and glucose (DC_S)
DC_eP = 150 # �um^2/s, based on Elyakova et al. 1994
DC_O = 300 # μum^2/s, approximated using Stokes-Einstein equation for oligosaccharides of ~6 glucose monomers
DC_S = 600 # μum^2/s, BioNumbers (ID: 104089)

# molecular weights of laminarin and its breakdown products
eP_mw = 30 * 180.156 - 29 * 18.015 # g/mol, molecular weight of extracellular polysaccharide with 30 glucose monomers
pO_mw = 6 * 180.156 - 5 * 18.015 # g/mol, molecular weight of periplasmic oligosaccharide with 6 glucose monomers
pdbO_mw = 5 * 180.156 - 4 * 18.015 # g/mol, molecular weight of periplasmic debranched oligosaccharide with 5 glucose monomers

# minima and maxima for laminarin and its breakdown products: upper limits should be maximum density * max volume of compartment / molecular weight * Avogadro's number
eO_min = 0 # numbers, minimum extracellular oligosaccharide concentration
eO_max = Dmax * 4/3 * pi * ((r_cy_max + r_pp_max + 0.005)^3 - (r_cy_max + r_pp_max)^3) * 1e-15 / pO_mw * Av # numbers, maximum extracellular oligosaccharide concentration
pO_min = 0 # numbers, minimum periplasmic oligosaccharide concentration
pO_max = Dmax * 4/3*pi*((r_cy_max + r_pp_max)^3 - (r_cy_max)^3) * 1e-15 / pO_mw * Av # numbers, maximum periplasmic oligosaccharide concentration
pdbO_min = 0 # numbers, minimum periplasmic debranched oligosaccharide concentration
pdbO_max = Dmax * 4/3*pi*((r_cy_max + r_pp_max)^3 - (r_cy_max)^3) * 1e-15 / pdbO_mw * Av # numbers, maximum periplasmic debranched oligosaccharide concentration

# stoichiometry matrix
                    # ENDO IMPRT DEBR EXO   U    C    G    P    M    
lam_stoich_matrix =   [-0.2  0    0    0    0    0    0    0    0       # extracellular polysaccharide ... eP
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
# KM given as M of monomeric sugar units
# extracellular endolaminarase GH16
GH16_kcat = 35  # 1/s
GH16_KM = 0.005/30 # M
GH16_mw = 50000 # g/mol

# TBDT in outer membrane
TBDT_kcat = 5 # 1/s, few examples in literature, see Silale and van den Berg 2023 https://doi.org/10.1146/annurev-micro-032421-111116
TBDT_KM = 0.00008/6  # M, Balhesteros et al. 2017 https://doi.org/10.1128/jb.00723-16
TBDT_mw = 190000 # g/mol, Silale and Van den Berg 2023 10.1146/annurev-micro-032421-111116 and Sverzhinsky et al. 2015 https://doi.org/10.1128%2FJB.00069-15 (Chimento et al. 2003 https://doi.org/10.1038/nsb914)

# periplasmic debranching laminarase GH30
GH30_kcat = 66  # 1/s
GH30_KM = 0.560/6 # M
GH30_mw = 55000 # g/mol

# periplasmic exolaminarase GH17
GH17_kcat = 80 # 1/s
GH17_KM = 0.015/5 # M
GH17_mw = 90000 # g/mol




