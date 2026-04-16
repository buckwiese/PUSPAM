# Constraints for cell growing on laminarin
using JuMP, Ipopt, Random, Distributions

# Set up model
model = Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "max_iter", 500000)

# Initialize extracellular macromolecular substrate concentration parameter
@NLparameter(model, eP_env == 1) # M,  extracellular substrate concentration

# Variables to optimize
@variables(model, begin
    logmu, (start = rand(Uniform(-2,-1))) # 1/s, growth rate

    eP_surf_rel >= 0, (start = rand()) # fraction, substrate polysaccharide at cell surface
    eO_rel >= 0, (start = rand())      # fraction, extracellular oligosaccharide 
    pO_rel >= 0, (start = rand())      # fraction, periplasmic oligosaccharide 
    pdbO_rel >= 0, (start = rand())    # fraction, periplasmic debranched oligosaccharide

    S_e_rel >= 0, (start = rand())     # fraction, periplasmic sugar
    S_c_rel >= 0, (start = rand())     # fraction, intracellular substrate
    C_c_rel >= 0, (start = rand())     # fraction, cellular carbon 
    P_c_rel >= 0, (start = rand())     # fraction, cellular protein 
    lm_rel >= 0, (start = rand())      # fraction, cellular membrane
    
    phi_GH16 >=0.0, (start = rand())    # fraction, proteome allocation to extracellular endolaminarase GH16
    phi_TBDT >=0.0, (start = rand())    # fraction, proteome allocation to TBDT in outer membrane
    phi_GH30 >=0.0, (start = rand())    # fraction, proteome allocation to periplasmic debranching laminarase GH30
    phi_GH17 >=0.0, (start = rand())    # fraction, proteome allocation to periplasmic exolaminarase GH17

    phi_C >= 0.0, (start = rand())      # fraction, proteome allocation to catabolism
    phi_P >= 0.0, (start = rand())      # fraction, proteome allocation to protein synthesis
    phi_M >= 0.0, (start = rand())      # fraction, proteome allocation to membrane synthesis
    phi_UT >= 0.0, (start = rand())     # fraction, proteome allocation to uptake

    lm_inner >= 0.1, (start = rand())   # fraction of lipids in inner membrane
    C_to_ATP >= 0, (start = rand())     # fraction, carbon to ATP investment ratio
end)
    
# Set the objective: maximize the log growth rate.
@objective(model, Max, logmu)

# Define intermediate calculations.
@NLexpressions(model, begin
    # metabolite numbers
    eP_surf, eP_surf_rel * eP_env   # number, substrate polysaccharide at cell surface
    eO, eO_rel * eO_max             # number, extracellular oligosaccharide
    pO, pO_rel * pO_max             # number, periplasmic oligosaccharide
    pdbO, pdbO_rel * pdbO_max       # number, periplasmic debranched oligosaccharide

    S_e, S_e_rel * S_e_max          # number, periplasmic sugar
    S_c, S_c_rel * S_c_max          # number, intracellular substrate
    C_c, C_c_rel * C_c_max          # number, cellular carbon
    P_c, P_c_rel * P_c_max          # number, cellular protein
    lm, lm_rel * lm_max             # number, cellular membrane

    # protein numbers
    GH16, phi_GH16 * P_c * C_c_mw/GH16_mw/Qcp   # number, extracellular endolaminarase GH16
    TBDT, phi_TBDT * P_c * C_c_mw/TBDT_mw/Qcp   # number, TBDT in outer membrane
    GH30, phi_GH30 * P_c * C_c_mw/GH30_mw/Qcp   # number, periplasmic debranching laminarase GH30
    GH17, phi_GH17 * P_c * C_c_mw/GH17_mw/Qcp   # number, periplasmic exolaminarase GH17

    UT, phi_UT * P_c * C_c_mw/UT_mw/Qcp         # number, uptake transporters
    CE, phi_C * P_c * C_c_mw/C_mw/Qcp           # number, catabolism enzymes
    PE, phi_P * P_c * C_c_mw/P_mw/Qcp           # number, protein synthesis enzymes
    ME, phi_M * P_c * C_c_mw/M_mw/Qcp           # number, membrane synthesis enzymes
    
    # cell parameters
    area_cy, sa_lm*lm/2*lm_inner + sa_UT*UT                             # μm2, cytoplasm surface area
    r_cy, sqrt(area_cy / (4*pi))                                        # μm, cytoplasm radius
    volume_cy,  4/3*pi*r_cy^3                                           # μm3, cytoplasm volume
    area_pp, sa_lm*lm/2*(1-lm_inner) + TBDT*sa_TBDT                     # μm2, periplasm surface area
    r_pp, sqrt(area_pp/(4*pi)) -r_cy                                    # μm, periplasm thickness
    volume_pp, 4/3*pi*(r_pp+r_cy)^3 - volume_cy                         # μm3, periplasm volume
    volume_excelrp, 4/3*pi*(r_pp+r_cy+delta)^3 - volume_cy - volume_pp  # μm3, extracellular reaction space volume
#    volume_btwcells, V_per_cell - volume_cy - volume_pp                 # μm3, volume between cells

    # metabolite concentrations
    eO_M, eO/(volume_excelrp*1e-15*Av)     # M, extracellular oligosaccharide concentration
    pO_M, pO/(volume_pp*1e-15*Av)          # M, periplasmic oligosaccharide concentration
    pdbO_M, pdbO/(volume_pp*1e-15*Av)      # M, periplasmic debranched oligosaccharide concentration
    S_e_M, S_e/(volume_pp*1e-15*Av)        # M, periplasmic sugar concentration
    S_c_M, S_c/(volume_cy*1e-15*Av)        # M, intracellular substrate concentration
    C_c_M, C_c/(volume_cy*1e-15*Av)        # M, cellular carbon concentration
    P_c_M, P_c/(volume_cy*1e-15*Av)        # M, cellular protein concentration
    lm_M, lm/(volume_cy*1e-15*Av)          # M, cellular membrane concentration
    
    # CAZyme concentrations and distributions
    GH16_M, GH16/(volume_excelrp*1e-15*Av)      # M, endolaminarase GH16 concentration in extracellular reaction space
    TBDT_M, TBDT/(volume_excelrp*1e-15*Av)      # M, TBDT concentration in outer membrane
    GH30_M, GH30/(volume_pp*1e-15*Av)           # M, debranching laminarase GH30 concentration in periplasm
    GH17_M, GH17/(volume_pp*1e-15*Av)           # M, exolaminarase GH17 concentration in periplasm
    
    # central enzyme concentrations and distributions
    UT_M, UT/(volume_pp*1e-15*Av)       # M, uptake transporter concentration in periplasm
    CE_M, CE/(volume_cy*1e-15*Av)       # M, catabolism enzyme concentration in cytoplasm
    PE_M, PE/(volume_cy*1e-15*Av)       # M, protein synthesis enzyme concentration in cytoplasm
    ME_M, ME/(volume_cy*1e-15*Av)       # M, membrane synthesis enzyme concentration in cytoplasm
    
    # reaction rates
    v_ENDO, GH16_kcat * GH16 * eP_surf / (GH16_KM + eP_surf)      # numbers/s, extracellular hydrolysis
    v_IMPRT, TBDT_kcat * TBDT * eO_M / (TBDT_KM + eO_M)           # numbers/s, uptake
    v_DEBR, GH30_kcat * GH30 * pO_M / (GH30_KM + pO_M)            # numbers/s, periplasmic debranching
    v_EXO, GH17_kcat * GH17 * pdbO_M / (GH17_KM + pdbO_M)         # numbers/s, periplasmic hydrolysis

    v_UT, UT_kcat * UT * S_e_M / (UT_KM + S_e_M)                  # numbers/s, uptake
    v_C, C_to_ATP * C_kcat * CE * S_c_M / (C_KM + S_c_M)          # numbers/s, catabolism
    v_G, (1 - C_to_ATP) * C_kcat * CE * S_c_M / (C_KM + S_c_M)    # numbers/s, ATP production
    v_P, P_kcat * PE * C_c_M / (P_KM + C_c_M)                     # numbers/s, protein synthesis
    v_M, M_kcat * ME * C_c_M / (M_KM + C_c_M)                     # numbers/s, membrane synthesis
    
    # cell density
    D_cy, UT_M*UT_mw + lm_M*lm_inner*lm_mw + CE_M*C_mw + PE_M*P_mw + ME_M*M_mw + C_c_M*C_c_mw*3 + S_c_M*S_mw  # g/L, cell density, assuming c makes up 1/3 of cellular carbon metabolite
    D_pp, TBDT_M*TBDT_mw + lm_M*(1-lm_inner)*lm_mw + GH30_M*GH30_mw + GH17_M*GH17_mw + pO_M*pO_mw + pdbO_M*pdbO_mw  # g/L, cell density in periplasm
    
    # diffusion processes
    eP_diff, DC_eP/(r_cy + r_pp) * (eP_env - eP_surf)*1e-15*Av  # numbers/um2 s, diffusion rate of substrate polysaccharide (eP) to cell surface per unit area
    eO_diff, DC_O/(r_cy + r_pp) * eO_M*1e-15*Av                   # numbers/um2 s, diffusion loss rate of extracellular oligosaccharides (eO) away from cell surface per unit area

    # production rates
    ATP_rate, lam_stoich_matrix[8, 2]*v_IMPRT + lam_stoich_matrix[8, 5]*v_UT + lam_stoich_matrix[8, 6]*v_C + lam_stoich_matrix[8, 7]*v_G + lam_stoich_matrix[8, 8]*v_P + lam_stoich_matrix[8, 9]*v_M # number/s, ATP production rate
    eO_rate, lam_stoich_matrix[2, 1]*v_ENDO + lam_stoich_matrix[2, 2]*v_IMPRT # number/s, extracellular oligosaccharide production rate
    pO_rate, lam_stoich_matrix[3, 2]*v_IMPRT + lam_stoich_matrix[3, 3]*v_DEBR # number/s, periplasmic oligosaccharide production rate
    pdbO_rate, lam_stoich_matrix[4, 3]*v_DEBR + lam_stoich_matrix[4, 4]*v_EXO # number/s, periplasmic debranched oligosaccharide production rate
    S_e_rate, lam_stoich_matrix[5, 3]*v_DEBR + lam_stoich_matrix[5, 4]*v_EXO + lam_stoich_matrix[5, 5]*v_UT # number/s, periplasmic sugar production rate
    S_c_rate, lam_stoich_matrix[6, 5]*v_UT + lam_stoich_matrix[6, 6]*v_C + lam_stoich_matrix[6, 7]*v_G # number/s, intracellular substrate production rate
    C_c_rate, lam_stoich_matrix[7, 6]*v_C + lam_stoich_matrix[7, 8]*v_P + lam_stoich_matrix[7, 9]*v_M # number/s, cellular carbon production rate
    P_c_rate, lam_stoich_matrix[9, 8]*v_P # number/s, cellular protein production rate
    lm_rate, lam_stoich_matrix[10, 9]*v_M # number/s, cellular membrane production rate

    # rates per area
    eP_degr, v_ENDO /area_pp # numbers/um2 s, extracellular polysaccharide degradation rate per unit area
    eO_prod, eO_rate/area_pp   # numbers/um2 s, net extracellular oligosaccharide production rate per unit area
end)
    
# Define the constraints.
# Linear constraints
@constraints(model, begin
    # proteome allocation constraint
    phi_C + phi_P + phi_M + phi_UT + phi_GH16 + phi_TBDT + phi_GH30 + phi_GH17 == 0.5 # proteome allocation constraint, reserves 50% of proteome for housekeeping functions
end)

# Nonlinear constraints
@NLconstraints(model, begin
    # diffusion constraints
    eP_diff - eP_degr >= 0 # mol/um2 s, diffusion supply constraint for substrate polysaccharides
    eO_prod >= 0
    # constraints on the cell membrane architechture and composition
    UT / (lm/2*lm_inner) <= Mmax            # upper boundary of transporter to lipid ratio in IM
    (TBDT)/ (lm/2*(1-lm_inner)) <= Mmax     # upper boundary of transporter to lipid ratio in OM
    sa_CAZY*GH16/(sa_lm*lm/2*(1-lm_inner)) <= 1  # upper boundary of enzyme on cell surface
    lm_inner + (1-lm_inner) == 1.0          # lipid partitioning between IM and OM
    # cell stoichiometry constraints
    eO_rate - exp(logmu) * eO == 0   # number/s, steady-state balanced, growth-rate constraint for extracellular oligosaccharide
    pO_rate - exp(logmu) * pO == 0                   # number/s, steady-state balanced, growth-rate constraint for periplasmic oligosaccharide
    pdbO_rate - exp(logmu) * pdbO == 0               # number/s, steady-state balanced, growth-rate constraint for periplasmic debranched oligosaccharide
    S_e_rate - exp(logmu) * S_e == 0                 # number/s, steady-state balanced, growth-rate constraint for periplasmic sugars
    S_c_rate - exp(logmu) * S_c == 0                 # number/s, steady-state balanced, growth-rate constraint for intracellular substrate
    C_c_rate - exp(logmu) * C_c == 0                 # number/s, steady-state balanced, growth-rate constraint for cellular carbon
    ATP_rate >= 0                                    # constraint for energy in ATP
    P_c_rate - exp(logmu) * P_c == 0                 # number/s, steady-state balanced, growth-rate constraint for cellular protein
    lm_rate - exp(logmu) * lm == 0                   # number/s, steady-state balanced, growth-rate constraint for cellular membrane
    # proteome allocation constraints
    phi_GH16 * v_P * C_c_mw/GH16_mw/Qcp - exp(logmu) * GH16 == 0     # proteome allocation constraint for extracellular endolaminarase GH16
    phi_TBDT * v_P * C_c_mw/TBDT_mw/Qcp - exp(logmu) * TBDT == 0     # proteome allocation constraint for TBDT in outer membrane
    phi_GH30 * v_P * C_c_mw/GH30_mw/Qcp - exp(logmu) * GH30 == 0     # proteome allocation constraint for periplasmic debranching laminarase GH30
    phi_GH17 * v_P * C_c_mw/GH17_mw/Qcp - exp(logmu) * GH17 == 0     # proteome allocation constraint for periplasmic exolaminarase GH17
    phi_UT * v_P * C_c_mw/UT_mw/Qcp - exp(logmu) * UT == 0   # proteome allocation constraint for uptake    
    phi_C * v_P * C_c_mw/C_mw/Qcp - exp(logmu) * CE == 0     # proteome allocation constraint for catabolism
    phi_P * v_P * C_c_mw/P_mw/Qcp - exp(logmu) * PE == 0     # proteome allocation constraint for protein synthesis
    phi_M * v_P * C_c_mw/M_mw/Qcp - exp(logmu) * ME == 0     # proteome allocation constraint for membrane synthesis
    # upper and lower boundaries of physical parameters
    D_cy <= Dmax                # maximum density constraint; g/L
    D_cy >= Dmin                # minimum density constraint; g/L
    D_pp <= Dmax                # maximum density constraint; g/L
    D_pp >= Dmin                # minimum density constraint; g/L
    r_cy >= r_cy_min            # minimum cell cytoplasm radius; μm
    r_cy <= r_cy_max            # maximum cell cytoplasm radius; μm
    r_pp >= r_pp_min            # minimum cell periplasm radius; μm
    r_pp <= r_pp_max            # maximum cell periplasm radius; μm
    # boundaries of metabolite pools
    C_to_ATP >= 0           # minimum carbon to ATP investment ratio
    C_to_ATP <= 1           # maximum carbon to ATP investment ratio; fraction of C_to_ATP
    eP_surf >= 0            # minimum substrate polysaccharide concentration at cell surface
    eP_surf_rel <= 1.0      # maximum substrate polysaccharide concentration at cell surface; fraction of eP_env
    eO_rel >= 0.0           # minimum extracellular oligosaccharide concentration
    eO_rel <= 1.0           # maximum extracellular oligosaccharide concentration; fraction of eO_max
    pO_rel <= 1.0           # maximum periplasmic oligosaccharide concentration; fraction of pO_max
    pdbO_rel <= 1.0         # maximum periplasmic debranched oligosaccharide concentration; fraction of pdbO_max
    S_e_rel <= 1.0          # maximum periplasmic sugar concentration; fraction of S_e_max
    S_c_rel <= 1.0          # maximum intracellular substrate concentration; fraction of S_c_max
    C_c_rel <= 1.0          # maximum cellular carbon concentration; fraction of C_c_max
    P_c_rel <= 1.0          # maximum cellular protein concentration; fraction of P_c_max
    lm_rel <= 1.0           # maximum cellular membrane concentration; fraction of lm_max
end)

