using JuMP, Ipopt

#eq for PAM optimizer
# set up model
model = Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "max_iter", 15000)

# Set the environment:
#@NLparameter(model, T   == 20.0 + 273.0)   # K,  temperature
@NLparameter(model, eP_env == 1) # M,  external substrate concentration
println("eP_env = ", value(eP_env))

#variables to optimize
#concentrations in M
@variables(model, begin
    logmu, (start = -7) # 1/s, growth rate
    eP_surf >= 0, (start = 1e-5) # M, substrate polysaccharide at cell surface
    eO >= 0, (start = 100000)   # number, extracellular oligosaccharide 
    pO >= 0, (start = 100000)   # number, periplasmic oligosaccharide 
    pdbO >= 0, (start = 100000) # number, periplasmic debranched oligosaccharide
    S_e >= 0, (start = 100000)  # number, periplasmic sugar
    S_c >= 1, (start = 100000) # number, intracellular substrate
    C_c >= 1, (start = 100000) # number, cellular carbon 
    P_c >= 1, (start = 100000)   # number, cellular protein 
    lm >= 1, (start = 100000)  # number, cellular membrane
    
    phi_porin >= 0.0, (start = 0.05) # relative proteome allocation to outer membrane porins
    phi_GH16 >=0.0, (start = 0.05) # relative proteome allocation to extracellular endolaminarase GH16
    phi_TBDT >=0.0, (start = 0.1) # relative proteome allocation to TBDT in outer membrane
    phi_GH30 >=0.0, (start = 0.01) # relative proteome allocation to periplasmic debranching laminarase GH30
    phi_GH17 >=0.0, (start = 0.01) # relative proteome allocation to periplasmic exolaminarase GH17
    phi_C >= 0.0, (start = 0.06) # relative proteome allocation to catabolism
    phi_P >= 0.0, (start = 0.06) # relative proteome allocation to protein synthesis
    phi_M >= 0.0, (start = 0.06) # relative proteome allocation to membrane synthesis
    phi_UT >= 0.0, (start = 0.03) # relative proteome allocation to uptake

    lm_inner >= 0.1, (start = 0.5) # fraction of lipids in inner membrane
    C_to_ATP >= 0, (start = 0.5) # mol/mol, carbon to ATP investment ratio
end)
    
# Set the objective: maximize the log growth rate.
@objective(model, Max, logmu)

# Define intermediate calculations.
@NLexpressions(model, begin
    # protein numbers
    porin, phi_porin * P_c * C_c_mw/porin_mw/Qcp # number, outer membrane porin channels
    GH16, phi_GH16 * P_c * C_c_mw/GH16_mw/Qcp # number, extracellular endolaminarase GH16
    TBDT, phi_TBDT * P_c * C_c_mw/TBDT_mw/Qcp # number, TBDT in outer membrane
    GH30, phi_GH30 * P_c * C_c_mw/GH30_mw/Qcp # number, periplasmic debranching laminarase GH30
    GH17, phi_GH17 * P_c * C_c_mw/GH17_mw/Qcp # number, periplasmic exolaminarase GH17
    UT, phi_UT * P_c * C_c_mw/UT_mw/Qcp # number, uptake transporters
    CE, phi_C * P_c * C_c_mw/C_mw/Qcp # number, catabolism enzymes
    PE, phi_P * P_c * C_c_mw/P_mw/Qcp # number, protein synthesis enzymes
    ME, phi_M * P_c * C_c_mw/M_mw/Qcp # number, membrane synthesis enzymes
    # cell parameters
    area_cy, sa_lm*lm*lm_inner + sa_UT*UT           # μm2, cytoplasm surface area
    r_cy, sqrt(area_cy / (4*pi))                    # μm, cytoplasm radius
    volume_cy,  4/3*pi*r_cy^3                       # μm3, cytoplasm volume
    area_pp, sa_lm*lm*(1-lm_inner) + TBDT*sa_TBDT + porin*sa_porin       # μm2, periplasm surface area
    r_pp, sqrt(area_pp/(4*pi)) -r_cy                      # μm, periplasm thickness
    volume_pp, 4/3*pi*(r_pp+r_cy)^3 - volume_cy     # μm3, periplasm volume
    volume_excelrp, 4/3*pi*(r_pp+r_cy+delta)^3 - volume_cy - volume_pp   # μm3, extracellular reaction space volume
    # metabolite concentrations
    eO_M, eO/(volume_excelrp*1e-15*Av)     # M, extracellular oligosaccharide concentration
    pO_M, pO/(volume_pp*1e-15*Av)          # M, periplasmic oligosaccharide concentration
    pdbO_M, pdbO/(volume_pp*1e-15*Av)      # M, periplasmic debranched oligosaccharide concentration
    S_e_M, S_e/(volume_pp*1e-15*Av)        # M, periplasmic sugar concentration
    S_c_M, S_c/(volume_cy*1e-15*Av)        # M, intracellular substrate concentration
    C_c_M, C_c/(volume_cy*1e-15*Av)        # M, cellular carbon concentration
    P_c_M, P_c/(volume_cy*1e-15*Av)        # M, cellular protein concentration
    lm_M, lm/(volume_cy*1e-15*Av)          # M, cellular membrane concentration
    # CAZyme numbers and distributions
    porin_M, porin/(volume_excelrp*1e-15*Av)     # M, porin concentration in outer membrane
    GH16_M, GH16/(volume_excelrp*1e-15*Av)       # M, endolaminarase GH16 concentration in extracellular reaction space
    TBDT_M, TBDT/(volume_excelrp*1e-15*Av)       # M, TBDT concentration in outer membrane
    GH30_M, GH30/(volume_pp*1e-15*Av)       # M, debranching laminarase GH30 concentration in periplasm
    GH17_M, GH17/(volume_pp*1e-15*Av)       # M, exolaminarase GH17 concentration in periplasm
    UT_M, UT/(volume_pp*1e-15*Av)       # M, uptake transporter concentration in periplasm
    CE_M, CE/(volume_cy*1e-15*Av)       # M, catabolism enzyme concentration in cytoplasm
    PE_M, PE/(volume_cy*1e-15*Av)       # M, protein synthesis enzyme concentration in cytoplasm
    ME_M, ME/(volume_cy*1e-15*Av)       # M, membrane synthesis enzyme concentration in cytoplasm
    # reaction rates
    v_ENDO, GH16_kcat * GH16_M * eP_surf / (GH16_km + eP_surf)  # M/s, extracellular hydrolysis
    v_IMPRT, TBDT_kcat * TBDT_M * eO_M / (TBDT_km + eO_M)       # M/s, uptake
    v_DEBR, GH30_kcat * GH30_M * pO_M / (GH30_km + pO_M)        # M/s, periplasmic debranching
    v_EXO, GH17_kcat * GH17_M * pdbO_M / (GH17_km + pdbO_M)     # M/s, periplasmic hydrolysis
    v_UT, UT_kcat * UT_M * S_e_M / (UT_km + S_e_M)              # M/s, uptake
    v_C, C_to_ATP * C_kcat * CE_M * S_c_M / (C_km + S_c_M)      # M/s, catabolism
    v_G, (1 - C_to_ATP) * C_kcat * CE_M * S_c_M / (C_km + S_c_M)       # M/s, ATP production
    v_P, P_kcat * PE_M * C_c_M / (P_km + C_c_M)       # M/s, protein synthesis
    v_M, M_kcat * ME_M * C_c_M / (M_km + C_c_M)       # M/s, membrane synthesis
    # cell density
    D_cy, UT_M*UT_mw + CE_M*C_mw + PE_M*P_mw + ME_M*M_mw + C_c_M*C_c_mw*3 + S_c_M*S_c_mw  # g/L, cell density, assuming c makes up 1/3 of cellular carbon metabolite
    D_pp, porin_M*porin_mw + GH16_M*GH16_mw + TBDT_M*TBDT_mw + GH30_M*GH30_mw + GH17_M*GH17_mw + pO_M*pO_mw + pdbO_M*pdbO_mw  # g/L, cell density, assuming c makes up 30% of cell mass
    # diffusion processes
    diff_rate, DC_eP * (eP_env-eP_surf)*1e-15*Av * (4*pi*(r_cy+r_pp+delta)^2) / 50 # number/s, diffusion rate of substrate polysaccharide (eP) to cell surface based on Karp-Boss et al. 1996 
    oligoloss_rate, DC_O * (eO_M)*volume_excelrp*1e-15*Av * 4*pi*(r_cy+r_pp+delta)^2 / 50 # number/s, diffusion rate of oligosaccharide (eO) from cell surface based on Karp-Boss et al. 1996
#    porin_O_diff, DC_O*porin_pore*1e-15*Av*volume_excelrp*porin_M*(eO_M-pO_M)/(0.008) # mol/s, diffusion rate of extracellular oligosaccharides (eO) through porins into the periplasm
#    porin_dbO_diff, DC_O*porin_pore*1e-15*Av*volume_excelrp*porin_M*(pdbO_M)/(0.008) # mol/s, diffusion rate of extracellular oligosaccharides (eO) through porins into the periplasm
    # production rates
    ATP_rate, (stoich_matrix[8, 2]*v_IMPRT*volume_excelrp + stoich_matrix[8, 5]*v_UT*volume_pp + stoich_matrix[8, 6]*v_C*volume_cy + stoich_matrix[8, 7]*v_G*volume_cy + stoich_matrix[8, 8]*v_P*volume_cy + stoich_matrix[8, 9]*v_M*volume_cy)*1e-15*Av # number/s, ATP production rate
    eO_rate, (stoich_matrix[2, 1]*v_ENDO + stoich_matrix[2, 2]*v_IMPRT)*volume_excelrp*1e-15*Av # number/s, extracellular oligosaccharide production rate
    pO_rate, (stoich_matrix[3, 2]*v_IMPRT + stoich_matrix[3, 3]*v_DEBR)*volume_pp*1e-15*Av # number/s, periplasmic oligosaccharide production rate
    pdbO_rate, (stoich_matrix[4, 3]*v_DEBR + stoich_matrix[4, 4]*v_EXO)*volume_pp*1e-15*Av # number/s, periplasmic debranched oligosaccharide production rate
    S_e_rate, (stoich_matrix[5, 3]*v_DEBR + stoich_matrix[5, 4]*v_EXO + stoich_matrix[5, 5]*v_UT)*volume_pp*1e-15*Av # number/s, periplasmic sugar production rate
    S_c_rate, (stoich_matrix[6, 5]*v_UT + stoich_matrix[6, 6]*v_C + stoich_matrix[6, 7]*v_G)*volume_cy*1e-15*Av # number/s, intracellular substrate production rate
    C_c_rate, (stoich_matrix[7, 6]*v_C + stoich_matrix[7, 8]*v_P + stoich_matrix[7, 9]*v_M)*volume_cy*1e-15*Av # number/s, cellular carbon production rate
    P_c_rate, (stoich_matrix[9, 8]*v_P)*volume_cy*1e-15*Av # number/s, cellular protein production rate
    lm_rate, (stoich_matrix[10, 9]*v_M)*volume_cy*1e-15*Av # number/s, cellular membrane production rate
end)
    
# Define the constraints.
# Linear constraints
@constraints(model, begin
    # proteome allocation constraint
    phi_C + phi_P + phi_M + phi_UT + phi_porin + phi_GH16 + phi_TBDT + phi_GH30 + phi_GH17 == 0.5 # proteome allocation constraint
end)

# Nonlinear constraints
@NLconstraints(model, begin
    # constraints on the cell membrane architechture and composition
    UT / (lm*lm_inner) <= Mmax                                          # upper boundary of transporter to lipid ratio in IM
    (TBDT + porin)/ (lm*(1-lm_inner)) <= Mmax                           # upper boundary of transporter to lipid ratio in OM
    lm_inner + (1-lm_inner) == 1.0                                      # lipid partitioning between IM and OM
    # cell stoichiometry constraints
    eO_rate - exp(logmu) * eO - oligoloss_rate == 0  # number/s, steady-state balanced, growth-rate constraint for extracellular oligosaccharide
    pO_rate - exp(logmu) * pO == 0  # number/s, steady-state balanced, growth-rate constraint for periplasmic oligosaccharide
    pdbO_rate - exp(logmu) * pdbO == 0  # number/s, steady-state balanced, growth-rate constraint for periplasmic debranched oligosaccharide
    S_e_rate - exp(logmu) * S_e == 0    # number/s, steady-state balanced, growth-rate constraint for periplasmic sugars
    S_c_rate - exp(logmu) * S_c == 0    # number/s, steady-state balanced, growth-rate constraint for intracellular substrate
    C_c_rate - exp(logmu) * C_c == 0    # number/s, steady-state balanced, growth-rate constraint for cellular carbon
    ATP_rate >= 0     # constraint for energy in ATP
    P_c_rate - exp(logmu) * P_c == 0    # number/s, steady-state balanced, growth-rate constraint for cellular protein
    lm_rate - exp(logmu) * lm == 0      # M/s, steady-state balanced, growth-rate constraint for cellular membrane
    # proteome allocation constraints
    phi_porin * v_P * C_c_mw/porin_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * porin == 0  # proteome allocation constraint for outer membrane porins
    phi_GH16 * v_P * C_c_mw/GH16_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * GH16 == 0     # proteome allocation constraint for extracellular endolaminarase GH16
    phi_TBDT * v_P * C_c_mw/TBDT_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * TBDT == 0     # proteome allocation constraint for TBDT in outer membrane
    phi_GH30 * v_P * C_c_mw/GH30_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * GH30 == 0     # proteome allocation constraint for periplasmic debranching laminarase GH30
    phi_GH17 * v_P * C_c_mw/GH17_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * GH17 == 0     # proteome allocation constraint for periplasmic exolaminarase GH17
    phi_UT * v_P * C_c_mw/UT_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * UT == 0   # proteome allocation constraint for uptake    
    phi_C * v_P * C_c_mw/C_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * CE == 0     # proteome allocation constraint for catabolism
    phi_P * v_P * C_c_mw/P_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * PE == 0     # proteome allocation constraint for protein synthesis
    phi_M * v_P * C_c_mw/M_mw/Qcp * volume_cy*1e-15*Av - exp(logmu) * ME == 0     # proteome allocation constraint for membrane synthesis
    # diffusion constraint
    stoich_matrix[1,1]*v_ENDO*volume_excelrp*1e-15*Av + diff_rate == 0 # number/s, substrate polysaccharide at cell surface
        # upper and lower boundaries of physical parameters
        D_cy <= Dmax                # maximum density constraint; g/L
        D_cy >= Dmin                # minimum density constraint; g/L
        D_pp <= Dmax                # maximum density constraint; g/L
        D_pp >= Dmin                # minimum density constraint; g/L
        r_cy >= 0.05             # minimum cell cytoplasm radius set to 0.05 μm
        r_pp >= 0.01            # minimum cell periplasm radius set to 0.01 μm
        C_to_ATP >= 0.01        # minimum carbon to ATP investment ratio
        C_to_ATP <= 0.99        # maximum carbon to ATP investment ratio
        eP_surf >= 0            # minimum substrate polysaccharide concentration at cell surface
end)
# Solve the optimization problem.
status = JuMP.optimize!(model)