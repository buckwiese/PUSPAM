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
    logmu, (start = -3.2) # 1/s, growth rate
    eO >= 0, (start = 800)   # number, extracellular oligosaccharide 
    pO >= 0, (start = 800000)   # number, periplasmic oligosaccharide 
    pdbO >= 0, (start = 500000) # number, periplasmic debranched oligosaccharide
    S_e >= 0, (start = 12000)  # number, periplasmic sugar
    S_c >= 1, (start = 2e7) # number, intracellular substrate
    C_c >= 1, (start = 8e7) # number, cellular carbon 
    P_c >= 1, (start = 1e8)   # number, cellular protein 
    lm >= 1, (start = 8e8)  # number, cellular membrane
    porin >= 0, (start = 1) # number, outer membrane porin channels
    GH16 >= 0, (start = 3e7)   # number, extracellular endolaminarase GH16
    TBDT >= 0, (start = 400000)   # number, TBDT in outer membrane
    GH30 >= 0, (start = 600000)   # number, periplasmic debranching laminarase GH30
    GH17 >= 0, (start = 180000)   # number, periplasmic exolaminarase GH17    
    UT >= 0, (start = 1e6)     # number, uptake transporters
    CE >= 0, (start = 150000)     # number, catabolism enzymes
    PE >= 0, (start = 190000)    # number, protein synthesis enzymes
    ME >= 0, (start = 60000)     # number, membrane synthesis enzymes
    phi_porin >= 0.0, (start = 0.000001) # relative proteome allocation to outer membrane porins
    phi_GH16 >=0.0, (start = 0.2) # relative proteome allocation to extracellular endolaminarase GH16
    phi_TBDT >=0.0, (start = 0.015) # relative proteome allocation to TBDT in outer membrane
    phi_GH30 >=0.0, (start = 0.006) # relative proteome allocation to periplasmic debranching laminarase GH30
    phi_GH17 >=0.0, (start = 0.003) # relative proteome allocation to periplasmic exolaminarase GH17
    phi_C >= 0.0, (start = 0.13) # relative proteome allocation to catabolism
    phi_P >= 0.0, (start = 0.08) # relative proteome allocation to protein synthesis
    phi_M >= 0.0, (start = 0.027) # relative proteome allocation to membrane synthesis
    phi_UT >= 0.0, (start = 0.04) # relative proteome allocation to uptake
    beta >= 0, (start = 0.5) # μm, cell volume to surface area ratio
    f_pp >= 0, (start = 0.2) # fraction of cell radius contributed by periplasm
    lm_inner >= 0.1, (start = 0.5) # fraction of lipids in inner membrane
    eP_surf >= 0, (start = 1e-7) # M, substrate polysaccharide at cell surface
end)
    
# Set the objective: maximize the log growth rate.
@objective(model, Max, logmu)

# Define intermediate calculations.
@NLexpressions(model, begin
    # cell parameters
    area_cy, 4*pi*(beta*3)^2       # μm2, cytoplasm surface area
    volume_cy, 4/3*pi*(beta*3)^3   # μm3, cytoplasm volume
    r_cy, beta*3                   # μm, cell radius
    r_pp, r_cy*(f_pp/(1-f_pp))     # μm, periplasm thickness
    area_pp, 4*pi*(r_pp+r_cy)^2    # μm2, periplasm surface area
    volume_pp, 4/3*pi*(r_pp+r_cy)^3 - volume_cy     # μm3, periplasm volume
    volume_excelrp, 4/3*pi*(r_pp+r_cy+0.05)^3 - volume_cy - volume_pp     # μm3, extracellular reaction space volume
    # metabolite concentrations
    eO_M, eO/volume_excelrp*1e15/Av     # M, extracellular oligosaccharide concentration
    pO_M, pO/volume_pp*1e15/Av     # M, periplasmic oligosaccharide concentration
    pdbO_M, pdbO/volume_pp*1e15/Av   # M, periplasmic debranched oligosaccharide concentration
    S_e_M, S_e/volume_pp*1e15/Av     # M, periplasmic sugar concentration
    S_c_M, S_c/volume_cy*1e15/Av     # M, intracellular substrate concentration
    C_c_M, C_c/volume_cy*1e15/Av     # M, cellular carbon concentration
    P_c_M, P_c/volume_cy*1e15/Av     # M, cellular protein concentration
    lm_M, lm/volume_cy*1e15/Av       # M, cellular membrane concentration
    # CAZyme numbers and distributions
    GH16_M, GH16/volume_excelrp*1e15/Av       # M, endolaminarase GH16 concentration in extracellular reaction space
    TBDT_M, TBDT/volume_excelrp*1e15/Av       # M, TBDT concentration in outer membrane
    GH30_M, GH30/volume_pp*1e15/Av       # M, debranching laminarase GH30 concentration in periplasm
    GH17_M, GH17/volume_pp*1e15/Av       # M, exolaminarase GH17 concentration in periplasm
    UT_M, UT/volume_pp*1e15/Av       # M, uptake transporter concentration in periplasm
    CE_M, CE/volume_cy*1e15/Av       # M, catabolism enzyme concentration in cytoplasm
    PE_M, PE/volume_cy*1e15/Av       # M, protein synthesis enzyme concentration in cytoplasm
    ME_M, ME/volume_cy*1e15/Av       # M, membrane synthesis enzyme concentration in cytoplasm
    # reaction rates
    v_ENDO, GH16_kcat * GH16_M * eP_surf / (GH16_km + eP_surf)    # M/s, extracellular hydrolysis
    v_IMPRT, TBDT_kcat * TBDT_M * eO_M / (TBDT_km + eO_M)    # M/s, uptake
    v_DEBR, GH30_kcat * GH30_M * pO_M / (GH30_km + pO_M)    # M/s, periplasmic debranching
    v_EXO, GH17_kcat * GH17_M * pdbO_M / (GH17_km + pdbO_M)    # M/s, periplasmic hydrolysis
    v_UT, UT_kcat * UT_M * S_e_M / (UT_km + S_e_M)    # M/s, uptake
    v_C, C_kcat * CE_M * S_c_M / (C_km + S_c_M)       # M/s, catabolism
    v_P, P_kcat * PE_M * C_c_M / (P_km + C_c_M)       # M/s, protein synthesis
    v_M, M_kcat * ME_M * C_c_M / (M_km + C_c_M)       # M/s, membrane synthesis
    D_cy, (UT*UT_mw + CE*C_mw + PE*P_mw + ME*M_mw + lm*lm_mw*lm_inner + C_c*C_c_mw*3 + S_c*S_c_mw + pO*pO_mw + pdbO*pdbO_mw)/(Av*volume_cy)*1e15  # g/L, cell density, assuming c makes up 30% of cell mass
    D_pp, (porin*porin_mw + GH16*GH16_mw + TBDT*TBDT_mw + GH30*GH30_mw + GH17*GH17_mw + (1-lm_inner)*lm*lm_mw + pO*pO_mw + pdbO*pdbO_mw)/(Av*volume_pp)*1e15  # g/L, cell density, assuming c makes up 30% of cell mass
    # diffusion processes
    diff_rate, DC_eP* (eP_env-eP_surf)*1e-15/50*(4*pi*(r_cy+r_pp+0.05)^2) # mol/s, diffusion rate of substrate polysaccharide (eP) to cell surface based on Karp-Boss et al. 1996 
    oligoloss_rate, DC_O*4*pi*(r_cy+r_pp+0.05)^2*(eO_M)*1e-15/50 # mol/s, diffusion rate of oligosaccharide (eO) from cell surface based on Karp-Boss et al. 1996
    porin_O_diff, DC_O*porin_pore*porin*(eO_M-pO_M)*1e-15/((r_pp+0.05)/2) # mol/s, diffusion rate of extracellular oligosaccharides (eO) through porins into the periplasm
    porin_dbO_diff, DC_O*sa_porin*porin*(pdbO_M)*1e-15/((r_pp+0.05)/2) # mol/s, diffusion rate of extracellular oligosaccharides (eO) through porins into the periplasm
end)
    
# Define the constraints.
# Linear constraints
@constraints(model, begin
    phi_C + phi_P + phi_M + phi_UT + phi_porin + phi_GH16 + phi_TBDT + phi_GH30 + phi_GH17 == 0.5 # proteome allocation constraint
end)

# Nonlinear constraints
@NLconstraints(model, begin
    UT / (lm*lm_inner) >= Mmin                                        # constraint on the transporter to lipid ratio in the membrane
    UT / (lm*lm_inner) <= Mmax                                        # constraint on the transporter to lipid ratio in the membrane
    TBDT/ (lm*(1-lm_inner)) >= Mmin                          # constraint on the transporter to lipid ratio in the membrane
    (TBDT+porin)/ (lm*(1-lm_inner)) <= Mmax                          # constraint on the transporter to lipid ratio in the membrane
    lm_inner + (1-lm_inner) == 1.0                              # constraint on the transporter to lipid ratio in the membrane
    beta * (sa_UT * UT/volume_cy + sa_lm * lm_inner*lm/volume_cy) == 1                  # cell volume to surface area ratio constraint
    area_cy - (sa_lm*lm*lm_inner + sa_UT*UT ) == 0.0    # constraint on the cell surface area made up of lipids and transporters
    area_pp - (sa_lm*lm*(1-lm_inner) + sa_TBDT*TBDT + sa_porin*porin) == 0.0    # constraint on the cell surface area made up of lipids and transporters
    sugar_matrix[2, 1]*v_ENDO*Av*volume_excelrp/1e15 + sugar_matrix[2, 2]*v_IMPRT*Av*volume_excelrp/1e15 + sugar_matrix[2, 3]*v_DEBR*Av*volume_excelrp/1e15 + sugar_matrix[2, 4]*v_EXO*Av*volume_excelrp/1e15 - exp(logmu)*eO - oligoloss_rate*Av - porin_O_diff*Av == 0  # steady-state balanced, growth-rate constraint for extracellular oligosaccharide
    sugar_matrix[3, 1]*v_ENDO*Av*volume_pp/1e15 + sugar_matrix[3, 2]*v_IMPRT*Av*volume_pp/1e15 + sugar_matrix[3, 3]*v_DEBR*Av*volume_pp/1e15 + sugar_matrix[3, 4]*v_EXO*Av*volume_pp/1e15 - exp(logmu) * pO + porin_O_diff*Av == 0  # steady-state balanced, growth-rate constraint for periplasmic oligosaccharide
    sugar_matrix[4, 1]*v_ENDO*Av*volume_pp/1e15 + sugar_matrix[4, 2]*v_IMPRT*Av*volume_pp/1e15 + sugar_matrix[4, 3]*v_DEBR*Av*volume_pp/1e15 + sugar_matrix[4, 4]*v_EXO*Av*volume_pp/1e15 - exp(logmu) * pdbO - porin_dbO_diff*Av == 0  # steady-state balanced, growth-rate constraint for periplasmic debranched oligosaccharide
    sugar_matrix[5, 1]*v_ENDO*Av*volume_pp/1e15 + sugar_matrix[5, 2]*v_IMPRT*Av*volume_pp/1e15 + sugar_matrix[5, 3]*v_DEBR*Av*volume_pp/1e15 + sugar_matrix[5, 4]*v_EXO*Av*volume_pp/1e15 + matrix[5,1]*v_UT*Av*volume_pp/1e15 - exp(logmu) * S_e == 0  # steady-state balanced, growth-rate constraint for periplasmic sugars
    matrix[1, 1]*v_UT*Av*volume_cy/1e15 + matrix[1, 2]*v_C*Av*volume_cy/1e15 + matrix[1, 3]*v_P*Av*volume_cy/1e15 + matrix[1, 4]*v_M*Av*volume_cy/1e15 - exp(logmu) * S_c == 0  # steady-state balanced, growth-rate constraint for intracellular substrate
    matrix[2, 1]*v_UT*Av*volume_cy/1e15 + matrix[2, 2]*v_C*Av*volume_cy/1e15 + matrix[2, 3]*v_P*Av*volume_cy/1e15 + matrix[2, 4]*v_M*Av*volume_cy/1e15 - exp(logmu) * C_c == 0  # steady-state balanced, growth-rate constraint for cellular carbon
    matrix[3, 1]*v_UT*Av*volume_cy/1e15 + matrix[3, 2]*v_C*Av*volume_cy/1e15 + matrix[3, 3]*v_P*Av*volume_cy/1e15 + matrix[3, 4]*v_M*Av*volume_cy/1e15 - exp(logmu) * P_c == 0  # steady-state balanced, growth-rate constraint for cellular protein
    matrix[4, 1]*v_UT*Av*volume_cy/1e15 + matrix[4, 2]*v_C*Av*volume_cy/1e15 + matrix[4, 3]*v_P*Av*volume_cy/1e15 + matrix[4, 4]*v_M*Av*volume_cy/1e15 - exp(logmu) * lm == 0  # steady-state balanced, growth-rate constraint for cellular membrane
    phi_porin * P_c *C_c_mw/porin_mw/Qcp - exp(logmu) * porin == 0
    phi_GH16 * P_c *C_c_mw/GH16_mw/Qcp - exp(logmu) * GH16 == 0     # proteome allocation constraint for extracellular endolaminarase GH16
    phi_TBDT * P_c *C_c_mw/TBDT_mw/Qcp - exp(logmu) * TBDT == 0     # proteome allocation constraint for TBDT in outer membrane
    phi_GH30 * P_c *C_c_mw/GH30_mw/Qcp - exp(logmu) * GH30 == 0     # proteome allocation constraint for periplasmic debranching laminarase GH30
    phi_GH17 * P_c *C_c_mw/GH17_mw/Qcp - exp(logmu) * GH17 == 0     # proteome allocation constraint for periplasmic exolaminarase GH17
    phi_UT * P_c *C_c_mw/UT_mw/Qcp - exp(logmu) * UT == 0   # proteome allocation constraint for uptake    
    phi_C * P_c *C_c_mw/C_mw/Qcp - exp(logmu) * CE == 0     # proteome allocation constraint for catabolism
    phi_P * P_c *C_c_mw/P_mw/Qcp - exp(logmu) * PE == 0     # proteome allocation constraint for protein synthesis
    phi_M * P_c *C_c_mw/M_mw/Qcp - exp(logmu) * ME == 0     # proteome allocation constraint for membrane synthesis
    diff_rate/(volume_excelrp*1e-15) - v_ENDO == 0 # steady-state balanced, diffusion-compensated hydrolysis of substrate polysaccharide (eP) at cell surface
    D_cy <= Dmax/2                # maximum density constraint; g/L
    D_pp <= Dmax/2                # maximum density constraint; g/L
    D_cy >= Dmin/2                # density constraint; g/L
    D_pp >= Dmin/2                # density constraint; g/L
    beta >= 0.01            # cell volume to surface area ratio constraint
    beta <= 0.5             # cell volume to surface area ratio constraint
    r_cy >= 0.05            # minimum cell cytoplasm radius set to 0.05 μm
    r_pp >= 0.01            # minimum cell periplasm radius set to 0.01 μm
end)

# Solve the optimization problem.
status = JuMP.optimize!(model)