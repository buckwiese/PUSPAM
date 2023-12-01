using JuMP, Ipopt

#eq for PAM optimizer
# set up model
model = Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "max_iter", 1000)

# Set the environment:
#@NLparameter(model, T   == 20.0 + 273.0)   # K,  temperature
@NLparameter(model, S_e_M == 1) # M,  external substrate concentration
println("S_e = ", value(S_e_M))

#variables to optimize
#concentrations in M
@variables(model, begin
    logmu, (start = 0.00001) # 1/s, growth rate
    S_c >= 1, (start = 1.5e5)   # number/μm3, intracellular substrate concentration
    C_c >= 1, (start = 1.5e5)   # number/μm3, cellular carbon concentration
    P_c >= 1, (start = 2e5)   # number/μm3, cellular protein concentration
    lm >= 1, (start = 1.5e6)    # number/μm3, cellular membrane concentration
#    UT >= 1, (start = 7000)   # number/μm3, uptake transporters
#    CE >= 1, (start = 400)   # number/μm3, catabolism enzymes
#    PE >= 1, (start = 1400)   # number/μm3, protein synthesis enzymes
#    ME >= 1, (start = 100)   # number/μm3, membrane synthesis enzymes
    phi_C >= 0.05, (start = 0.25) # relative proteome allocation to catabolism
    phi_P >= 0.05, (start = 0.25) # relative proteome allocation to protein synthesis
    phi_M >= 0.05, (start = 0.1) # relative proteome allocation to membrane synthesis
    phi_UT >= 0.05, (start = 0.4) # relative proteome allocation to uptake
    beta >= 0, (start = 0.04) # μm, cell volume to surface area ratio
end)
    
# Set the objective: maximize the log growth rate.
@objective(model, Max, logmu)

# Define intermediate calculations.
@NLexpressions(model, begin
    area, 4*pi*(beta*3)^2       # μm2, cell surface area
    volume, 4/3*pi*(beta*3)^3   # μm3, cell volume
    r, sqrt(area/(4*pi))        # μm, cell radius
    S_c_M, S_c*1e15/Av     # M, intracellular substrate concentration
    C_c_M, C_c*1e15/Av     # M, cellular carbon concentration
    P_c_M, P_c*1e15/Av     # M, cellular protein concentration
    UT, phi_UT*P_c*C_c_mw/UT_mw/Qcp   # number/μm3, uptake transporter concentration
    CE, phi_C*P_c*C_c_mw/C_mw/Qcp     # number/μm3, catabolism enzyme concentration
    PE, phi_P*P_c*C_c_mw/P_mw/Qcp     # number/μm3, protein synthesis enzyme concentration
    ME, phi_M*P_c*C_c_mw/M_mw/Qcp     # number/μm3, membrane synthesis enzyme concentration
    lm_M, lm*1e15/Av       # M, cellular membrane concentration
    UT_M, UT*1e15/Av       # M, uptake transporter concentration
    CE_M, CE*1e15/Av       # M, catabolism enzyme concentration
    PE_M, PE*1e15/Av       # M, protein synthesis enzyme concentration
    ME_M, ME*1e15/Av       # M, membrane synthesis enzyme concentration    
    v_UT, UT_kcat * UT_M * S_e_M / (UT_km + S_e_M)    # M/s, uptake
    v_C, C_kcat * CE_M * S_c_M / (C_km + S_c_M)       # M/s, catabolism
    v_P, P_kcat * PE_M * C_c_M / (P_km + C_c_M)       # M/s, protein synthesis
    v_M, M_kcat * ME_M * C_c_M / (M_km + C_c_M)       # M/s, membrane synthesis
    D, (P_c_M * C_c_mw + lm_M * lm_mw + S_c_M * S_c_mw + C_c_M * C_c_mw)/0.3  # g/L, cell density, assuming c makes up 30% of cell mass
end)
    
# Define the constraints.
# Linear constraints
@constraints(model, begin
    phi_C + phi_P + phi_M + phi_UT == 1 # proteome allocation constraint
end)

# Nonlinear constraints
@NLconstraints(model, begin
    UT / lm >= Mmin                                        # constraint on the transporter to lipid ratio in the membrane
    UT / lm <= Mmax                                        # constraint on the transporter to lipid ratio in the membrane
    beta * (sa_UT * UT + sa_lm * lm ) == 1.0               # cell volume to surface area ratio constraint
    area-sa_lm*lm*volume-sa_UT*UT*volume == 0.0            # constraint on the cell surface area made up of lipids and transporters
    matrix[1, 1] * v_UT + matrix[1, 2] * v_C + matrix[1, 3] * v_P + matrix[1, 4] * v_M - exp(logmu) * S_c_M == 0  # steady-state balanced, growth-rate constraint for intracellular substrate
    matrix[2, 1] * v_UT + matrix[2, 2] * v_C + matrix[2, 3] * v_P + matrix[2, 4] * v_M - exp(logmu) * C_c_M == 0  # steady-state balanced, growth-rate constraint for cellular carbon
    matrix[3, 1] * v_UT + matrix[3, 2] * v_C + matrix[3, 3] * v_P + matrix[3, 4] * v_M - exp(logmu) * P_c_M == 0  # steady-state balanced, growth-rate constraint for cellular proteins
    matrix[4, 1] * v_UT + matrix[4, 2] * v_C + matrix[4, 3] * v_P + matrix[4, 4] * v_M - exp(logmu) * lm_M == 0   # steady-state balanced, growth-rate constraint for membrane units
#    phi_C * v_P *C_c_mw/C_mw/Qcp - exp(logmu) * CE_M == 0     # proteome allocation constraint for catabolism
#    phi_P * v_P *C_c_mw/P_mw/Qcp - exp(logmu) * PE_M == 0     # proteome allocation constraint for protein synthesis
#    phi_M * v_P *C_c_mw/M_mw/Qcp - exp(logmu) * ME_M == 0     # proteome allocation constraint for membrane synthesis
#    phi_UT * v_P *C_c_mw/UT_mw/Qcp - exp(logmu) * UT_M == 0   # proteome allocation constraint for uptake
    P_c_M * C_c_mw - CE_M * C_mw * Qcp - PE_M * P_mw * Qcp - ME_M * M_mw * Qcp - UT_M * UT_mw * Qcp <= 0 # proteome mass balance constraint    
    D ≤ Dmax                # maximum density constraint; g/L
    D ≥ Dmin                # density constraint; g/L
    beta >= 0.01            # cell volume to surface area ratio constraint
    beta <= 0.5             # cell volume to surface area ratio constraint
    r >= 0.1                # minimum cell radius set to 0.1 μm
end)

# Solve the optimization problem.
status = JuMP.optimize!(model)