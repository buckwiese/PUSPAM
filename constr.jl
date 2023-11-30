using JuMP, Ipopt

#eq for PAM optimizer
# set up model
model = Model(Ipopt.Optimizer)
#set_optimizer_attribute(model, "max_iter", 100)

# Set the environment:
#@NLparameter(model, T   == 20.0 + 273.0)   # K,  temperature
@NLparameter(model, S_e == 1) # M,  external substrate concentration
println("S_e = ", value(S_e))

#variables to optimize
#concentrations in M
@variables(model, begin
    logmu, (start = 0.00001) # 1/s, growth rate
    S_c >= 1, (start = 1e2)   # number/μm3, intracellular substrate concentration
    C_c >= 1, (start = 1e2)   # number/μm3, cellular carbon concentration
    P_c >= 1, (start = 1e5)   # number/μm3, cellular protein concentration
    lm >= 1, (start = 1e6)    # number/μm3, cellular membrane concentration
    UT >= 1, (start = 10)   # number/μm3, uptake transporters
    CE >= 1, (start = 10)   # number/μm3, catabolism enzymes
    PE >= 1, (start = 10)   # number/μm3, protein synthesis enzymes
    ME >= 1, (start = 10)   # number/μm3, membrane synthesis enzymes
    phi_C >= 0.05, (start = 0.25) # relative proteome allocation to catabolism
    phi_P >= 0.05, (start = 0.25) # relative proteome allocation to protein synthesis
    phi_M >= 0.05, (start = 0.25) # relative proteome allocation to membrane synthesis
    phi_UT >= 0.05, (start = 0.25) # relative proteome allocation to uptake
    beta >= 0, (start = 0.1) # μm, cell volume to surface area ratio
end)
    
# Set the objective: maximize the log growth rate.
@objective(model, Max, logmu)

# Define intermediate calculations.
@NLexpressions(model, begin
    v_UT, UT_kcat *UT* S_e / (UT_km + S_e) # M/s, uptake
    v_C, C_kcat * CE * S_c / (C_km + S_c) # M/s, catabolism
    v_P, P_kcat*PE * C_c / (P_km + C_c) # M/s, protein synthesis
    v_M, M_kcat * ME * C_c / (M_km + C_c) # M/s, membrane synthesis
    D, P_c*C_c_mw + lm*lm_mw + S_c*S_c_mw + C_c*C_c_mw # fg/μm3, cell density
    area, 4*pi*(beta*3)^2 # μm2, cell surface area
    volume, 4/3*pi*(beta*3)^3 # μm3, cell volume
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
#    area-lm*sa_lm-UT*sa_UT == 0.0                          # cell volume to surface area ratio constraint
    matrix[1, 1] * v_UT + matrix[1, 2] * v_C + matrix[1, 3] * v_P + matrix[1, 4] * v_M - exp(logmu) * S_c == 0 # steady-state balanced, growth-rate constraint for intracellular substrate
    matrix[2, 1] * v_UT + matrix[2, 2] * v_C + matrix[2, 3] * v_P + matrix[2, 4] * v_M - exp(logmu) * C_c == 0 # steady-state balanced, growth-rate constraint for cellular carbon
    matrix[3, 1] * v_UT + matrix[3, 2] * v_C + matrix[3, 3] * v_P + matrix[3, 4] * v_M - exp(logmu) * P_c == 0 # steady-state balanced, growth-rate constraint for cellular proteins
    matrix[4, 1] * v_UT + matrix[4, 2] * v_C + matrix[4, 3] * v_P + matrix[4, 4] * v_M - exp(logmu) * lm == 0 # steady-state balanced, growth-rate constraint for membrane units
    phi_C * v_P *C_c_mw/C_mw - exp(logmu) * CE == 0 # proteome allocation constraint for catabolism
    phi_P * v_P *C_c_mw/P_mw - exp(logmu) * PE == 0 # proteome allocation constraint for protein synthesis
    phi_M * v_P *C_c_mw/M_mw - exp(logmu) * ME == 0 # proteome allocation constraint for membrane synthesis
    phi_UT * v_P *C_c_mw/UT_mw - exp(logmu) * UT == 0 # proteome allocation constraint for uptake
    D ≤ Dmax                                                  # density constraint; Da/μm3
    D ≥ Dmin                                                  # density constraint; Da/μm3
    beta >= 0.001                                             # cell volume to surface area ratio constraint
    beta <= 10.0                                              # cell volume to surface area ratio constraint
    sqrt((lm*sa_lm+sa_UT * UT)/(4*pi)) >= 0.02                 # minimum cell radius set to 0.01 μm
end)

# Solve the optimization problem.
status = JuMP.optimize!(model)