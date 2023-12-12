# PAM optimizer

using Random
using JuMP
using Ipopt
using JLD
using Dates
using CSV
using DataFrames


Random.seed!(41)  # Set a random seed for reproducibility

now = Dates.now()
datetime_str = Dates.format(now, "yyyymmdd_HH")
file_name = "results_$(datetime_str).csv"
df = DataFrame(eP =  Any[], mu_per_hour = Float64[], optmsr_response = Any[], attempts=Int64[],
               density_CYTO= Float64[], density_PERI = Float64[], r_cyt = Float64[], r_pp = Float64[], vol_cy = Float64[], vol_pp = Float64[], vol_exrp=Float64[],
               lm_inner = Float64[], diffusion_rate=Float64[], eP_surf = Float64[], oligoloss =Float64[], porin_O_influx = Float64[], porin_dbO_outflux = Float64[],
               v_ENDO=Float64[], v_IMPRT=Float64[], v_DEBR=Float64[], v_EXO=Float64[],
               v_C=Float64[], v_P=Float64[], v_M=Float64[], v_UT=Float64[],
               eO = Float64[], pO = Float64[], pdbO = Float64[], S_e = Float64[],
               S_c = Float64[], C_c = Float64[], P_c = Float64[], lm = Float64[],
               phi_C = Float64[], phi_P = Float64[], phi_M = Float64[], phi_UT = Float64[], phi_porin = Float64[],
               phi_GH16 = Float64[], phi_TBDT = Float64[], phi_GH30 = Float64[], phi_GH17=Float64[],
               UT = Float64[], C = Float64[], P = Float64[], M = Float64[], porin = Float64[],
               GH16 = Float64[], TBDT = Float64[], GH30 = Float64[], GH17 = Float64[])                  

# environmental concentration of substrate polysaccharide (eP)
stop_log =   log10(0.00000002)
start_log =  log10(0.002)
range_eP_env = 10 .^ range(start_log, stop = stop_log, length = 50) # M, external substrate concentration

for i in 1:length(range_eP_env) 
    attempts = 0  # Initialize the attempts counter
    term_status = MOI.OPTIMIZE_NOT_CALLED  # Initialize the termination status
    optimal_values = zeros(1, 31) # Initialize the optimal values
    while term_status != MOI.LOCALLY_SOLVED && attempts < 8
        mean_for_eP_env = range_eP_env[i] # M, external substrate concentration
        set_value(eP_env, mean_for_eP_env + mean_for_eP_env*0.01*randn()) # M, external substrate concentration
        attempts += 1  # Increment the attempts counter
        # set up optimizer
        JuMP.optimize!(model)

        # Set the callback function
        println("with ext. substrate concentration ", [value(eP_M)], " got growth rate ", [exp.(JuMP.value.(logmu))], 
        " and phis of ", [value(phi_C), value(phi_P), value(phi_M), value(phi_UT)])
        term_status = JuMP.termination_status(model)

        # Save results
        optimal_values_temp =[exp.(JuMP.value.(logmu))*60*60, term_status, attempts,
        value(D_cy), value(D_pp), value(r_cy), value(r_pp), value(volume_cy), value(volume_pp), value(volume_excelrp),
        value(lm_inner), value(diff_rate), value(eP_surf), value(oligoloss_rate), value(porin_O_diff), value(porin_dbO_diff),
        value(v_ENDO), value(v_IMPRT), value(v_DEBR), value(v_EXO),
        value(v_C), value(v_P), value(v_M), value(v_UT),
        value(eO), value(pO), value(pdbO), value(S_e),
        value(S_c), value(C_c), value(P_c), value(lm),
        value(phi_C), value(phi_P), value(phi_M), value(phi_UT), value(phi_porin),
        value(phi_GH16), value(phi_TBDT), value(phi_GH30), value(phi_GH17),
        value(UT), value(CE), value(PE), value(ME), value(porin),
        value(GH16), value(TBDT), value(GH30), value(GH17)]
        if optimal_values_temp[1] >= optimal_values[1]
            optimal_values = optimal_values_temp
        end
    end
    push!(df, vcat(value(eP_env), optimal_values[:]))
end

# Write the results to a CSV file
CSV.write(file_name, df, header=true)
