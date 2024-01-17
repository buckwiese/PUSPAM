# PAM optimizer

using Random
using JuMP
using Ipopt
using JLD
using Dates
using CSV
using DataFrames


Random.seed!(42)  # Set a random seed for reproducibility

now = Dates.now()
datetime_str = Dates.format(now, "yyyymmdd_HH")
file_name = "mixedmodel_results_$(datetime_str).csv"


df = DataFrame(ePlam =  Any[], ePalg = Any[], mu_per_hour = Float64[], optmsr_response = Any[], attempts=Int64[],
            density_CYTO= Float64[], density_PERI = Float64[], r_cyt = Float64[], r_pp = Float64[], vol_cy = Float64[], vol_pp = Float64[], vol_exrp=Float64[],
            lm_inner = Float64[], ePlam_surf = Float64[], ePalg_surf = Float64[], ePlam_surf_rel = Float64[], ePalg_surf_rel = Float64[],
            phi_GH16 = Float64[], phi_TBDTlam = Float64[], phi_PL7 = Float64[], phi_TBDTalg = Float64[], phi_PL5 = Float64[], phi_PL38 = Float64[],
            phi_GH30 = Float64[], phi_GH17 = Float64[], phi_PL6 = Float64[], phi_PL17 = Float64[], 
            phi_UGlc = Float64[], phi_UUro = Float64[], phi_C = Float64[], phi_P = Float64[], phi_M = Float64[], phi_ACB = Float64[], 
            v_UGlc = Float64[], v_UUro = Float64[], v_CGlc = Float64[], v_CPyr = Float64[], v_GGlc = Float64[], v_GPyr = Float64[],
            C_to_ATP = Float64[], oligoloss_rate_lam = Float64[], oligoloss_rate_LOalg = Float64[], oligoloss_rate_SOalg = Float64[]
            )                  

# environmental concentration of laminarin polysaccharide (ePlam)
stop_log =   log10(0.0000001)
start_log =  log10(0.01)
range_ePalg_env = 10 .^ range(start_log, stop = stop_log, length = 30) # M, external substrate concentration
ePlam = 0.0000001 # M, external alginate concentration

for i in 1:length(range_ePalg_env) 
    attempts = 0  # Initialize the attempts counter
    term_status = MOI.OPTIMIZE_NOT_CALLED  # Initialize the termination status
    optimal_values = zeros(1, 31) # Initialize the optimal values
    while term_status != MOI.LOCALLY_SOLVED && attempts < 5
        mean_for_ePalg_env = range_ePalg_env[i] # M, external substrate concentration
        set_value(ePalg_env, mean_for_ePalg_env + mean_for_ePalg_env*0.01*randn()) # M, external substrate concentration
        set_value(ePlam_env, ePlam) # M, external alginate concentration
        attempts += 1  # Increment the attempts counter
        # set up optimizer
        JuMP.optimize!(model)

        # Set the callback function
        println("with ext. substrate concentration ", [value(mean_for_ePalg_env)], " got growth rate ", [exp.(JuMP.value.(logmu))], 
        " and phis of ", [value(phi_C), value(phi_P), value(phi_M)])
        term_status = JuMP.termination_status(model)

        # Save results
        optimal_values_temp =[exp.(JuMP.value.(logmu))*60*60, term_status, attempts,
        value(D_cy), value(D_pp), value(r_cy), value(r_pp), value(volume_cy), value(volume_pp), value(volume_excelrp),
        value(lm_inner), value(ePlam_surf), value(ePalg_surf), value(ePlam_surf_rel), value(ePalg_surf_rel),
        value(phi_GH16), value(phi_TBDTlam), value(phi_PL7), value(phi_TBDTalg),
        value(phi_PL5), value(phi_PL38), value(phi_GH30), value(phi_GH17), value(phi_PL6), value(phi_PL17),
        value(phi_UGlc), value(phi_UUro), value(phi_C), value(phi_P), value(phi_M), value(phi_ACB),
        value(v_UGlc), value(v_UUro), value(v_CGlc), value(v_CPyr), value(v_GGlc), value(v_GPyr),
        value(C_to_ATP), value(oligoloss_rate_lam), value(oligoloss_rate_LM)+value(oligoloss_rate_LG), value(oligoloss_rate_SM)+value(oligoloss_rate_SG)
        ]
        if optimal_values_temp[1] >= optimal_values[1]
            optimal_values = optimal_values_temp
        end
    end
    push!(df, vcat(value(ePlam_env), value(ePalg_env), optimal_values[:]))
end

# Write the results to a CSV file
CSV.write(file_name, df, header=true)
    

