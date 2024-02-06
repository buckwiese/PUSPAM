# PAM optimizer

using Random
using JuMP
using Ipopt
using JLD
using Dates
using CSV
using DataFrames


#Random.seed!(20)  # Set a random seed for reproducibility

now = Dates.now()
datetime_str = Dates.format(now, "yyyymmdd_HH")
file_name = "mixedmodel_results_$(datetime_str).csv"


df = DataFrame(ePlam =  Any[], ePalg = Any[], distance = Float64[], mu_per_hour = Float64[], optmsr_response = Any[], attempts=Int64[],
            density_CYTO= Float64[], density_PERI = Float64[], r_cyt = Float64[], r_pp = Float64[], vol_cy = Float64[], vol_pp = Float64[], vol_exrp=Float64[],
            lm_inner = Float64[], ePlam_surf = Float64[], ePalg_surf = Float64[], ePlam_surf_rel = Float64[], ePalg_surf_rel = Float64[],
            phi_GH16A = Float64[], phi_GH16B = Float64[], phi_GH16C = Float64[], phi_TBDTlam = Float64[],
            phi_PL7A = Float64[], phi_PL7B = Float64[], phi_PL7C = Float64[],  
            phi_PL5A = Float64[], phi_PL5B = Float64[], phi_PL5C = Float64[], phi_PL38A = Float64[], phi_PL38B = Float64[], phi_PL38C = Float64[], phi_TBDTalg = Float64[],
            phi_GH30A = Float64[], phi_GH30B = Float64[], phi_GH17A = Float64[], phi_GH17B = Float64[], phi_PL6A = Float64[], phi_PL6B = Float64[], phi_PL17A = Float64[], phi_PL17B = Float64[],
            phi_UGlc = Float64[], phi_UUro = Float64[], phi_C = Float64[], phi_P = Float64[], phi_M = Float64[], phi_ACB = Float64[], 
            v_UGlc = Float64[], v_UUro = Float64[], v_CGlc = Float64[], v_CPyr = Float64[], v_GGlc = Float64[], v_GPyr = Float64[],
            C_to_ATP = Float64[], oligoloss_rate_lam = Float64[], oligoloss_rate_LOalg = Float64[], oligoloss_rate_SOalg = Float64[],
            eOlam_rel = Float64[], pOlam_rel = Float64[], pdbOlam_rel = Float64[], pGlc_rel = Float64[], 
            eLMO_rel = Float64[], eLGO_rel = Float64[], eSMO_rel = Float64[], eSGO_rel = Float64[], 
            pSMO_rel = Float64[], pSGO_rel = Float64[], pU_rel = Float64[], 
            cGlc_rel = Float64[], cU_rel = Float64[], cPyr_rel = Float64[], cC_rel = Float64[], cP_rel = Float64[], lm_rel = Float64[],
            kcat_GH16A = Float64[], km_GH16A = Float64[], kcat_GH16B = Float64[], km_GH16B = Float64[], kcat_GH16C = Float64[], km_GH16C = Float64[],
            kcat_PL7A = Float64[],  km_PL7A = Float64[],  kcat_PL7B = Float64[],  km_PL7B = Float64[],  kcat_PL7C = Float64[],  km_PL7C = Float64[],
            kcat_PL5A = Float64[],  km_PL5A = Float64[],  kcat_PL5B = Float64[],  km_PL5B = Float64[],  kcat_PL5C = Float64[],  km_PL5C = Float64[],
            kcat_PL38A = Float64[], km_PL38A = Float64[], kcat_PL38B = Float64[], km_PL38B = Float64[], kcat_PL38C = Float64[], km_PL38C = Float64[],
            kcat_GH30A = Float64[], km_GH30A = Float64[], kcat_GH30B = Float64[], km_GH30B = Float64[],
            kcat_GH17A = Float64[], km_GH17A = Float64[], kcat_GH17B = Float64[], km_GH17B = Float64[],
            kcat_PL6A = Float64[],  km_PL6A = Float64[],  kcat_PL6B = Float64[],  km_PL6B = Float64[],
            kcat_PL17A = Float64[], km_PL17A = Float64[], kcat_PL17B = Float64[], km_PL17B = Float64[]
            )                  

# environmental distance between bacterial cells (micrometers)
start_log =   log10(0.5)
stop_log =  log10(5000)
range_cell_dist = 10 .^ range(start_log, stop = stop_log, length = 15) # M, external substrate concentration
ePlam = 0.000005 # M, external laminarin concentration
ePalg = 0.0 # M, external alginate concentration

for i in 1:length(range_cell_dist) 
    attempts = 0  # Initialize the attempts counter
    term_status = MOI.OPTIMIZE_NOT_CALLED  # Initialize the termination status
    optimal_values = zeros(1, 31) # Initialize the optimal values
    while term_status != MOI.LOCALLY_SOLVED && attempts < 2
        set_value(ePlam_env, ePlam) # M, external substrate concentration
        set_value(ePalg_env, ePalg) # M, external alginate concentration
        set_value(cell_dist, range_cell_dist[i]) # micrometers, distance between bacterial cells
        attempts += 1  # Increment the attempts counter
        # set up optimizer
        JuMP.optimize!(model)

        # Set the callback function
        println("with cell distance of ", [value(range_cell_dist[i])], " got growth rate ", [exp.(JuMP.value.(logmu))], 
        " and phis of ", [value(phi_C), value(phi_P), value(phi_M)])
        term_status = JuMP.termination_status(model)

        # Save results
        optimal_values_temp =[exp.(JuMP.value.(logmu))*60*60, term_status, attempts,
        value(D_cy), value(D_pp), value(r_cy), value(r_pp), value(volume_cy), value(volume_pp), value(volume_excelrp),
        value(lm_inner), value(ePlam_surf), value(ePalg_surf), value(ePlam_surf_rel), value(ePalg_surf_rel),
        value(phi_GH16A), value(phi_GH16B), value(phi_GH16C), value(phi_TBDTlam), 
        value(phi_PL7A), value(phi_PL7B), value(phi_PL7C),
        value(phi_PL5A), value(phi_PL5B), value(phi_PL5C), value(phi_PL38A), value(phi_PL38B), value(phi_PL38C), value(phi_TBDTalg),
        value(phi_GH30A), value(phi_GH30B), value(phi_GH17A), value(phi_GH17B), value(phi_PL6A), value(phi_PL6B), value(phi_PL17A), value(phi_PL17B),
        value(phi_UGlc), value(phi_UUro), value(phi_C), value(phi_P), value(phi_M), value(phi_ACB),
        value(v_UGlc), value(v_UUro), value(v_CGlc), value(v_CPyr), value(v_GGlc), value(v_GPyr),
        value(C_to_ATP), value(oligoloss_rate_lam), value(oligoloss_rate_LM)+value(oligoloss_rate_LG), value(oligoloss_rate_SM)+value(oligoloss_rate_SG),
        value(eOlam_rel), value(pOlam_rel), value(pdbOlam_rel), value(pGlc_rel),
        value(eLMO_rel), value(eLGO_rel), value(eSMO_rel), value(eSGO_rel),
        value(pSMO_rel), value(pSGO_rel), value(pU_rel),
        value(cGlc_rel), value(cU_rel), value(cPyr_rel), value(cC_rel), value(cP_rel), value(lm_rel),
        value(GH16A_kcat), value(GH16A_km), value(GH16B_kcat), value(GH16B_km), value(GH16C_kcat), value(GH16C_km),
        value(PL7A_kcat), value(PL7A_km), value(PL7B_kcat), value(PL7B_km), value(PL7C_kcat), value(PL7C_km),
        value(PL5A_kcat), value(PL5A_km), value(PL5B_kcat), value(PL5B_km), value(PL5C_kcat), value(PL5C_km),
        value(PL38A_kcat), value(PL38A_km), value(PL38B_kcat), value(PL38B_km), value(PL38C_kcat), value(PL38C_km),
        value(GH30A_kcat), value(GH30A_km), value(GH30B_kcat), value(GH30B_km),
        value(GH17A_kcat), value(GH17A_km), value(GH17B_kcat), value(GH17B_km),
        value(PL6A_kcat), value(PL6A_km), value(PL6B_kcat), value(PL6B_km),
        value(PL17A_kcat), value(PL17A_km), value(PL17B_kcat), value(PL17B_km)
        ]
        if optimal_values_temp[1] >= optimal_values[1]
            optimal_values = optimal_values_temp
        end
    end
    push!(df, vcat(value(ePlam_env), value(ePalg_env), value(range_cell_dist[i]), optimal_values[:]))
end

# Write the results to a CSV file
CSV.write(file_name, df, header=true)
    

