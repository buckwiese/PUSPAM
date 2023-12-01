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
file_name = "results_$(datetime_str).csv"
df = DataFrame(S_e_M =  Any[], mu_per_hour = Float64[], S_c_M = Float64[], C_c_M = Float64[], P_c_M = Float64[], lm_M = Float64[], 
                        S_c_per_cell =Float64[], C_c_per_cell =Float64[], P_c_per_cell =Float64[], lm_per_cell =Float64[],  
                        phi_C = Float64[], phi_P = Float64[], phi_M = Float64[], phi_UT = Float64[], beta = Float64[], 
#                        UT_per_cell = Float64[], CE_per_cell = Float64[], PE_per_cell = Float64[], ME_per_cell = Float64[], 
                        v_UT = Float64[], v_C = Float64[], v_P = Float64[], v_M = Float64[], 
                        D = Float64[], area = Float64[], volume = Float64[], r=Float64[])

# Load values specifying proteome allocation problem
#PAM_parameters=load(PAM_parameter_dict) # this can go, Julia remembers!

start_log = log10(0.0000001)
stop_log = log10(0.1)
S_e_change = 10 .^ range(start_log, stop = stop_log, length = 20) # M, external substrate concentration       

for i in 1:length(S_e_change)
    set_value(S_e_M, S_e_change[i])
    optimal_values = zeros(1, 26)
    for j in 1:10
        # set up optimizer
        JuMP.optimize!(model)

        # Set the callback function
        println("with ext. substrate concentration ", [value(S_e_M)], " got growth rate ", [exp.(JuMP.value.(logmu))], 
        " and phis of ", [value(phi_C), value(phi_P), value(phi_M), value(phi_UT)])

        # Save results
        optimal_values_temp =   [exp.(JuMP.value.(logmu))*60*60, value(S_c_M), value(C_c_M), value(P_c_M), value(lm_M), 
                            value(S_c*volume), value(C_c*volume), value(P_c*volume), value(lm*volume),
                            value(phi_C), value(phi_P), value(phi_M), value(phi_UT), value(beta), 
#                            value(UT*volume), value(CE*volume), value(PE*volume), value(ME*volume), 
                            value(v_UT), value(v_C), value(v_P), value(v_M), 
                            value(D), value(area), value(volume), value(r)]
        if optimal_values_temp[1] >= optimal_values[1]
            optimal_values = optimal_values_temp
        end
    end
    push!(df, vcat(value(S_e_M), optimal_values))
end

# Write the results to a CSV file
CSV.write(file_name, df, header=true)
