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
df = DataFrame(eP_M =  Any[], mu_per_hour = Float64[],
               phi_C = Float64[], phi_P = Float64[], phi_M = Float64[], phi_UT = Float64[],
               phi_GH16 = Float64[], phi_TBDT = Float64[], phi_GH30 = Float64[], phi_GH17=Float64[],
               density = Float64[], r_cyt = Float64[], r_pp = Float64[], vol_cy = Float64[], vol_pp = Float64[],
               eO = Float64[], pO = Float64[], pdbO = Float64[], S_e = Float64[],
               S_c = Float64[], C_c = Float64[], P_c = Float64[], lm = Float64[],
               UT = Float64[], C = Float64[], P = Float64[], M = Float64[],
               GH16 = Float64[], TBDT = Float64[], GH30 = Float64[], GH17 = Float64[])                

# Load values specifying proteome allocation problem
#PAM_parameters=load(PAM_parameter_dict) # this can go, Julia remembers!

start_log = log10(0.00000001)
stop_log = log10(0.001)
eP_change = 10 .^ range(start_log, stop = stop_log, length = 30) # M, external substrate concentration       

for i in 1:length(eP_change)
    set_value(eP_M, eP_change[i])
    optimal_values = zeros(1, 9)
    for j in 1:1
        # set up optimizer
        JuMP.optimize!(model)

        # Set the callback function
        println("with ext. substrate concentration ", [value(eP_M)], " got growth rate ", [exp.(JuMP.value.(logmu))], 
        " and phis of ", [value(phi_C), value(phi_P), value(phi_M), value(phi_UT)])

        # Save results
        optimal_values_temp =[exp.(JuMP.value.(logmu))*60*60, value(phi_C), value(phi_P), value(phi_M), value(phi_UT),
        value(phi_GH16), value(phi_TBDT), value(phi_GH30), value(phi_GH17),
        value(D), value(r_cy), value(r_pp), value(volume_cy), value(volume_pp),
        value(eO), value(pO), value(pdbO), value(S_e),
        value(S_c), value(C_c), value(P_c), value(lm),
        value(UT), value(CE), value(PE), value(ME),
        value(GH16), value(TBDT), value(GH30), value(GH17)]
        if optimal_values_temp[1] >= optimal_values[1]
            optimal_values = optimal_values_temp
        end
    end
    push!(df, vcat(value(eP_M), optimal_values))
end

# Write the results to a CSV file
CSV.write(file_name, df, header=true)
