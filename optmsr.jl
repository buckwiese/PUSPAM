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
df = DataFrame(S_e = Any[], logmu = Float64[], S_c = Float64[], C_c = Float64[], P_c = Float64[], lm = Float64[], phi_C = Float64[], phi_P = Float64[], phi_M = Float64[], phi_UT = Float64[], beta = Float64[], UT = Float64[], CE = Float64[], PE = Float64[], ME = Float64[], v_UT = Float64[], v_C = Float64[], v_P = Float64[], v_M = Float64[], D = Float64[], area = Float64[], volume = Float64[])

# Load values specifying proteome allocation problem
#PAM_parameters=load(PAM_parameter_dict) # this can go, Julia remembers!

start_log = log10(1)
stop_log = log10(100000)
S_e_change = 10 .^ range(start_log, stop = stop_log, length = 20) # μM, external substrate concentration       

for i in 1:length(S_e_change)
    set_value(S_e, S_e_change[i])

    # set up optimizer
    JuMP.optimize!(model)

    # Set the callback function
    println("with ext. substrate concentration ", [value(S_e)], " got growth rate ", [exp.(JuMP.value.(logmu))], 
    " and phis of ", [value(phi_C), value(phi_P), value(phi_M), value(phi_UT)])

    # Save results
    optimal_values = [exp.(JuMP.value.(logmu)), value(S_c), value(C_c), value(P_c), value(lm), value(phi_C), value(phi_P), value(phi_M), value(phi_UT), value(beta), value(UT), value(CE), value(PE), value(ME), value(v_UT), value(v_C), value(v_P), value(v_M), value(D), value(area), value(volume)]
    push!(df, vcat(value(S_e), optimal_values))
end

# Write the results to a CSV file
CSV.write(file_name, df, header=true)
