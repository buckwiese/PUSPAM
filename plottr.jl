using Loess
using StatsPlots
using Cairo
using Fontconfig
using CSV
using DataFrames
using Plots
using ColorSchemes
using StatsBase

# Read the results from the CSV file
df = CSV.read(file_name, DataFrame)
#df = filter(row -> row[:optmsr_response] in ["LOCALLY_SOLVED", "ALMOST_LOCALLY_SOLVED"], df)

# Extract the p values and final growth rates
p_values = df[:, "eP"]
growth_rates = df[:, "mu_per_hour"]

# Add a new column with the log-transformed values
df[!, :log_eP] = log10.(df[!, :eP])

# Define the color palette
my_palette = palette(:Set3_11)

# Plot the growth rate over the environmental concentrations with loess approximation
yticks = 0.0:0.1:0.6
Fig_growth_rate = scatter(df.log_eP, df.mu_per_hour, color = my_palette[1],
    xlabel = "", ylabel = "growth rate (1/h)", yticks = yticks, ylims = (0, 0),
    legend = :outerright,
    label = "microbial cell divisions per hour")
lo_gr=loess(df.log_eP, df.mu_per_hour)
range_all=range(minimum(df.log_eP), maximum(df.log_eP))
pred_lo_gr=predict(lo_gr, range_all)
plot!(Fig_growth_rate, range_all, pred_lo_gr, line = my_palette[1], label = "loess approximation")

# Add a circle representing r_cy = 0.2
scatter!(df.log_eP, df.mu_per_hour .+ 0.2, color = :red, markersize = 50*df.r_cyt, label = "cytoplasm with periplasm", markerstrokewidth = df.r_pp*50)

# Plot the core proteome fractions over the environmental concentrations with loess approximation
Fig_core_phis = scatter(df.log_eP, df.phi_C, color = my_palette[2],
    xlabel = "", ylabel = "Core proteome fractions",
    ylims = (-0.05, 0.5), legend = :outerright, label = "catabolism")
lo_phi_C=loess(df.log_eP, df.phi_C)
pred_lo_phi_C=predict(lo_phi_C, range_all)
plot!(Fig_core_phis, range_all, pred_lo_phi_C, line = my_palette[2], label = "")
scatter!(df.log_eP, df.phi_P, color = my_palette[3], label = "ribosomes")
lo_phi_P=loess(df.log_eP, df.phi_P)
pred_lo_phi_P=predict(lo_phi_P, range_all)
plot!(Fig_core_phis, range_all, pred_lo_phi_P, line = my_palette[3], label = "")
scatter!(df.log_eP, df.phi_M, color = my_palette[4], label = "membrane synthesis")
lo_phi_M=loess(df.log_eP, df.phi_M)
pred_lo_phi_M=predict(lo_phi_M, range_all)
plot!(Fig_core_phis, range_all, pred_lo_phi_M, line = my_palette[4], label = "")
scatter!(df.log_eP, df.phi_UT, color = my_palette[5], label = "inner membrane transporters ")
lo_phi_UT=loess(df.log_eP, df.phi_UT)
pred_lo_phi_UT=predict(lo_phi_UT, range_all)
plot!(Fig_core_phis, range_all, pred_lo_phi_UT, line = my_palette[5], label = "")
#vline!(Fig_core_phis, [log10.(2e-6)], linestyle=:dot, linecolor=:grey, label="")

# Plot the CAZyme proteome fractions over the environmental concentrations with loess approximation
Fig_CAZY_phis = scatter(df.log_eP, df.phi_GH16, color = my_palette[6],
    xlabel = "log10(env. laminarin conc., M)", ylabel = "CAZyme fractions",
    ylims = (-0.05, 0.5), legend = :outerright, label = "endo-laminarase")
lo_phi_GH16=loess(df.log_eP, df.phi_GH16)
pred_lo_phi_GH16=predict(lo_phi_GH16, range_all)
plot!(Fig_CAZY_phis, range_all, pred_lo_phi_GH16, line = my_palette[6], label = "")
scatter!(df.log_eP, df.phi_TBDT, color = my_palette[7], label = "Ton-B dependent transporters ")
lo_phi_TBDT=loess(df.log_eP, df.phi_TBDT)
pred_lo_phi_TBDT=predict(lo_phi_TBDT, range_all)
plot!(Fig_CAZY_phis, range_all, pred_lo_phi_TBDT, line = my_palette[7], label = "")
scatter!(df.log_eP, df.phi_GH30, color = my_palette[8], label = "debranching laminarase")
lo_phi_GH30=loess(df.log_eP, df.phi_GH30)
pred_lo_phi_GH30=predict(lo_phi_GH30, range_all)
plot!(Fig_CAZY_phis, range_all, pred_lo_phi_GH30, line = my_palette[8], label = "")
scatter!(df.log_eP, df.phi_GH17, color = my_palette[9], label = "exo-laminarase")
lo_phi_GH17=loess(df.log_eP, df.phi_GH17)
pred_lo_phi_GH17=predict(lo_phi_GH17, range_all)
plot!(Fig_CAZY_phis, range_all, pred_lo_phi_GH17, line = my_palette[9], label = "")
scatter!(df.log_eP, df.phi_porin, color = my_palette[10], label = "outer membrane porin")
lo_phi_porin=loess(df.log_eP, df.phi_porin)
pred_lo_phi_porin=predict(lo_phi_porin, range_all)
plot!(Fig_CAZY_phis, range_all, pred_lo_phi_porin, line = my_palette[10], label = "")
#vline!(Fig_CAZY_phis, [log10.(2e-6)], linestyle=:dot, linecolor=:grey, label="")

# Display all plots together
Figures_PAM = plot(Fig_growth_rate, Fig_core_phis, Fig_CAZY_phis,
    layout = (3, 1), size = (800, 600))

# Save to file
savefig(Figures_PAM, "output_$(datetime_str).png")
