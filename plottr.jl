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
#df = filter(row -> row[:optmsr_response] in ["LOCALLY_SOLVED", "ALMOST_LOCALLY_SOLVED", "ITERATION_LIMIT"], df)
df = filter(row -> row[:optmsr_response] in ["LOCALLY_SOLVED", "ALMOST_LOCALLY_SOLVED"], df)

# Extract the p values and final growth rates
p_values = df[:, "eP"]
growth_rates = df[:, "mu_per_hour"]

# Add a new column with the log-transformed values
df[!, :log_eP] = log10.(df[!, :eP])

# Define the color palette
my_palette = palette(:Set3_11)

# Plot the growth rate over the environmental concentrations with loess approximation
yticks = 0.0:0.1:0.4
Fig_growth_rate = scatter(df.log_eP, df.mu_per_hour, color = my_palette[1],
    xlabel = "", ylabel = "growth rate (1/h)", yticks = yticks, ylims = (0, 0.45),
    legend = :outerright,
    label = "microbial cell divisions per hour", dpi = 300)
lo_gr=loess(df.log_eP, df.mu_per_hour)
range_all=range(minimum(df.log_eP), maximum(df.log_eP))
pred_lo_gr=predict(lo_gr, range_all)
plot!(Fig_growth_rate, range_all, pred_lo_gr, line = my_palette[1], label = "loess approximation")

#strokecolors = [RGBA{Float64}(1, 1, 1, (alpha/1000-0.1)*5) for alpha in df.density_PERI]
#scatter!(Fig_growth_rate, df.log_eP, fill(0.5, length(df.log_eP)), markersize = 50*df.r_cyt, 
#label = "cytoplasm (red) with periplasm (black)\nand opacity indicating density", markercolor = :transparent,
#markerstrokewidth = df.r_pp*50, markerstrokecolor = strokecolors, alpha = strokecolors)

colors = [RGBA{Float64}(1, 0, 0, (alpha/1000-0.1)*5) for alpha in df.density_CYTO]
scatter!(Fig_growth_rate, df.log_eP, fill(0.35, length(df.log_eP)), markersize = 10* sqrt.(df.r_cyt), 
label = "cytoplasm and periplasm size\nwith opacity indicating density", markerstrokewidth = df.r_pp*50, markercolor = colors, alpha = colors)

# Plot the core proteome fractions over the environmental concentrations with loess approximation
Fig_core_phis = scatter(df.log_eP, df.phi_C, color = my_palette[2],
    xlabel = "", ylabel = "Core proteome fractions",
    ylims = (-0.05, 0.5), legend = :outerright, label = "catabolism", dpi = 300)
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
    xlabel = "", ylabel = "CAZyme fractions",
    ylims = (-0.05, 0.5), legend = :outerright, label = "endo-laminarase", dpi = 300)
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
#vline!(Fig_CAZY_phis, [log10.(2e-6)], linestyle=:dot, linecolor=:grey, label="")

# Plot cytoplasm and periplasm density over the environmental concentrations with loess approximation
yticks_CELLs = 0.0:0.2:1.0
Fig_CELLs = scatter(df.log_eP, df.C_to_ATP, color = my_palette[1], label = "fraction of C used for ATP",
    xlabel = "log10(env. laminarin conc., M)", ylabel = "fraction of total", yticks = yticks_CELLs, 
    ylims = (-0.1, 1.1), legend = :outerright, dpi = 300)
scatter!(df.log_eP, df.eP_surf_rel, color = my_palette[11], label = "fraction of substrate at surface")
scatter!(df.log_eP, df.oligoloss_rate./(df.v_ENDO .* df.vol_exrp .* 1e-15 .* Av), color = my_palette[6], label = "lost oligosaccharide fraction")

# Display all plots together
Figures_PAM = plot(Fig_growth_rate, Fig_core_phis, Fig_CAZY_phis, Fig_CELLs,
    layout = (4, 1), size = (800, 800), dpi = 300)

# Save to file
savefig(Figures_PAM, "output_$(datetime_str).png")
