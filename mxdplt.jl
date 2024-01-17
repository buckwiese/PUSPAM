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

# Add a new column with the log-transformed values
df[!, :log_ePalg] = log10.(df[!, :ePalg])

# Define the color palette
my_palette = palette(:Set3_11)

# Plot the growth rate over the environmental concentrations with loess approximation
yticks = 0.0:0.1:0.4
Fig_growth_rate = scatter(df.log_ePalg, df.mu_per_hour, color = my_palette[1],
    xlabel = "", ylabel = "growth rate (1/h)", yticks = yticks, ylims = (0, 0.25),
    legend = :outerright,
    label = "microbial cell divisions per hour", dpi = 300)
lo_gr=loess(df.log_ePalg, df.mu_per_hour)
range_all=range(minimum(df.log_ePalg), maximum(df.log_ePalg))
pred_lo_gr=predict(lo_gr, range_all)
plot!(Fig_growth_rate, range_all, pred_lo_gr, line = my_palette[1], label = "loess approximation")

#strokecolors = [RGBA{Float64}(1, 1, 1, (alpha/1000-0.1)*5) for alpha in df.density_PERI]
#scatter!(Fig_growth_rate, df.log_eP, fill(0.5, length(df.log_eP)), markersize = 50*df.r_cyt, 
#label = "cytoplasm (red) with periplasm (black)\nand opacity indicating density", markercolor = :transparent,
#markerstrokewidth = df.r_pp*50, markerstrokecolor = strokecolors, alpha = strokecolors)

colors = [RGBA{Float64}(1, 0, 0, (alpha/1000-0.1)*5) for alpha in df.density_CYTO]
scatter!(Fig_growth_rate, df.log_ePalg, fill(0.20, length(df.log_ePalg)), markersize = 10* sqrt.(df.r_cyt), 
label = "cytoplasm and periplasm size\nwith opacity indicating density", markerstrokewidth = df.r_pp*50, markercolor = colors, alpha = colors)

# Plot the core proteome fractions over the environmental concentrations with loess approximation
Fig_core_phis = scatter(df.log_ePalg, df.phi_C, color = my_palette[2],
    xlabel = "", ylabel = "Core proteome fractions",
    ylims = (-0.025, 0.25), legend = :outerright, label = "glycolysis and related enzymes", dpi = 300)
lo_phi_C=loess(df.log_ePalg, df.phi_C)
pred_lo_phi_C=predict(lo_phi_C, range_all)
plot!(Fig_core_phis, range_all, pred_lo_phi_C, line = my_palette[2], label = "")
scatter!(df.log_ePalg, df.phi_P, color = my_palette[3], label = "ribosomal proteins")
lo_phi_P=loess(df.log_ePalg, df.phi_P)
pred_lo_phi_P=predict(lo_phi_P, range_all)
plot!(Fig_core_phis, range_all, pred_lo_phi_P, line = my_palette[3], label = "")
scatter!(df.log_ePalg, df.phi_M, color = my_palette[4], label = "membrane synthesis")
lo_phi_M=loess(df.log_ePalg, df.phi_M)
pred_lo_phi_M=predict(lo_phi_M, range_all)
plot!(Fig_core_phis, range_all, pred_lo_phi_M, line = my_palette[4], label = "")
#vline!(Fig_core_phis, [log10.(2e-6)], linestyle=:dot, linecolor=:grey, label="")

# Plot the laminarin CAZyme proteome fractions over the environmental concentrations with loess approximation
Fig_lamCAZY_phis = scatter(df.log_ePalg, df.phi_GH16, color = my_palette[6],
    xlabel = "", ylabel = "CAZyme fractions",
    ylims = (-0.025, 0.25), legend = :outerright, label = "endo-laminarase", dpi = 300)
lo_phi_GH16=loess(df.log_ePalg, df.phi_GH16)
pred_lo_phi_GH16=predict(lo_phi_GH16, range_all)
plot!(Fig_lamCAZY_phis, range_all, pred_lo_phi_GH16, line = my_palette[6], label = "")
scatter!(df.log_ePalg, df.phi_TBDTlam, color = my_palette[7], label = "Ton-B dependent transporters ")
lo_phi_TBDTlam=loess(df.log_ePalg, df.phi_TBDTlam)
pred_lo_phi_TBDTlam=predict(lo_phi_TBDTlam, range_all)
plot!(Fig_lamCAZY_phis, range_all, pred_lo_phi_TBDTlam, line = my_palette[7], label = "")
scatter!(df.log_ePalg, df.phi_GH30, color = my_palette[8], label = "debranching laminarase")
lo_phi_GH30=loess(df.log_ePalg, df.phi_GH30)
pred_lo_phi_GH30=predict(lo_phi_GH30, range_all)
plot!(Fig_lamCAZY_phis, range_all, pred_lo_phi_GH30, line = my_palette[8], label = "")
scatter!(df.log_ePalg, df.phi_GH17, color = my_palette[9], label = "exo-laminarase")
lo_phi_GH17=loess(df.log_ePalg, df.phi_GH17)
pred_lo_phi_GH17=predict(lo_phi_GH17, range_all)
plot!(Fig_lamCAZY_phis, range_all, pred_lo_phi_GH17, line = my_palette[9], label = "")
scatter!(df.log_ePalg, df.phi_UGlc, color = my_palette[10], label = "glucose uptake")
lo_phi_UGlc=loess(df.log_ePalg, df.phi_UGlc)
pred_lo_phi_UGlc=predict(lo_phi_UGlc, range_all)
plot!(Fig_lamCAZY_phis, range_all, pred_lo_phi_UGlc, line = my_palette[10], label = "")

# Plot the alginate CAZyme proteome fractions over the environmental concentrations with loess approximation
Fig_algCAZY_phis = scatter(df.log_ePalg, df.phi_PL7, color = my_palette[6],
    xlabel = "", ylabel = "CAZyme fractions",
    ylims = (-0.025, 0.25), legend = :outerright, label = "sentinel alginate lyase", dpi = 300)
lo_phi_PL7=loess(df.log_ePalg, df.phi_PL7)
pred_lo_phi_PL7=predict(lo_phi_PL7, range_all)
plot!(Fig_algCAZY_phis, range_all, pred_lo_phi_PL7, line = my_palette[6], label = "")
scatter!(df.log_ePalg, df.phi_PL5, color = my_palette[5], label = "M-specific endo-alginate lyase")
lo_phi_PL5=loess(df.log_ePalg, df.phi_PL5)
pred_lo_phi_PL5=predict(lo_phi_PL5, range_all)
plot!(Fig_algCAZY_phis, range_all, pred_lo_phi_PL5, line = my_palette[5], label = "")
scatter!(df.log_ePalg, df.phi_PL38, color = my_palette[11], label = "G-specific endo-alginate lyase")
lo_phi_PL38=loess(df.log_ePalg, df.phi_PL38)
pred_lo_phi_PL38=predict(lo_phi_PL38, range_all)
plot!(Fig_algCAZY_phis, range_all, pred_lo_phi_PL38, line = my_palette[11], label = "")
scatter!(df.log_ePalg, df.phi_TBDTalg, color = my_palette[7], label = "Ton-B dependent transporters ")
lo_phi_TBDTalg=loess(df.log_ePalg, df.phi_TBDTalg)
pred_lo_phi_TBDTalg=predict(lo_phi_TBDTalg, range_all)
plot!(Fig_algCAZY_phis, range_all, pred_lo_phi_TBDTalg, line = my_palette[7], label = "")
scatter!(df.log_ePalg, df.phi_PL17, color = my_palette[8], label = "M-specific exo-alginate lyase")
lo_phi_PL17=loess(df.log_ePalg, df.phi_PL17)
pred_lo_phi_PL17=predict(lo_phi_PL17, range_all)
plot!(Fig_algCAZY_phis, range_all, pred_lo_phi_PL17, line = my_palette[8], label = "")
scatter!(df.log_ePalg, df.phi_PL6, color = my_palette[9], label = "G-specific exo-alginate lyase")
lo_phi_PL6=loess(df.log_ePalg, df.phi_PL6)
pred_lo_phi_PL6=predict(lo_phi_PL6, range_all)
plot!(Fig_algCAZY_phis, range_all, pred_lo_phi_PL6, line = my_palette[9], label = "")
scatter!(df.log_ePalg, df.phi_UUro, color = my_palette[10], label = "uronate uptake")
lo_phi_UUro=loess(df.log_ePalg, df.phi_UUro)
pred_lo_phi_UUro=predict(lo_phi_UUro, range_all)
plot!(Fig_algCAZY_phis, range_all, pred_lo_phi_UUro, line = my_palette[10], label = "")
#vline!(Fig_algCAZY_phis, [log10.(2e-6)], linestyle=:dot, linecolor=:grey, label="")

# Plot cytoplasm and periplasm density over the environmental concentrations with loess approximation
yticks_CELLs = 0.0:0.2:1.0
Fig_CELLs = scatter(df.log_ePalg, df.C_to_ATP, color = my_palette[1], label = "fraction of C used for ATP",
    xlabel = "log10(env. alginate conc., M)", ylabel = "fraction of total", yticks = yticks_CELLs, 
    ylims = (-0.1, 1.1), legend = :outerright, dpi = 300)
scatter!(df.log_ePalg, df.ePlam_surf_rel, color = my_palette[10], label = "fraction of laminarin at surface")
scatter!(df.log_ePalg, df.ePalg_surf_rel, color = my_palette[11], label = "fraction of alginate at surface")
#scatter!(df.log_ePalg, df.oligoloss_rate_lam./maximum(df.oligoloss_rate_lam), color = my_palette[6], label = "lost laminarin oligosaccharides")
#scatter!(df.log_ePalg, df.oligoloss_rate_LOalg./maximum(df.oligoloss_rate_LOalg), color = my_palette[7], label = "lost large alginate oligos")
#scatter!(df.log_ePalg, df.oligoloss_rate_SOalg./maximum(df.oligoloss_rate_SOalg), color = my_palette[8], label = "lost small alginate oligos")

# Display all plots together
title!(Fig_growth_rate, "with laminarin concentration of $(df[1, :ePlam]) (M)", titlefont = font(12))
Figures_PAM = plot(Fig_growth_rate, Fig_core_phis, Fig_lamCAZY_phis,  Fig_algCAZY_phis, Fig_CELLs,
    layout = (5, 1), size = (1000, 800), dpi = 300)

# Save to file
output_file_name = replace(file_name, "results_" => "output_")
output_file_name = replace(output_file_name, ".csv" => ".png")
savefig(Figures_PAM, output_file_name)
