# visualize laminarin-alginate preferences
using Loess
using StatsPlots
using Cairo
using Fontconfig
using CSV
using DataFrames
using Plots
using ColorSchemes
using StatsBase

my_palette = palette(:Set3_11)

df = CSV.read(file_name, DataFrame)
df = filter(row -> row[:optmsr_response] in ["LOCALLY_SOLVED", "ALMOST_LOCALLY_SOLVED", "ITERATION_LIMIT"], df)
df = filter(row -> row[:optmsr_response] in ["LOCALLY_SOLVED", "ALMOST_LOCALLY_SOLVED"], df)

df[!, :log_ePalg] = log10.(df[!, :ePalg])
df[!, :rel_growth_rate] = df[!, :mu_per_hour] ./ maximum(df[!, :mu_per_hour])
df[!, :rel_lam_invest_enz] = df[!, :phi_GH16] ./ maximum(df[!, :phi_GH16])
df[!, :rel_lam_invest_transp] = df[!, :phi_TBDTlam] ./ maximum(df[!, :phi_TBDTlam])
df[!, :rel_alg_invest_enz] = df[!, :phi_PL17] ./ maximum(df[!, :phi_PL17])
df[!, :rel_alg_invest_transp] = df[!, :phi_TBDTalg] ./ maximum(df[!, :phi_TBDTalg])


yticks = 0:0.2:1.0
Fig_growth_rate = scatter(df.log_ePalg, df.rel_lam_invest_enz, 
    color = my_palette[2], label = "investment in laminarinase", 
    legend = :outertop, legend_columns = 2, ylabel = "relative to max", 
    ylims = (-0.05, 1.1), yticks = 0:0.2:1.0, alpha = 0.7)

scatter!(Fig_growth_rate, df.log_ePalg, df.rel_lam_invest_transp,
    color = my_palette[4], label = "investment in laminarin transporter", 
    legend = :outertop, ylabel = "relative to max", 
    ylims = (-0.05, 1.1), yticks = 0:0.2:1.0, alpha = 0.7)

scatter!(Fig_growth_rate, df.log_ePalg, df.rel_alg_invest_enz,
    color = my_palette[3], label = "investment in alginate lyase", 
    legend = :outertop, ylabel = "relative to max", 
    ylims = (-0.05, 1.1), yticks = 0:0.2:1.0, alpha = 0.7)

scatter!(Fig_growth_rate, df.log_ePalg, df.rel_alg_invest_transp,
    color = my_palette[5], label = "investment in alginate transporter", 
    legend = :outertop, ylabel = "relative to max", 
    ylims = (-0.05, 1.1), yticks = 0:0.2:1.0, alpha = 0.7)

scatter!(Fig_growth_rate, df.log_ePalg, df.rel_growth_rate, color = my_palette[1],
    xlabel = "Alginate concentration (M, log10)", ylabel = "relative to max", yticks = yticks, ylims = (-0.05, 1.05),
    legend = :outertop, legend_columns = 2,
    label = "divisions per hour", dpi = 300, alpha = 0.7)

title!(Fig_growth_rate, "with background laminarin concentration of $(df[1, :ePlam]) (M)", titlefont = font(12))

plot(Fig_growth_rate, size = (800, 600), dpi = 300, layout = (1, 1), legend=:outertop, ncol=2)

savefig(Fig_growth_rate, "Fig_mixed_model_laminarin$(df[1, :ePlam])M.png")