#!/usr/bin/env julia

ENV["GKSwstype"] = "100"
using XLSX
using Plots

const HERE = @__DIR__
const PROJECT_ROOT = normpath(joinpath(HERE, "..", ".."))
const DELTA_FILE = joinpath(PROJECT_ROOT, "data", "processed", "hydrogen_bond_temperature.xlsx")
const OCC_FILE = joinpath(PROJECT_ROOT, "data", "derived", "water_occupancy_by_host.xlsx")
const OUTPUT_DIR = joinpath(PROJECT_ROOT, "figures")

# Order and redundant encodings are shared between panels so the figure remains
# legible in grayscale, photocopying, and common forms of color-vision deficiency.
const SERIES = [
    ("CB7",               :circle,    :solid,   0.05, 5, 6, "CB7"),
    ("OA",                :utriangle, :dash,    0.28, 3, 4, "OA"),
    ("β-CD",              :square,    :dashdot, 0.68, 7, 8, "beta-CD"),
    ("restrained β-CD",   :diamond,   :dot,     0.48, 9,10, "restrained beta-CD"),
]

temperature = Float64[]
delta = Dict{String,Vector{Float64}}()
delta_se = Dict{String,Vector{Float64}}()
XLSX.openxlsx(DELTA_FILE) do xf
    sh = xf["Summar"]
    append!(temperature, Float64[sh[r,2] for r in 5:14])
    for (name, _, _, _, value_col, error_col, _) in SERIES
        delta[name] = Float64[sh[r,value_col] for r in 5:14]
        delta_se[name] = Float64[sh[r,error_col] for r in 5:14]
    end
end

occupancy = Dict{String,Vector{Float64}}()
occupancy_se = Dict{String,Vector{Float64}}()
XLSX.openxlsx(OCC_FILE) do xf
    for (name, _, _, _, _, _, sheetname) in SERIES
        sh = xf[sheetname]
        occupancy[name] = Float64[sh[r,6] for r in 2:11]    # established means
        occupancy_se[name] = Float64[sh[r,4] for r in 2:11] # trajectory block SEM
    end
end

default(
    fontfamily="Helvetica",
    guidefontsize=17,
    tickfontsize=13,
    legendfontsize=12,
    foreground_color_subplot=:black,
    background_color=:white,
    linewidth=1.8,
)

common = (
    xlims=(0,100), xticks=5:10:95,
    xlabel="Temperature (°C)",
    framestyle=:box, grid=false,
    minorgrid=false, minorticks=false,
    margin=6Plots.mm,
)

p_delta = plot(; common...,
    ylabel="Hydrogen-bond balance, δ",
    ylims=(-3.2,3.5), yticks=-3:1:3,
    legend=:bottomleft,
)
hline!(p_delta, [0.0]; color=RGB(0.65,0.65,0.65), linewidth=1.0, label="")

p_occ = plot(; common...,
    ylabel="Mean cavity occupancy, ⟨N⟩",
    ylims=(2.5,7.0), yticks=2.5:0.5:7.0,
    legend=false,
)

for (name, marker, style, gray, _, _, _) in SERIES
    color = RGB(gray,gray,gray)
    plot!(p_delta, temperature, delta[name];
        yerror=delta_se[name], label=name, color=color,
        linestyle=style, marker=marker, markersize=5.2,
        markercolor=:white, markerstrokecolor=:black,
        markerstrokewidth=1.0, linewidth=1.8,
        errorbar_color=color, errorbar_linewidth=1.1,
    )
    plot!(p_occ, temperature, occupancy[name];
        yerror=occupancy_se[name], label="", color=color,
        linestyle=style, marker=marker, markersize=5.2,
        markercolor=:white, markerstrokecolor=:black,
        markerstrokewidth=1.0, linewidth=1.8,
        errorbar_color=color, errorbar_linewidth=1.1,
    )
end

annotate!(p_delta, 2.2, 3.27, text("a", 18, :bold, :left))
annotate!(p_occ,   2.2, 6.86, text("b", 18, :bold, :left))

fig = plot(p_delta, p_occ;
    layout=(1,2), size=(1600,680), dpi=600,
    left_margin=11Plots.mm, right_margin=6Plots.mm,
    bottom_margin=10Plots.mm, top_margin=5Plots.mm,
)

for ext in ("png", "pdf", "svg")
    savefig(fig, joinpath(OUTPUT_DIR, "temperature_response.$ext"))
end

println("Wrote figures/temperature_response.png, .pdf, and .svg")
