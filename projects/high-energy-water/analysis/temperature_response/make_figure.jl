#!/usr/bin/env julia

ENV["GKSwstype"] = "100"
using XLSX
using Plots

const HERE = @__DIR__
const PROJECT_ROOT = normpath(joinpath(HERE, "..", ".."))
const INPUT = joinpath(PROJECT_ROOT, "data", "derived", "water_occupancy_by_host.xlsx")
const OUTPUT_DIR = joinpath(PROJECT_ROOT, "figures")
const HOSTS = [
    ("CB7", "CB7", :circle, :solid, 0.05),
    ("OA", "OA", :utriangle, :dash, 0.30),
    ("restrained β-CD", "restrained beta-CD", :diamond, :dot, 0.52),
    ("β-CD", "beta-CD", :square, :dashdot, 0.72),
]

values(sheet, column) = Float64[sheet[row, column] for row in 2:11]

data = Dict{String,NamedTuple}()
XLSX.openxlsx(INPUT) do workbook
    for (host, sheet_name, _, _, _) in HOSTS
        sheet = workbook[sheet_name]
        data[host] = (
            temperature = values(sheet, 1),
            mean = values(sheet, 6),
            sem = values(sheet, 4),
        )
    end
end

default(
    fontfamily = "Helvetica",
    guidefontsize = 15,
    tickfontsize = 12,
    legendfontsize = 11,
    linewidth = 1.6,
    foreground_color_subplot = :black,
    background_color = :white,
)

figure = plot(
    xlabel = "Temperature (°C)",
    ylabel = "Mean cavity occupancy, ⟨N⟩",
    xlims = (0, 100),
    ylims = (2.5, 7.0),
    xticks = 5:10:95,
    yticks = 2.5:0.5:7.0,
    legend = :topright,
    framestyle = :box,
    grid = false,
    size = (920, 650),
    dpi = 600,
    margin = 7Plots.mm,
)

for (host, _, marker, style, gray) in HOSTS
    series = data[host]
    color = RGB(gray, gray, gray)
    # Draw capped error bars explicitly because most SEMs are smaller than
    # the markers and would otherwise be visually hidden.
    for (temperature, mean, sem) in zip(series.temperature, series.mean, series.sem)
        halfcap = 0.75
        plot!(figure, [temperature, temperature], [mean-sem, mean+sem];
            color=color, linewidth=1.25, label="")
        plot!(figure, [temperature-halfcap, temperature+halfcap], [mean-sem, mean-sem];
            color=color, linewidth=1.25, label="")
        plot!(figure, [temperature-halfcap, temperature+halfcap], [mean+sem, mean+sem];
            color=color, linewidth=1.25, label="")
    end
    plot!(figure, series.temperature, series.mean;
        label=host, color=color, linestyle=style,
        marker=marker, markersize=5.5,
        markercolor=:white, markerstrokecolor=:black,
        markerstrokewidth=1.1, linewidth=1.7)
end

for extension in ("png", "pdf", "svg")
    savefig(figure, joinpath(OUTPUT_DIR, "figure3_temperature_response.$extension"))
end
println("Wrote figures/figure3_temperature_response.png, .pdf, and .svg")
