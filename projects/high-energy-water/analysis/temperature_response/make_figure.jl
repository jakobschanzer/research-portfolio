#!/usr/bin/env julia

ENV["GKSwstype"] = "100"
using XLSX
using Plots

const HERE = @__DIR__
const PROJECT_ROOT = normpath(joinpath(HERE, "..", ".."))
const MEAN_INPUT = joinpath(PROJECT_ROOT, "data", "processed", "unbiased_occupancy.xlsx")
const UNCERTAINTY_INPUT = joinpath(PROJECT_ROOT, "data", "derived", "water_occupancy_by_host.xlsx")
const OUTPUT_DIR = joinpath(PROJECT_ROOT, "figures")
const HOSTS = [
    ("CB7", "CB7", 7, :circle, :solid, 0.05),
    ("OA", "OA", 3, :utriangle, :dash, 0.30),
    ("restrained β-CD", "restrained beta-CD", 5, :diamond, :dot, 0.52),
    ("β-CD", "beta-CD", 9, :square, :dashdot, 0.72),
]

processed_values(sheet, column) = Float64[sheet[row, column] for row in 4:13]
uncertainty_values(sheet, column) = Union{Missing,Float64}[
    sheet[row, column] for row in 2:11
]

data = Dict{String,NamedTuple}()
XLSX.openxlsx(MEAN_INPUT) do mean_workbook
    mean_sheet = mean_workbook[1]
    XLSX.openxlsx(UNCERTAINTY_INPUT) do uncertainty_workbook
        for (host, sheet_name, mean_column, _, _, _) in HOSTS
            uncertainty_sheet = uncertainty_workbook[sheet_name]
            data[host] = (
                temperature = processed_values(mean_sheet, 2),
                mean = processed_values(mean_sheet, mean_column),
                sem = uncertainty_values(uncertainty_sheet, 4),
            )
        end
    end
end

function show_error_bar(host, temperature)
    # The revised unrestrained β-CD 25 °C run has an uncertainty smaller
    # than its marker. Do not reuse the obsolete SEM from the prior run.
    return !(host == "β-CD" && temperature == 25.0)
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

for (host, _, _, marker, style, gray) in HOSTS
    series = data[host]
    color = RGB(gray, gray, gray)
    # Draw capped error bars explicitly because most SEMs are smaller than
    # the markers and would otherwise be visually hidden.
    for (temperature, mean, sem) in zip(series.temperature, series.mean, series.sem)
        show_error_bar(host, temperature) || continue
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
