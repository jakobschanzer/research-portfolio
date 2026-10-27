#!/usr/bin/env julia

ENV["GKSwstype"] = "100" # headless GR output; avoids opening a GUI workstation

using XLSX
using Plots

const R_KJ = 8.31446261815324e-3
const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const MAPS_FILE = joinpath(PROJECT_ROOT, "data", "processed", "spatial_maps.xlsx")
const OUTPUT_DIR = joinpath(PROJECT_ROOT, "figures")

# Nature-compatible compact sizing. Helvetica falls back to Arial/sans-serif
# if it is unavailable on a given system.
default(
    fontfamily = "Helvetica",
    guidefontsize = 14,
    tickfontsize = 12,
    titlefontsize = 15,
    annotationfontsize = 13,
    linewidth = 1.0,
    foreground_color_subplot = :black,
    background_color = :white,
)

function numeric_block(sheet, rows, columns)
    data = Matrix{Float64}(undef, length(rows), length(columns))
    for (i, row) in enumerate(rows), (j, column) in enumerate(columns)
        value = sheet[row, column]
        data[i, j] = value isa Number ? Float64(value) : NaN
    end
    data
end

function bilinear_upsample(x, y, values; factor = 4)
    """Smooth display interpolation only; no new physical resolution is implied."""
    nx = (length(x) - 1) * factor + 1
    ny = (length(y) - 1) * factor + 1
    new_x = range(first(x), last(x), length = nx)
    new_y = range(first(y), last(y), length = ny)
    output = Matrix{Float64}(undef, ny, nx)
    for (iy, _) in enumerate(new_y), (ix, _) in enumerate(new_x)
        source_x = (ix - 1) / factor + 1
        source_y = (iy - 1) / factor + 1
        x0 = clamp(floor(Int, source_x), 1, length(x))
        y0 = clamp(floor(Int, source_y), 1, length(y))
        x1 = min(x0 + 1, length(x))
        y1 = min(y0 + 1, length(y))
        wx = source_x - x0
        wy = source_y - y0
        corners = (values[y0, x0], values[y0, x1], values[y1, x0], values[y1, x1])
        if any(isnan, corners)
            output[iy, ix] = values[round(Int, source_y), round(Int, source_x)]
        else
            output[iy, ix] =
                (1 - wx) * (1 - wy) * corners[1] +
                wx * (1 - wy) * corners[2] +
                (1 - wx) * wy * corners[3] +
                wx * wy * corners[4]
        end
    end
    collect(new_x), collect(new_y), output
end

function host_maps(sheet)
    r = vec(numeric_block(sheet, 3:3, 3:110))

    # Rows in Excel run from positive z to negative z. Reverse both the z axis
    # and each data matrix so plotting uses the conventional bottom-to-top axis.
    z_density = vec(numeric_block(sheet, 4:111, 2:2))
    z_hbonds = vec(numeric_block(sheet, 115:222, 2:2))
    z_enthalpy = vec(numeric_block(sheet, 227:334, 2:2))
    density = reverse(numeric_block(sheet, 4:111, 3:110), dims = 1)
    hbonds = reverse(numeric_block(sheet, 115:222, 3:110), dims = 1)
    h_over_r = reverse(numeric_block(sheet, 227:334, 3:110), dims = 1)

    # Hank's stored coefficient is H/R in kelvin. The sentinel ±22222.222...
    # denotes undefined/excluded locations and should not be color-mapped.
    enthalpy = h_over_r .* R_KJ
    enthalpy[abs.(h_over_r) .> 20_000.0] .= NaN

    return (
        r = r,
        z_density = reverse(z_density),
        z_hbonds = reverse(z_hbonds),
        z_enthalpy = reverse(z_enthalpy),
        density = density,
        hbonds = hbonds,
        enthalpy = enthalpy,
    )
end

function map_panel(x, y, values; limits, palette, title = "", ylabel = "",
                   xlabel = "", colorbar = false, colorbar_title = "",
                   colorbar_ticks = nothing, letter = "")
    kwargs = (
        aspect_ratio = :equal,
        xlims = (-13.5, 13.5),
        ylims = (-13.5, 13.5),
        xticks = [-10, 0, 10],
        yticks = [-10, 0, 10],
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        clims = limits,
        color = palette,
        colorbar = colorbar,
        colorbar_title = colorbar_title,
        framestyle = :box,
        grid = false,
        widen = false,
        left_margin = isempty(ylabel) ? 1.5Plots.mm : 5Plots.mm,
        right_margin = 1.5Plots.mm,
        top_margin = 1.5Plots.mm,
        bottom_margin = 1.5Plots.mm,
    )
    panel = heatmap(x, y, values; kwargs...)
    if colorbar && colorbar_ticks !== nothing
        plot!(panel, colorbar_ticks = colorbar_ticks)
    end
    annotate!(panel, -12.6, 12.4, text(letter, 13, :left, :top, :black))
    panel
end

workbook = XLSX.readxlsx(MAPS_FILE)
hosts = [
    (sheet = "CB7", title = "CB7"),
    (sheet = "OA", title = "OA"),
    (sheet = "beta", title = "restrained β-CD"),
]
maps = [host_maps(workbook[item.sheet]) for item in hosts]
# Bulk reference values at 25 C extracted from the manuscript analysis
# workbooks. Keeping the three scalars here avoids distributing large raw
# trajectory-derived workbooks solely to recover one reference cell each.
bulk_hbonds = [
    3.63633882228599, # CB7
    3.61521237780236, # OA
    3.63661387156408, # restrained beta-CD
]

# Orient restrained beta-CD with its wider rim at positive z (top of panels).
beta = maps[3]
maps[3] = merge(beta, (
    density = reverse(beta.density, dims = 1),
    hbonds = reverse(beta.hbonds, dims = 1),
    enthalpy = reverse(beta.enthalpy, dims = 1),
))

# Sequential grayscale: white means zero/depleted; black means larger value.
sequential_gray = cgrad(:grays, rev = true)
# H-bond deviation from bulk: black is deficient, middle gray is bulk-like,
# and white is enhanced. Zero lies at 75% of the asymmetric -1.5 to +0.5 scale.
hbond_gray = cgrad(
    ["#101010", "#808080", "#f0f0f0"],
    [0.0, 0.75, 1.0],
)
# Signed monotonic grayscale: black is negative, mid-gray is zero, white is positive.
signed_gray = cgrad(:grays)

panels = Plots.Plot[]
letters = collect('a':'i')

# H-bond maps are shown first. Bilinear interpolation improves display
# smoothness while preserving the original numerical grid and scale.
for (column, item) in enumerate(hosts)
    m = maps[column]
    hb_delta = m.hbonds .- bulk_hbonds[column]
    hb_delta[m.density .< 0.05] .= NaN
    hb_r, hb_z, hb_values = bilinear_upsample(m.r, m.z_hbonds, hb_delta)
    panel = map_panel(
        hb_r, hb_z, hb_values;
        limits = (-1.5, 0.5), palette = hbond_gray,
        title = item.title,
        ylabel = column == 1 ? "z (Å)" : "",
        colorbar = column == 3,
        colorbar_title = "Δ H-bonds",
        colorbar_ticks = [-1.5, -1.0, -0.5, 0.0, 0.5],
        letter = string(letters[column]),
    )
    # Mark the boundary below which per-water H-bond coordination is not
    # reported because essentially no water samples the bin.
    contour!(panel, m.r, m.z_density, m.density;
             levels = [0.05], color = :black, linewidth = 0.8,
             label = false)
    push!(panels, panel)
end

for column in 1:3
    m = maps[column]
    push!(panels, map_panel(
        m.r, m.z_density, m.density;
        limits = (0.0, 1.3), palette = sequential_gray,
        ylabel = column == 1 ? "z (Å)" : "",
        colorbar = column == 3,
        colorbar_title = "g(r,z)",
        colorbar_ticks = [0.0, 0.5, 1.0, 1.3],
        letter = string(letters[3 + column]),
    ))
end

for column in 1:3
    m = maps[column]
    push!(panels, map_panel(
        m.r, m.z_enthalpy, m.enthalpy;
        limits = (-4.0, 4.0), palette = signed_gray,
        ylabel = column == 1 ? "z (Å)" : "",
        xlabel = "r (Å)",
        colorbar = column == 3,
        colorbar_title = "ΔH (kJ mol^-1)",
        colorbar_ticks = [-4, -2, 0, 2, 4],
        letter = string(letters[6 + column]),
    ))
end

figure = plot(
    panels...,
    # The third column is wider so its shared colorbar does not reduce the map
    # area relative to the first two columns.
    layout = grid(3, 3, widths = [0.29, 0.29, 0.42]),
    size = (1600, 1350),
    dpi = 600,
    link = :both,
    margin = 1.5Plots.mm,
)

savefig(figure, joinpath(OUTPUT_DIR, "spatial_maps.pdf"))
savefig(figure, joinpath(OUTPUT_DIR, "spatial_maps.svg"))
savefig(figure, joinpath(OUTPUT_DIR, "spatial_maps.png"))

println("Wrote spatial maps to: ", OUTPUT_DIR)
