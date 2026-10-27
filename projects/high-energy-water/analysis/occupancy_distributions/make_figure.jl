#!/usr/bin/env julia

using XLSX
using Base64
using Statistics
using Printf

const HERE = @__DIR__
const PROJECT_ROOT = normpath(joinpath(HERE, "..", ".."))
const INPUT = joinpath(PROJECT_ROOT, "data", "processed", "occupancy_along_binding_coordinate.xlsx")
const SVG_OUT = joinpath(PROJECT_ROOT, "figures", "figure6_occupancy_distributions.svg")
const AUDIT_OUT = joinpath(PROJECT_ROOT, "data", "derived", "pn_fattail_fits.xlsx")

struct OccupancyData
    z::Vector{Float64}
    n::Vector{Float64}
    p::Matrix{Float64}
end

struct FitResult
    coefficients::Vector{Float64}
    fit_indices::Vector{Int}
    mode_index::Int
    r_squared::Float64
    delta_ln_p0::Float64
    has_wet_basin::Bool
end

xml_escape(value) = replace(string(value), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;")
parse_z(value) = parse(Float64, replace(string(value), "z=" => ""))

function read_host(sheet; header_row, first_col, last_col, n_rows, explicit_n_col=nothing)
    z = Float64[parse_z(sheet[header_row, column]) for column in first_col:last_col]
    n = explicit_n_col === nothing ? collect(0.0:n_rows-1) :
        Float64[sheet[header_row + row, explicit_n_col] for row in 1:n_rows]
    p = zeros(n_rows, length(z))
    for (j, column) in enumerate(first_col:last_col), row in 1:n_rows
        value = sheet[header_row + row, column]
        p[row, j] = ismissing(value) ? 0.0 : Float64(value)
    end
    p ./= sum(p, dims=1)
    OccupancyData(z, n, p)
end

function probability_mean(n, p)
    normalized = p ./ sum(p)
    sum(n .* normalized)
end

function read_oa(sheet)
    curve_z = Float64[sheet[row, 1] for row in 4:37]
    curve_mean = Float64[sheet[row, 2] for row in 4:37]
    order = sortperm(curve_z)
    curve_z, curve_mean = curve_z[order], curve_mean[order]

    n = Float64[sheet[row, 7] for row in 4:13]
    probability_columns = Vector{Vector{Float64}}()
    for column in 8:39
        values = [sheet[row, column] for row in 4:13]
        if all(!ismissing, values) && sum(Float64.(values)) > 0
            push!(probability_columns, Float64.(values))
        end
    end
    p = hcat(probability_columns...)
    p ./= sum(p, dims=1)
    derived_mean = Float64[probability_mean(n, p[:, j]) for j in axes(p, 2)]
    rms = Float64[]
    for first in 1:(length(curve_z)-length(derived_mean)+1)
        difference = derived_mean .- curve_mean[first:first+length(derived_mean)-1]
        push!(rms, sqrt(mean(difference .^ 2)))
    end
    first = argmin(rms)
    actual_z = curve_z[first:first+length(derived_mean)-1]
    maximum(abs.(derived_mean .- curve_mean[first:first+length(derived_mean)-1])) < 5e-4 ||
        error("OA P(N|z) columns could not be aligned")
    OccupancyData(actual_z, n, p)
end

cb7 = nothing; oa = nothing; cd = nothing
XLSX.openxlsx(INPUT) do workbook
    global cb7 = read_host(workbook["CB7"];
        header_row=2, first_col=8, last_col=33, n_rows=10)
    global oa = read_oa(workbook["OA"])
    global cd = read_host(workbook["CD_Restrained"];
        header_row=3, first_col=9, last_col=59, n_rows=10, explicit_n_col=8)
end

function selected_distribution(data, z_value)
    index = argmin(abs.(data.z .- z_value))
    abs(data.z[index] - z_value) < 1e-8 || error("No distribution at z=$z_value")
    data.p[:, index] ./ sum(data.p[:, index])
end

function quadratic_fit(n, probability)
    local_maxima = Int[]
    maximum_probability = maximum(probability)
    for i in eachindex(probability)
        left = i == firstindex(probability) ? -Inf : probability[i - 1]
        right = i == lastindex(probability) ? -Inf : probability[i + 1]
        if probability[i] >= left && probability[i] >= right &&
                probability[i] >= 0.10maximum_probability && probability[i] > 0
            push!(local_maxima, i)
        end
    end
    wet_candidates = [i for i in local_maxima if n[i] >= 2]
    mode_index = isempty(wet_candidates) ? argmax(probability) : maximum(wet_candidates)

    if n[mode_index] >= 2
        first_fit = max(firstindex(n), mode_index - 1)
        last_fit = min(lastindex(n), mode_index + 2)
    else
        first_fit = firstindex(n)
        last_fit = min(lastindex(n), firstindex(n) + 3)
    end
    fit_indices = [i for i in first_fit:last_fit if probability[i] > 0]
    length(fit_indices) >= 3 || error("Fewer than three positive-probability fit points")
    x = n[fit_indices]
    y = log.(probability[fit_indices])
    design = hcat(x .^ 2, x, ones(length(x)))
    coefficients = design \ y
    predicted = design * coefficients
    total = sum((y .- mean(y)) .^ 2)
    residual = sum((y .- predicted) .^ 2)
    r_squared = total > 0 ? 1 - residual / total : 1.0
    fitted_zero = coefficients[3]
    delta_zero = probability[1] > 0 ? log(probability[1]) - fitted_zero : NaN
    FitResult(coefficients, fit_indices, mode_index, r_squared, delta_zero,
        !isempty(wet_candidates) && coefficients[1] < 0)
end

function image_uri(path)
    "data:image/png;base64," * base64encode(read(path))
end

const HOST_IMAGES = Dict(
    "CB7" => image_uri(joinpath(HERE, "assets", "cb7_host_red.png")),
    "OA" => image_uri(joinpath(HERE, "assets", "oa_host_red.png")),
    "restrained β-CD" => image_uri(joinpath(HERE, "assets", "betacd_host_red.png")),
)

const NODE_W = 294.0
const NODE_GAP = 16.0
const LEFT = 64.0
const TOP = 78.0
const NODE_H = 392.0
const BLOCK_GAP = 62.0
const WIDTH = LEFT + 6NODE_W + 5NODE_GAP + 30
const CB7_Y = TOP
const OA_Y = CB7_Y + NODE_H + BLOCK_GAP
const BETA_NEG_Y = OA_Y + NODE_H + BLOCK_GAP + 38
const BETA_POS_Y = BETA_NEG_Y + NODE_H + 54
const BETA_SHARED_Y = (BETA_NEG_Y + BETA_POS_Y) / 2
const HEIGHT = BETA_POS_Y + NODE_H + 55

node_x(column) = LEFT + (column - 1) * (NODE_W + NODE_GAP)

function line_path(x_values, y_values, map_x, map_y)
    points = [(@sprintf("%.2f", map_x(x)), @sprintf("%.2f", map_y(y))) for (x, y) in zip(x_values, y_values)]
    isempty(points) && return ""
    "M " * join(["$(point[1]) $(point[2])" for point in points], " L ")
end

function chart_svg(parts, x, y, width, height, n, probability, fit, kind, clip_id)
    plot_left = x + 37
    plot_right = x + width - 8
    plot_top = y + 6
    plot_bottom = y + height - 22
    map_x(value) = plot_left + (value / 9) * (plot_right - plot_left)
    if kind == :probability
        ymin, ymax = 0.0, 1.0
        values = probability
        ylabel = "P(N|z)"
        yticks = [(0.0, "0"), (0.5, "0.5"), (1.0, "1")]
    else
        ymin, ymax = -12.0, 0.0
        values = log.(probability)
        ylabel = "ln P(N|z)"
        yticks = [(-12.0, "−12"), (-6.0, "−6"), (0.0, "0")]
    end
    map_y(value) = plot_bottom - ((value - ymin) / (ymax - ymin)) * (plot_bottom - plot_top)

    push!(parts, "<defs><clipPath id=\"$clip_id\"><rect x=\"$plot_left\" y=\"$plot_top\" width=\"$(plot_right-plot_left)\" height=\"$(plot_bottom-plot_top)\"/></clipPath></defs>")
    push!(parts, "<rect x=\"$plot_left\" y=\"$plot_top\" width=\"$(plot_right-plot_left)\" height=\"$(plot_bottom-plot_top)\" fill=\"white\" stroke=\"#222\" stroke-width=\"0.8\"/>")
    for (tick, label) in yticks
        yy = map_y(tick)
        push!(parts, "<line x1=\"$plot_left\" y1=\"$yy\" x2=\"$plot_right\" y2=\"$yy\" stroke=\"#e2e2e2\" stroke-width=\"0.7\"/>")
        push!(parts, "<text x=\"$(plot_left-5)\" y=\"$(yy+3)\" class=\"tick\" text-anchor=\"end\">$label</text>")
    end
    for tick in (0, 3, 6, 9)
        xx = map_x(tick)
        push!(parts, "<line x1=\"$xx\" y1=\"$plot_bottom\" x2=\"$xx\" y2=\"$(plot_bottom+4)\" stroke=\"#222\" stroke-width=\"0.8\"/>")
        push!(parts, "<text x=\"$xx\" y=\"$(plot_bottom+15)\" class=\"tick\" text-anchor=\"middle\">$tick</text>")
    end
    push!(parts, "<text x=\"$(x+9)\" y=\"$((plot_top+plot_bottom)/2)\" class=\"axis\" text-anchor=\"middle\" transform=\"rotate(-90 $(x+9) $((plot_top+plot_bottom)/2))\">$ylabel</text>")
    kind == :logarithm && push!(parts, "<text x=\"$((plot_left+plot_right)/2)\" y=\"$(y+height-2)\" class=\"axis italic\" text-anchor=\"middle\">N</text>")

    if kind == :probability
        path = line_path(n, values, map_x, map_y)
        push!(parts, "<path d=\"$path\" fill=\"none\" stroke=\"#111\" stroke-width=\"1.4\" clip-path=\"url(#$clip_id)\"/>")
        for i in eachindex(n)
            push!(parts, "<circle cx=\"$(map_x(n[i]))\" cy=\"$(map_y(values[i]))\" r=\"2.7\" fill=\"white\" stroke=\"#111\" stroke-width=\"1.1\" clip-path=\"url(#$clip_id)\"/>")
        end
    else
        valid = [i for i in eachindex(n) if probability[i] > 0 && values[i] >= ymin]
        if !isempty(valid)
            path = line_path(n[valid], values[valid], map_x, map_y)
            push!(parts, "<path d=\"$path\" fill=\"none\" stroke=\"#111\" stroke-width=\"1.4\" clip-path=\"url(#$clip_id)\"/>")
        end
        a, b, c = fit.coefficients
        fitted_set = Set(fit.fit_indices)
        minimum_fit = minimum(fit.fit_indices)
        if fit.has_wet_basin
            dense_n = collect(range(0, 9, length=121))
            fitted = a .* dense_n .^ 2 .+ b .* dense_n .+ c
            fit_path = line_path(dense_n, fitted, map_x, map_y)
            push!(parts, "<path d=\"$fit_path\" fill=\"none\" stroke=\"#D13C2F\" stroke-width=\"1.8\" stroke-dasharray=\"6 4\" clip-path=\"url(#$clip_id)\"/>")
        else
            push!(parts, "<text x=\"$(plot_right-4)\" y=\"$(plot_top+12)\" class=\"tick\" text-anchor=\"end\" fill=\"#666\">wet basin absent</text>")
        end
        for i in valid
            fitted_value = a*n[i]^2 + b*n[i] + c
            fat_tail = fit.has_wet_basin && i < minimum_fit && values[i] - fitted_value > 0.75
            fill = fat_tail ? "#1769AA" : "#111"
            radius = fat_tail ? 3.5 : 2.7
            push!(parts, "<circle cx=\"$(map_x(n[i]))\" cy=\"$(map_y(values[i]))\" r=\"$radius\" fill=\"$fill\" stroke=\"#111\" stroke-width=\"0.7\" clip-path=\"url(#$clip_id)\"/>")
            if fit.has_wet_basin && i in fitted_set
                push!(parts, "<rect x=\"$(map_x(n[i])-4.2)\" y=\"$(map_y(values[i])-4.2)\" width=\"8.4\" height=\"8.4\" fill=\"none\" stroke=\"#D13C2F\" stroke-width=\"1.0\" clip-path=\"url(#$clip_id)\"/>")
            end
        end
    end
end

function guest_svg(cx, cy)
    pieces = String[]
    for angle in range(0, 2pi, length=7)[1:6]
        x = cx + 10cos(angle)
        y = cy + 10sin(angle)
        push!(pieces, "<circle cx=\"$x\" cy=\"$y\" r=\"5.8\" fill=\"#1769AA\"/>")
    end
    "<g>" * join(pieces) * "</g>"
end

function schematic_svg(parts, host, x, y, z_value; show_axis_labels=false)
    # One convention for every host: negative z is left and positive z is
    # right. The guest approaches the centered host from the indicated side.
    portal_on_right = z_value >= 0
    center_x = x + NODE_W * 0.50
    center_y = y + 60
    displacement = (host == "restrained β-CD" ? 100 : 120) * min(abs(z_value), 10) / 10
    direction = portal_on_right ? 1 : -1
    guest_x = abs(z_value) < 1e-10 ? center_x : center_x + direction * displacement
    push!(parts, "<line x1=\"$(x+38)\" y1=\"$center_y\" x2=\"$(x+NODE_W-28)\" y2=\"$center_y\" stroke=\"#777\" stroke-width=\"0.8\" stroke-dasharray=\"3 3\"/>")
    if show_axis_labels
        push!(parts, "<text x=\"$(x+8)\" y=\"$(center_y+16)\" class=\"coordinate\" text-anchor=\"start\">−z</text>")
        push!(parts, "<text x=\"$(x+NODE_W-8)\" y=\"$(center_y+16)\" class=\"coordinate\" text-anchor=\"end\">+z</text>")
    end
    push!(parts, "<image href=\"$(HOST_IMAGES[host])\" x=\"$(center_x-48)\" y=\"$(center_y-34)\" width=\"96\" height=\"68\" preserveAspectRatio=\"xMidYMid meet\" transform=\"rotate(90 $center_x $center_y)\"/>")
    push!(parts, guest_svg(guest_x, center_y))
    if displacement >= 55
        arrow_start = guest_x - direction * 15
        arrow_end = center_x + direction * 52
        push!(parts, "<line x1=\"$arrow_start\" y1=\"$center_y\" x2=\"$arrow_end\" y2=\"$center_y\" stroke=\"#555\" stroke-width=\"1.2\" marker-end=\"url(#arrow)\"/>")
    end
end

function node_svg(parts, records, host, data, column, y, source_z, display_z, tag)
    x = node_x(column)
    probability = selected_distribution(data, source_z)
    fit = quadratic_fit(data.n, probability)
    label = abs(display_z) < 1e-10 ? "z = 0 Å" : @sprintf("z = %+.1f Å", display_z)
    if display_z > 0 && (host == "CB7" || host == "OA")
        label = @sprintf("z = %.1f Å", display_z)
    end
    push!(parts, "<text x=\"$(x+NODE_W/2)\" y=\"$(y+15)\" class=\"distance\" text-anchor=\"middle\">$label</text>")
    schematic_svg(parts, host, x, y+20, display_z; show_axis_labels=column == 1)
    chart_svg(parts, x, y+125, NODE_W, 116, data.n, probability, fit, :probability, "p_$tag")
    chart_svg(parts, x, y+254, NODE_W, 132, data.n, probability, fit, :logarithm, "l_$tag")
    push!(records, (
        host=host, display_z=display_z, source_z=source_z,
        mode_n=data.n[fit.mode_index],
        fit_n=join(string.(Int.(data.n[fit.fit_indices])), ","),
        a=fit.coefficients[1], b=fit.coefficients[2], c=fit.coefficients[3],
        r_squared=fit.r_squared, delta_ln_p0=fit.delta_ln_p0,
        wet_basin_fit_shown=fit.has_wet_basin,
    ))
end

parts = String[
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>",
    "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"$(WIDTH)\" height=\"$(HEIGHT)\" viewBox=\"0 0 $(WIDTH) $(HEIGHT)\">",
    "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>",
    "<style>text{font-family:Arial,Helvetica,sans-serif;fill:#111}.panel{font-size:24px;font-weight:700}.host{font-size:22px;font-weight:600}.branch{font-size:15px;font-weight:600}.distance{font-size:13px;font-weight:600}.axis{font-size:10px}.coordinate{font-size:9px;fill:#555}.tick{font-size:8.5px}.italic{font-style:italic}.legend{font-size:12px}</style>",
    "<defs><marker id=\"arrow\" markerWidth=\"8\" markerHeight=\"8\" refX=\"7\" refY=\"4\" orient=\"auto\"><path d=\"M0,0 L8,4 L0,8 Z\" fill=\"#555\"/></marker></defs>",
]

push!(parts, "<text x=\"18\" y=\"38\" class=\"panel\">a</text><text x=\"$LEFT\" y=\"38\" class=\"host\">CB7</text>")
push!(parts, "<text x=\"$(WIDTH-700)\" y=\"36\" class=\"legend\">observed ln P(N|z)</text><circle cx=\"$(WIDTH-716)\" cy=\"32\" r=\"3\" fill=\"#111\"/>")
push!(parts, "<line x1=\"$(WIDTH-500)\" y1=\"32\" x2=\"$(WIDTH-460)\" y2=\"32\" stroke=\"#D13C2F\" stroke-width=\"2\" stroke-dasharray=\"6 4\"/><text x=\"$(WIDTH-450)\" y=\"36\" class=\"legend\">quadratic wet-basin fit</text>")
push!(parts, "<circle cx=\"$(WIDTH-712)\" cy=\"53\" r=\"4\" fill=\"#1769AA\"/><text x=\"$(WIDTH-700)\" y=\"57\" class=\"legend\">low-N excess above fit</text>")

records = NamedTuple[]
cb7_windows = [10.0, 4.4, 4.0, 3.6, 3.2, 0.0]
oa_sources = [9.6, 4.0, 3.6, 3.2, 2.8, 0.0]
oa_display = [10.0, 4.0, 3.6, 3.2, 2.8, 0.0]
for (column, z) in enumerate(cb7_windows)
    node_svg(parts, records, "CB7", cb7, column, CB7_Y, z, z, "cb7_$column")
end

oa_title_y = OA_Y - 22
push!(parts, "<line x1=\"18\" y1=\"$(oa_title_y-22)\" x2=\"$(WIDTH-20)\" y2=\"$(oa_title_y-22)\" stroke=\"#888\" stroke-width=\"0.8\"/>")
push!(parts, "<text x=\"18\" y=\"$oa_title_y\" class=\"panel\">b</text><text x=\"$LEFT\" y=\"$oa_title_y\" class=\"host\">OA</text>")
for column in 1:6
    node_svg(parts, records, "OA", oa, column, OA_Y, oa_sources[column], oa_display[column], "oa_$column")
end

beta_title_y = BETA_NEG_Y - 42
push!(parts, "<line x1=\"18\" y1=\"$(beta_title_y-25)\" x2=\"$(WIDTH-20)\" y2=\"$(beta_title_y-25)\" stroke=\"#888\" stroke-width=\"0.8\"/>")
push!(parts, "<text x=\"18\" y=\"$beta_title_y\" class=\"panel\">c</text><text x=\"$LEFT\" y=\"$beta_title_y\" class=\"host\">restrained β-CD</text>")
push!(parts, "<text x=\"$(LEFT+210)\" y=\"$(BETA_NEG_Y-13)\" class=\"branch\">negative-z, narrower primary portal</text>")
push!(parts, "<text x=\"$(LEFT+210)\" y=\"$(BETA_POS_Y-13)\" class=\"branch\">positive-z, wider secondary portal</text>")

negative_windows = [-10.0, -4.4, -3.6, -2.8, -2.0]
positive_windows = [10.0, 4.4, 3.6, 2.8, 2.0]
for column in 1:5
    node_svg(parts, records, "restrained β-CD", cd, column, BETA_NEG_Y,
        negative_windows[column], negative_windows[column], "cdn_$column")
    node_svg(parts, records, "restrained β-CD", cd, column, BETA_POS_Y,
        positive_windows[column], positive_windows[column], "cdp_$column")
end
node_svg(parts, records, "restrained β-CD", cd, 6, BETA_SHARED_Y, 0.0, 0.0, "cd0")
shared_x = node_x(6)
source_x = node_x(5) + NODE_W - 5
target_x = shared_x + 8
# node_svg offsets each schematic by 20 px, and schematic_svg places its
# molecular centerline another 60 px below that origin. Anchor the shared
# connector to those true centerlines so the brace is vertically centered.
upper_y = BETA_NEG_Y + 80
lower_y = BETA_POS_Y + 80
center_y = BETA_SHARED_Y + 80
push!(parts, "<path d=\"M $source_x $upper_y C $(source_x+25) $upper_y $(target_x-25) $center_y $target_x $center_y\" fill=\"none\" stroke=\"#555\" stroke-width=\"1.4\" marker-end=\"url(#arrow)\"/>")
push!(parts, "<path d=\"M $source_x $lower_y C $(source_x+25) $lower_y $(target_x-25) $center_y $target_x $center_y\" fill=\"none\" stroke=\"#555\" stroke-width=\"1.4\" marker-end=\"url(#arrow)\"/>")

push!(parts, "</svg>")
write(SVG_OUT, join(parts, "\n"))

XLSX.openxlsx(AUDIT_OUT, mode="w") do workbook
    sheet = workbook[1]
    XLSX.rename!(sheet, "Quadratic fits")
    headers = ["host", "display_z_A", "source_z_A", "wet_mode_N", "fit_N_values",
        "quadratic_a", "quadratic_b", "quadratic_c", "fit_R_squared", "lnP0_minus_fit",
        "wet_basin_fit_shown"]
    sheet["A1:K1"] = reshape(headers, 1, :)
    for (row, record) in enumerate(records)
        values = [record.host, record.display_z, record.source_z, record.mode_n,
            record.fit_n, record.a, record.b, record.c, record.r_squared, record.delta_ln_p0,
            record.wet_basin_fit_shown]
        for (column, value) in enumerate(values)
            sheet[row + 1, column] = value
        end
    end
    notes = XLSX.addsheet!(workbook, "Method")
    notes["A1:B1"] = reshape(["item", "definition"], 1, :)
    method_rows = [
        "fit target" => "ln P(N|z)",
        "wet basin" => "highest-N local maximum at N >= 2 with probability at least 10% of the global maximum and a concave-down quadratic fit",
        "fit window" => "wet-basin mode minus one occupancy through mode plus two occupancies; no fit is displayed when a wet basin is absent",
        "model" => "ln P_fit = a N^2 + b N + c, unweighted least squares",
        "fat-tail marker" => "low-N observed point at least 0.75 natural-log units above the extrapolated quadratic fit",
        "OA z=10 display" => "uses the equivalent separated-state P(N|z) distribution sampled at z=9.6 A, as previously authorized",
    ]
    for (row, pair) in enumerate(method_rows)
        notes[row + 1, 1] = first(pair)
        notes[row + 1, 2] = last(pair)
    end
end

println("Wrote ", SVG_OUT)
println("Wrote ", AUDIT_OUT)
