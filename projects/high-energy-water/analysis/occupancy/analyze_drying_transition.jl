#!/usr/bin/env julia

using XLSX
using DelimitedFiles
using Random
using Statistics

const HERE = @__DIR__
const PROJECT_ROOT = normpath(joinpath(HERE, "..", ".."))
const INPUT = joinpath(PROJECT_ROOT, "data", "processed", "occupancy_along_binding_coordinate.xlsx")
const OUTPUT = joinpath(PROJECT_ROOT, "data", "derived", "drying_transition_analysis.xlsx")

parse_z(value) = parse(Float64, replace(string(value), "z=" => ""))

function moments(n, probability)
    p = probability ./ sum(probability)
    sum(n .* p)
end

function read_probability_curve(sheet; header_row, first_col, last_col, first_row, n_values)
    z = Float64[parse_z(sheet[header_row, col]) for col in first_col:last_col]
    mean_n = Float64[]
    for col in first_col:last_col
        probability = Float64[
            ismissing(sheet[first_row + row - 1, col]) ? 0.0 : sheet[first_row + row - 1, col]
            for row in eachindex(n_values)
        ]
        push!(mean_n, moments(n_values, probability))
    end
    z, mean_n
end

function read_sem(path)
    raw = readdlm(path, '\t', Float64; skipstart=1)
    Dict(raw[row, 1] => raw[row, 2] for row in axes(raw, 1))
end

function crossing(distance, occupancy, target; nominal=nothing)
    candidates = Float64[]
    for i in 1:length(distance)-1
        y1, y2 = occupancy[i], occupancy[i + 1]
        if (target - y1) * (target - y2) <= 0 && y1 != y2
            fraction = (target - y1) / (y2 - y1)
            if 0 <= fraction <= 1
                push!(candidates, distance[i] + fraction * (distance[i + 1] - distance[i]))
            end
        end
    end
    isempty(candidates) && return NaN
    nominal === nothing ? first(candidates) : candidates[argmin(abs.(candidates .- nominal))]
end

struct Branch
    host::String
    portal::String
    z::Vector{Float64}
    distance::Vector{Float64}
    mean_n::Vector{Float64}
    sem_n::Vector{Float64}
    wet_mask::BitVector
    dry_mask::BitVector
    sign::Float64
end

function summarize(branch::Branch)
    wet = mean(branch.mean_n[branch.wet_mask])
    dry = mean(branch.mean_n[branch.dry_mask])
    span = wet - dry
    target25 = wet - 0.25span
    target50 = wet - 0.50span
    target75 = wet - 0.75span
    order = sortperm(branch.distance)
    d = branch.distance[order]
    n = branch.mean_n[order]
    d25 = crossing(d, n, target25)
    d50 = crossing(d, n, target50)
    d75 = crossing(d, n, target75)
    width = d25 - d75
    (
        wet=wet,
        dry=dry,
        target25=target25,
        target50=target50,
        target75=target75,
        d25=d25,
        d50=d50,
        d75=d75,
        z25=branch.sign * d25,
        z50=branch.sign * d50,
        z75=branch.sign * d75,
        width=width,
        normalized_slope=0.5 / width,
        occupancy_slope=0.5span / width,
    )
end

function propagated_widths(branch::Branch, nominal; trials=50_000, seed=20260917)
    rng = MersenneTwister(seed)
    order = sortperm(branch.distance)
    d = branch.distance[order]
    widths = Float64[]
    midpoints = Float64[]
    for _ in 1:trials
        sampled = max.(0.0, branch.mean_n .+ branch.sem_n .* randn(rng, length(branch.mean_n)))
        wet = mean(sampled[branch.wet_mask])
        dry = mean(sampled[branch.dry_mask])
        wet > dry || continue
        span = wet - dry
        sampled_ordered = sampled[order]
        d25 = crossing(d, sampled_ordered, wet - 0.25span; nominal=nominal.d25)
        d50 = crossing(d, sampled_ordered, wet - 0.50span; nominal=nominal.d50)
        d75 = crossing(d, sampled_ordered, wet - 0.75span; nominal=nominal.d75)
        if all(isfinite, (d25, d50, d75)) && d25 > d75
            push!(widths, d25 - d75)
            push!(midpoints, branch.sign * d50)
        end
    end
    length(widths) > 0.9trials || error("Too few valid propagated transitions for $(branch.host) $(branch.portal)")
    sort!(widths)
    sort!(midpoints)
    q(v, p) = quantile(v, p)
    (
        width_se=std(widths),
        width_low=q(widths, 0.025),
        width_high=q(widths, 0.975),
        midpoint_se=std(midpoints),
        midpoint_low=q(midpoints, 0.025),
        midpoint_high=q(midpoints, 0.975),
        accepted=length(widths),
        trials=trials,
    )
end

cb7_z = Float64[]; cb7_mean = Float64[]
oa_z = Float64[]; oa_mean = Float64[]; oa_sem = Float64[]
cd_z = Float64[]; cd_mean = Float64[]

XLSX.openxlsx(INPUT) do workbook
    global cb7_z, cb7_mean = read_probability_curve(
        workbook["CB7"];
        header_row=2, first_col=8, last_col=33, first_row=3,
        n_values=collect(0.0:9.0),
    )
    oa = workbook["OA"]
    global oa_z = Float64[oa[row, 1] for row in 4:37]
    global oa_mean = Float64[oa[row, 2] for row in 4:37]
    global oa_sem = Float64[oa[row, 4] for row in 4:37]
    order = sortperm(oa_z)
    global oa_z, oa_mean, oa_sem = oa_z[order], oa_mean[order], oa_sem[order]

    cd = workbook["CD_Restrained"]
    n_values = Float64[cd[row, 8] for row in 4:13]
    global cd_z, cd_mean = read_probability_curve(
        cd;
        header_row=3, first_col=9, last_col=59, first_row=4,
        n_values=n_values,
    )
end

cb7_sem_lookup = read_sem(joinpath(PROJECT_ROOT, "data", "derived", "cb7_cluster_bootstrap_sem.tsv"))
cd_sem_lookup = read_sem(joinpath(PROJECT_ROOT, "data", "derived", "betacd_cluster_bootstrap_sem.tsv"))
cb7_sem = Float64[cb7_sem_lookup[z] for z in cb7_z]
cd_sem = Float64[cd_sem_lookup[z] for z in cd_z]

branches = Branch[
    Branch(
        "CB7", "symmetric portal", cb7_z, cb7_z, cb7_mean, cb7_sem,
        cb7_z .>= 6.0, cb7_z .<= 0.8, 1.0,
    ),
    Branch(
        "OA", "open portal", oa_z[oa_z .>= 0], oa_z[oa_z .>= 0],
        oa_mean[oa_z .>= 0], oa_sem[oa_z .>= 0],
        oa_z[oa_z .>= 0] .>= 6.0, oa_z[oa_z .>= 0] .<= 0.8, 1.0,
    ),
    Branch(
        "restrained β-CD", "positive-z wider secondary portal",
        cd_z[cd_z .>= 0], cd_z[cd_z .>= 0], cd_mean[cd_z .>= 0], cd_sem[cd_z .>= 0],
        cd_z[cd_z .>= 0] .>= 6.0, cd_z[cd_z .>= 0] .<= 0.8, 1.0,
    ),
    Branch(
        "restrained β-CD", "negative-z narrower primary portal",
        cd_z[cd_z .<= 0], -cd_z[cd_z .<= 0], cd_mean[cd_z .<= 0], cd_sem[cd_z .<= 0],
        cd_z[cd_z .<= 0] .<= -6.0, cd_z[cd_z .<= 0] .>= -0.8, -1.0,
    ),
]

summaries = [summarize(branch) for branch in branches]
uncertainties = [
    propagated_widths(branch, summary; seed=20260917 + i)
    for (i, (branch, summary)) in enumerate(zip(branches, summaries))
]
oa_width = summaries[2].width

XLSX.openxlsx(OUTPUT, mode="w") do workbook
    summary_sheet = workbook[1]
    XLSX.rename!(summary_sheet, "Summary")
    headers = [
        "host", "branch", "N_wet", "N_dry", "z_25pct_empty_A", "z_50pct_empty_A",
        "z_75pct_empty_A", "transition_width_A", "width_relative_to_OA",
        "normalized_slope_abs_per_A", "occupancy_slope_abs_waters_per_A",
        "approx_width_SE_A", "approx_width_95pct_low_A", "approx_width_95pct_high_A",
        "approx_midpoint_SE_A", "propagation_trials_accepted", "propagation_trials_total",
    ]
    summary_sheet["A1:Q1"] = reshape(headers, 1, :)
    for (index, (branch, summary, uncertainty)) in enumerate(zip(branches, summaries, uncertainties))
        row = index + 1
        values = Any[
            branch.host, branch.portal, summary.wet, summary.dry, summary.z25, summary.z50,
            summary.z75, summary.width, summary.width / oa_width,
            summary.normalized_slope, summary.occupancy_slope,
            uncertainty.width_se, uncertainty.width_low, uncertainty.width_high,
            uncertainty.midpoint_se, uncertainty.accepted, uncertainty.trials,
        ]
        for (column, value) in enumerate(values)
            summary_sheet[row, column] = value
        end
    end

    curve_sheet = XLSX.addsheet!(workbook, "Normalized curves")
    curve_sheet["A1:G1"] = reshape(
        ["host", "branch", "z_A", "mean_N", "SEM_N", "fraction_empty", "plateau_definition"],
        1, :,
    )
    row = 2
    for (branch, summary) in zip(branches, summaries)
        fraction_empty = (summary.wet .- branch.mean_n) ./ (summary.wet - summary.dry)
        for i in eachindex(branch.z)
            values = Any[
                branch.host, branch.portal, branch.z[i], branch.mean_n[i], branch.sem_n[i],
                fraction_empty[i], "wet: |z| >= 6 A; dry: |z| <= 0.8 A",
            ]
            for (column, value) in enumerate(values)
                curve_sheet[row, column] = value
            end
            row += 1
        end
    end

    notes = XLSX.addsheet!(workbook, "Definitions")
    notes["A1:B1"] = reshape(["quantity", "definition"], 1, :)
    definitions = [
        "fraction empty" => "(N_wet - <N>(z)) / (N_wet - N_dry)",
        "25-75% transition width" => "absolute separation in z between fraction-empty = 0.25 and 0.75, obtained by linear interpolation between adjacent 0.4-A bins",
        "normalized slope" => "0.50 / transition width; larger values indicate a sharper spatial transition",
        "occupancy slope" => "0.50*(N_wet - N_dry) / transition width, in waters per A",
        "wet plateau" => "arithmetic mean over |z| >= 6 A on the relevant branch",
        "dry plateau" => "arithmetic mean over |z| <= 0.8 A on the relevant branch",
        "uncertainty" => "approximate pointwise-SEM propagation with 50,000 independent Gaussian draws; covariance between z bins was unavailable",
        "interpretation" => "sharp/gradual refers to change with guest position, not kinetic drying rate",
    ]
    for (index, pair) in enumerate(definitions)
        row = index + 1
        notes[row, 1] = first(pair)
        notes[row, 2] = last(pair)
    end
end

for (branch, summary, uncertainty) in zip(branches, summaries, uncertainties)
    println(
        branch.host, " | ", branch.portal,
        ": width = ", round(summary.width, digits=3), " A (approx. 95% ",
        round(uncertainty.width_low, digits=3), "–", round(uncertainty.width_high, digits=3),
        "), |df/dz| = ", round(summary.normalized_slope, digits=3), " A^-1, midpoint z = ",
        round(summary.z50, digits=3), " A",
    )
end
println("Wrote ", OUTPUT)
