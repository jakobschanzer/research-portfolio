#!/usr/bin/env julia

ENV["GKSwstype"] = "100"
using XLSX
using Plots
using LaTeXStrings

const HERE = @__DIR__
const PROJECT_ROOT = normpath(joinpath(HERE, "..", ".."))
const THERMO = joinpath(PROJECT_ROOT, "data", "processed", "pmf_thermodynamics.xlsx")
const FORCES = joinpath(PROJECT_ROOT, "data", "processed", "direct_residual_decomposition.xlsx")
const OUTPUT_DIR = joinpath(PROJECT_ROOT, "figures")

const HOSTS = [
    (name="CB7", sheet="CB7", force_cols=(1,2,3,4), zfactor=10.0, xlim=(-12,12), ylim=(-120,60), yticks=-120:20:60),
    (name="OA", sheet="OA", force_cols=(7,8,9,10), zfactor=1.0, xlim=(-4,12), ylim=(-80,120), yticks=-80:20:120),
    (name="β-CD", sheet="BetaCD", force_cols=(15,16,17,18), zfactor=10.0, xlim=(-12,12), ylim=(-80,40), yticks=-80:20:40),
]

const BLACK = RGB(0.08,0.08,0.08)
const VERMILION = RGB(0.80,0.24,0.12)
const BLUE = RGB(0.05,0.35,0.70)
const GRAY = RGB(0.42,0.42,0.42)

function numeric_columns(sh, columns; firstrow=2)
    output = [Float64[] for _ in columns]
    row = firstrow
    while true
        vals = [sh[row,c] for c in columns]
        if any(ismissing, vals)
            break
        end
        try
            for (dest,val) in zip(output,vals)
                push!(dest, Float64(val))
            end
        catch
            break
        end
        row += 1
    end
    output
end

thermo = Dict{String,NamedTuple}()
XLSX.openxlsx(THERMO) do xf
    for h in HOSTS
        z,g,enthalpy,minus_ts = numeric_columns(xf[h.sheet], (1,2,3,5); firstrow=2)
        thermo[h.name] = (z=z, g=g, h=enthalpy, mts=minus_ts)
    end
end

forces = Dict{String,NamedTuple}()
XLSX.openxlsx(FORCES) do xf
    sh=xf["Solvent_Effect"]
    for h in HOSTS
        z,g,direct,indirect = numeric_columns(sh, h.force_cols; firstrow=3)
        forces[h.name] = (z=h.zfactor .* z, g=g, direct=direct, indirect=indirect)
    end
end

default(
    fontfamily="DejaVu Sans",
    guidefontsize=15,
    tickfontsize=11,
    legendfontsize=10,
    linewidth=1.8,
    foreground_color_subplot=:black,
    background_color=:white,
    background_color_legend=:transparent,
    foreground_color_legend=:transparent,
)

panels = Any[]
for (i,h) in enumerate(HOSTS)
    p=plot(
        xlims=h.xlim, ylims=h.ylim, yticks=h.yticks,
        xticks=h.xlim[1]:4:h.xlim[2],
        xlabel=L"z\ (\mathrm{\AA})", ylabel=i==1 ? L"\mathrm{Energy}\ (\mathrm{kJ\ mol}^{-1})" : "",
        title=h.name, titlefontsize=14, titlelocation=:right,
        framestyle=:box, grid=false,
        legend=i==1 ? :bottomleft : false,
        margin=4Plots.mm,
    )
    d=thermo[h.name]
    plot!(p,d.z,d.g; color=BLACK, linestyle=:solid, label=L"\Delta G(z)")
    plot!(p,d.z,d.h; color=VERMILION, linestyle=:dash, label=L"\Delta H(z)")
    plot!(p,d.z,d.mts; color=BLUE, linestyle=:dot, label=L"-T\Delta S(z)")
    hline!(p,[0.0]; color=RGB(.72,.72,.72), linewidth=.8, label="")
    push!(panels,p)
end

for (i,h) in enumerate(HOSTS)
    p=plot(
        xlims=h.xlim, ylims=h.ylim, yticks=h.yticks,
        xticks=h.xlim[1]:4:h.xlim[2],
        xlabel=L"z\ (\mathrm{\AA})", ylabel=i==1 ? L"\mathrm{Energy}\ (\mathrm{kJ\ mol}^{-1})" : "",
        framestyle=:box, grid=false,
        legend=i==1 ? :bottomleft : false,
        margin=4Plots.mm,
    )
    d=forces[h.name]
    plot!(p,d.z,d.g; color=BLACK, linestyle=:solid, label=L"\Delta G(z)")
    plot!(p,d.z,d.direct; color=VERMILION, linestyle=:dash, label=L"\mathrm{Direct}")
    plot!(p,d.z,d.indirect; color=BLUE, linestyle=:dot, label=L"\mathrm{Non\!\!\!-direct\ residual}")
    hline!(p,[0.0]; color=RGB(.72,.72,.72), linewidth=.8, label="")
    push!(panels,p)
end

annotate!(panels[1], HOSTS[1].xlim[1]+0.6, HOSTS[1].ylim[2]-9, text("a",17,:bold,:left))
annotate!(panels[4], HOSTS[1].xlim[1]+0.6, HOSTS[1].ylim[2]-9, text("b",17,:bold,:left))

fig=plot(panels...;
    layout=(2,3), size=(1500,930), dpi=600,
    left_margin=6Plots.mm, right_margin=4Plots.mm,
    bottom_margin=6Plots.mm, top_margin=4Plots.mm,
)

for ext in ("png","pdf","svg")
    savefig(fig,joinpath(OUTPUT_DIR,"figure2_pmf_decomposition.$ext"))
end
println("Wrote figures/figure2_pmf_decomposition.png, .pdf, and .svg")
