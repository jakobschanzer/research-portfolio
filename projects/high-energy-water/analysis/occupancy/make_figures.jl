#!/usr/bin/env julia

ENV["GKSwstype"]="100"
using XLSX
using Plots
using DelimitedFiles
using LaTeXStrings

const HERE=@__DIR__
const PROJECT_ROOT=normpath(joinpath(HERE,"..",".."))
const INPUT=joinpath(PROJECT_ROOT,"data","processed","occupancy_along_binding_coordinate.xlsx")
const MEAN_INPUT=joinpath(PROJECT_ROOT,"data","processed","mean_occupancy_vs_z.xlsx")
const OUTPUT_DIR=joinpath(PROJECT_ROOT,"figures")
const DERIVED_DIR=joinpath(PROJECT_ROOT,"data","derived")

struct OccupancyData
    z::Vector{Float64}
    n::Vector{Float64}
    p::Matrix{Float64} # N × z
    mean::Vector{Float64}
    variance::Vector{Float64}
    sd::Vector{Float64}
end

struct CurveData
    z::Vector{Float64}
    mean::Vector{Float64}
    sd::Vector{Float64}
    sem::Vector{Float64}
end

function parse_z(header)
    parse(Float64,replace(string(header),"z="=>""))
end

function moments(n,p)
    normalized=p ./ sum(p)
    μ=sum(n .* normalized)
    variance=sum((n .- μ).^2 .* normalized)
    μ,variance
end

function read_host(sh; header_row, first_col, last_col, n_rows, explicit_n_col=nothing)
    z=Float64[parse_z(sh[header_row,c]) for c in first_col:last_col]
    n=explicit_n_col===nothing ? collect(0.0:n_rows-1) : Float64[sh[header_row+r,explicit_n_col] for r in 1:n_rows]
    p=zeros(n_rows,length(z))
    for (j,c) in enumerate(first_col:last_col), r in 1:n_rows
        value=sh[header_row+r,c]
        p[r,j]=ismissing(value) ? 0.0 : Float64(value)
    end
    mean=Float64[]; variance=Float64[]
    for j in eachindex(z)
        μ,v=moments(n,p[:,j]); push!(mean,μ); push!(variance,v)
    end
    OccupancyData(z,n,p,mean,variance,sqrt.(variance))
end

function read_curve(sh,first_row,last_row;sem_col,variance_col=nothing,sd_col=nothing)
    curve_z=Float64[sh[r,1] for r in first_row:last_row]
    curve_mean=Float64[sh[r,2] for r in first_row:last_row]
    curve_sd=sd_col===nothing ? fill(NaN,length(curve_z)) : Float64[sh[r,sd_col] for r in first_row:last_row]
    curve_sem=Float64[sh[r,sem_col] for r in first_row:last_row]
    order=sortperm(curve_z)
    CurveData(curve_z[order],curve_mean[order],curve_sd[order],curve_sem[order])
end

"""Read OA P(N|z), correcting its shifted headers by matching first moments."""
function read_oa(sh,curve)
    n=Float64[sh[r,2] for r in 3:12]
    populated=Int[]
    raw_headers=String[]
    probability_columns=Vector{Vector{Float64}}()
    for c in 3:34
        values=[sh[r,c] for r in 3:12]
        if all(!ismissing,values) && sum(Float64.(values))>0
            push!(populated,c)
            push!(raw_headers,string(sh[2,c]))
            push!(probability_columns,Float64.(values))
        end
    end
    p=hcat(probability_columns...)
    derived_mean=Float64[moments(n,p[:,j])[1] for j in axes(p,2)]

    # The OA P(N|z) headers are one 0.4-A window behind the data. Determine the
    # actual contiguous coordinate span from the independently tabulated <N>(z).
    candidate_rmse=Float64[]
    for first in 1:(length(curve.z)-length(derived_mean)+1)
        push!(candidate_rmse,sqrt(sum((derived_mean .- curve.mean[first:first+length(derived_mean)-1]).^2)/length(derived_mean)))
    end
    first=argmin(candidate_rmse)
    actual_z=curve.z[first:first+length(derived_mean)-1]
    matched_mean=curve.mean[first:first+length(derived_mean)-1]
    maximum(abs.(derived_mean .- matched_mean)) < 5e-4 || error("OA P(N|z) columns cannot be aligned to the tabulated means")

    variance=Float64[moments(n,p[:,j])[2] for j in axes(p,2)]
    data=OccupancyData(actual_z,n,p,derived_mean,variance,sqrt.(variance))
    data,raw_headers,matched_mean
end

cb7=nothing; cb7_curve=nothing; oa=nothing; oa_curve=nothing; oa_raw_headers=nothing; oa_matched_mean=nothing; cd=nothing; cd_curve=nothing
XLSX.openxlsx(MEAN_INPUT) do xf
    global cb7_curve=read_curve(xf["CB7"],2,52;sem_col=3)
    global oa_curve=read_curve(xf["OA"],2,35;sem_col=5,variance_col=3,sd_col=4)
    global cd_curve=read_curve(xf["restrained beta-CD"],2,52;sem_col=5,variance_col=3,sd_col=4)
end
XLSX.openxlsx(INPUT) do xf
    global cb7=read_host(xf["CB7"];header_row=2,first_col=2,last_col=27,n_rows=10)
    global oa,oa_raw_headers,oa_matched_mean=read_oa(xf["OA"],oa_curve)
    global cd=read_host(xf["Beta"];header_row=2,first_col=3,last_col=53,n_rows=11,explicit_n_col=2)
end

positive=findall(>=(0.0),cb7_curve.z)
maximum(abs.(cb7_curve.mean[positive] .- cb7.mean)) < 5e-4 ||
    error("CB7 P(N|z) moments do not match the corrected mean-occupancy curve")
maximum(abs.(cd_curve.mean .- cd.mean)) < 5e-4 ||
    error("beta-CD P(N|z) moments do not match the corrected mean-occupancy curve")
cb7_curve=CurveData(
    cb7_curve.z,
    cb7_curve.mean,
    Float64[cb7.sd[argmin(abs.(cb7.z .- abs(value)))] for value in cb7_curve.z],
    cb7_curve.sem,
)

function read_transitions(path)
    transitions=Dict{Tuple{String,String},NamedTuple}()
    XLSX.openxlsx(path) do xf
        sh=xf["Summary"]
        for row in 2:5
            host=string(sh[row,1]); branch=string(sh[row,2])
            wet=Float64(sh[row,3]); dry=Float64(sh[row,4]); span=wet-dry
            z25=Float64(sh[row,5]); z50=Float64(sh[row,6]); z75=Float64(sh[row,7])
            transitions[(host,branch)]=(
                lower=min(z25,z75),upper=max(z25,z75),z25=z25,z75=z75,
                z50=z50,n25=wet-0.25span,n75=wet-0.75span,
                slope=Float64(sh[row,11]),
            )
        end
    end
    transitions
end

transitions=read_transitions(joinpath(DERIVED_DIR,"drying_transition_analysis.xlsx"))
cb7_positive=transitions[("CB7","symmetric portal")]
cb7_negative=(lower=-cb7_positive.upper,upper=-cb7_positive.lower,
    z25=-cb7_positive.z25,z50=-cb7_positive.z50,z75=-cb7_positive.z75,n25=cb7_positive.n25,
    n75=cb7_positive.n75,slope=cb7_positive.slope)
cb7_transitions=[cb7_negative,cb7_positive]
oa_transitions=[transitions[("OA","open portal")]]
cd_transitions=[
    transitions[("restrained β-CD","negative-z narrower primary portal")],
    transitions[("restrained β-CD","positive-z wider secondary portal")],
]

const CB7_SELECTED=[2.8,3.2,3.6,4.0,4.4,4.8]
const OA_SELECTED=[9.6,4.0,3.6,3.2,2.8,0.0]
const OA_DISPLAY_Z=Dict(9.6=>10.0)
const CD_SELECTED=[-10.0,-5.2,-4.8,-3.6,-2.0,0.0,2.0,2.8,3.6,4.4,10.0]

default(fontfamily="DejaVu Sans",foreground_color_subplot=:black,
    background_color=:white,grid=false,linewidth=1.5,
    guidefontsize=14,tickfontsize=10,titlefontsize=11)

function selected_index(data,z)
    i=argmin(abs.(data.z .- z))
    abs(data.z[i]-z)<1e-8 || error("No z=$z")
    i
end

function nz_panel(z,mean,sem,title_text,letter,xlims,xticks,transition_data)
    p=plot(;xlabel=L"z\ (\mathrm{\AA})",ylabel=L"\langle N\rangle",
        xlims=xlims,xticks=xticks,ylims=(0,7),yticks=0:1:7,
        title=title_text,framestyle=:box,margin=6Plots.mm)
    for (index,transition) in enumerate(transition_data)
        vspan!(p,[transition.lower,transition.upper];color=RGB(.72,.72,.72),alpha=.32,label="")
        vline!(p,[transition.lower,transition.upper];color=RGB(.42,.42,.42),linestyle=:dot,
            linewidth=.8,label="")
        vline!(p,[transition.z50];color=RGB(.25,.25,.25),linestyle=:dash,
            linewidth=1.0,label="")
    end
    plot!(p,z,mean;yerror=sem,color=:black,marker=:circle,
        markercolor=:white,markerstrokecolor=:black,markersize=4,
        errorbar_color=RGB(.30,.30,.30),errorbar_linewidth=1.0,
        label=L"\langle N\rangle\pm\mathrm{SEM}",legend=:bottomright,
    )
    annotate!(p,xlims[1]+0.035*(xlims[2]-xlims[1]),6.68,text(letter,17,"DejaVu Sans",:left))
    p
end

p1=nz_panel(cb7_curve.z,cb7_curve.mean,cb7_curve.sem,"CB7","a",(-10,10),-10:5:10,cb7_transitions)
p2=nz_panel(oa_curve.z,oa_curve.mean,oa_curve.sem,"OA","b",(-3,11),-2:2:10,oa_transitions)
p3=nz_panel(cd_curve.z,cd_curve.mean,cd_curve.sem,"restrained β-CD","c",(-10,10),-10:5:10,cd_transitions)
nz=plot(p1,p2,p3;layout=(1,3),size=(1800,560),dpi=600,
    left_margin=6Plots.mm,right_margin=5Plots.mm,bottom_margin=7Plots.mm,top_margin=5Plots.mm)
for ext in ("png","pdf","svg"); savefig(nz,joinpath(OUTPUT_DIR,"figure5_mean_occupancy_vs_z.$ext")); end

letters=[string(Char('a'+i)) for i in 0:25]
function distribution_figure(data,selected,host;logscale=false,layout=(3,4),stem="",display_z=Dict{Float64,Float64}())
    panels=Any[]
    for (k,zvalue) in enumerate(selected)
        j=selected_index(data,zvalue); probability=data.p[:,j] ./ sum(data.p[:,j])
        if logscale
            keep=probability .> 0
            x=data.n[keep]; y=log.(probability[keep])
            p=plot(x,y;color=:black,marker=:circle,markercolor=:black,markersize=3,
                label="",ylims=(-12,0.3),yticks=-12:2:0)
        else
            p=plot(data.n,probability;color=:black,marker=:circle,
                markercolor=:white,markerstrokecolor=:black,markersize=3,
                label="",ylims=(0,1.02),yticks=0:0.2:1.0)
        end
        zlabel=get(display_z,zvalue,zvalue)
        plot!(p;xlims=(-0.4,9.4),xticks=0:2:8,framestyle=:box,
            title=latexstring("\\mathbf{$(letters[k])}\\quad z = $(round(zlabel,digits=1))\\ \\mathrm{\\AA}"),
            xlabel=k>length(selected)-layout[2] ? L"N" : "",
            ylabel=(k-1)%layout[2]==0 ? (logscale ? L"\ln P(N\mid z)" : L"P(N\mid z)") : "",
            margin=3Plots.mm)
        push!(panels,p)
    end
    fig=plot(panels...;layout=layout,
        size=(layout[2]*310,layout[1]*280),dpi=600,
        plot_title=host,plot_titlefontsize=16,
        left_margin=4Plots.mm,right_margin=3Plots.mm,
        bottom_margin=4Plots.mm,top_margin=3Plots.mm)
    for ext in ("png","pdf","svg"); savefig(fig,joinpath(OUTPUT_DIR,"$stem.$ext")); end
end

distribution_figure(cb7,CB7_SELECTED,"CB7";layout=(3,4),stem="cb7_selected_pn")
distribution_figure(cb7,CB7_SELECTED,"CB7";logscale=true,layout=(3,4),stem="cb7_selected_lnpn")
distribution_figure(oa,OA_SELECTED,"OA";layout=(3,4),stem="oa_selected_pn",display_z=OA_DISPLAY_Z)
distribution_figure(oa,OA_SELECTED,"OA";logscale=true,layout=(3,4),stem="oa_selected_lnpn",display_z=OA_DISPLAY_Z)
distribution_figure(cd,CD_SELECTED,"restrained β-CD";layout=(3,5),stem="betacd_selected_pn")
distribution_figure(cd,CD_SELECTED,"restrained β-CD";logscale=true,layout=(3,5),stem="betacd_selected_lnpn")

function paired_distribution_figure(data,selected,host;stem="",display_z=Dict{Float64,Float64}())
    panels=Any[]
    for (k,zvalue) in enumerate(selected)
        j=selected_index(data,zvalue)
        probability=data.p[:,j] ./ sum(data.p[:,j])
        zlabel=get(display_z,zvalue,zvalue)
        p_linear=plot(data.n,probability;color=:black,marker=:circle,
            markercolor=:white,markerstrokecolor=:black,markersize=3,label="",
            xlims=(-0.4,9.4),xticks=0:2:8,ylims=(0,1.02),yticks=0:0.2:1.0,
            title=latexstring("\\mathbf{$(letters[k])}\\quad z = $(round(zlabel,digits=1))\\ \\mathrm{\\AA}"),
            xlabel=L"N",ylabel=L"P(N\mid z)",framestyle=:box,margin=3Plots.mm)
        keep=probability .> 0
        p_log=plot(data.n[keep],log.(probability[keep]);color=:black,marker=:circle,
            markercolor=:black,markersize=3,label="",xlims=(-0.4,9.4),xticks=0:2:8,
            ylims=(-12,0.3),yticks=-12:2:0,title=" ",
            xlabel=L"N",ylabel=L"\ln P(N\mid z)",framestyle=:box,margin=3Plots.mm)
        push!(panels,p_linear,p_log)
    end
    rows=cld(length(selected),2)
    fig=plot(panels...;layout=(rows,4),size=(1560,rows*300),dpi=600,
        plot_title=host,plot_titlefontsize=16,left_margin=4Plots.mm,
        right_margin=3Plots.mm,bottom_margin=4Plots.mm,top_margin=3Plots.mm)
    for ext in ("png","pdf","svg"); savefig(fig,joinpath(OUTPUT_DIR,"$stem.$ext")); end
end

paired_distribution_figure(cb7,CB7_SELECTED,"CB7";stem="cb7_pn_lnpn_paired")
paired_distribution_figure(oa,OA_SELECTED,"OA";stem="oa_pn_lnpn_paired",display_z=OA_DISPLAY_Z)
paired_distribution_figure(cd,CD_SELECTED,"restrained β-CD";stem="betacd_pn_lnpn_paired")

# Machine-readable audit table.
XLSX.openxlsx(joinpath(DERIVED_DIR,"occupancy_moments.xlsx"),mode="w") do xf
    records=(("CB7",cb7_curve.z,cb7_curve.mean,cb7_curve.sd,cb7_curve.sem,"rerun direct tabulation"),
             ("OA",oa_curve.z,oa_curve.mean,oa_curve.sd,oa_curve.sem,"source std err column"),
             ("restrained beta-CD",cd_curve.z,cd_curve.mean,cd_curve.sd,cd_curve.sem,
                "corrected direct tabulation from matched z/N frames"))
    for (idx,(name,z,mean,sd,sem,sem_source)) in enumerate(records)
        sh=idx==1 ? xf[1] : XLSX.addsheet!(xf,name)
        idx==1 && XLSX.rename!(sh,name)
        sh["A1:F1"]=reshape(["z_A","mean_N","variance_N","SD_N","SEM_N","SEM_source"],1,:)
        for i in eachindex(z)
            sh[i+1,1]=z[i]; sh[i+1,2]=mean[i]
            sh[i+1,3]=sd[i]^2; sh[i+1,4]=sd[i]
            sh[i+1,5]=sem[i]; sh[i+1,6]=sem_source
        end
    end

    audit=XLSX.addsheet!(xf,"OA Pn alignment")
    audit["A1:F1"]=reshape(["raw_workbook_header","source_z_A","figure_z_A","mean_from_Pn","tabulated_mean_N","absolute_difference"],1,:)
    for i in eachindex(oa.z)
        audit[i+1,1]=oa_raw_headers[i]
        audit[i+1,2]=oa.z[i]
        audit[i+1,3]=get(OA_DISPLAY_Z,oa.z[i],oa.z[i])
        audit[i+1,4]=oa.mean[i]
        audit[i+1,5]=oa_matched_mean[i]
        audit[i+1,6]=abs(oa.mean[i]-oa_matched_mean[i])
    end
end
println("OA P(N|z) aligned to z=$(first(oa.z)) to $(last(oa.z)) A; maximum mean mismatch=$(maximum(abs.(oa.mean.-oa_matched_mean)))")
println("Wrote N(z), selected P(N|z), selected ln P(N|z), and occupancy_moments.xlsx")
