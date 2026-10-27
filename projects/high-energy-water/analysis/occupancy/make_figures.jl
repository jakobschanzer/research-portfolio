#!/usr/bin/env julia

ENV["GKSwstype"]="100"
using XLSX
using Plots
using DelimitedFiles
using LaTeXStrings

const HERE=@__DIR__
const PROJECT_ROOT=normpath(joinpath(HERE,"..",".."))
const INPUT=joinpath(PROJECT_ROOT,"data","processed","occupancy_along_binding_coordinate.xlsx")
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

"""Read OA while correcting the shifted P(N|z) header row by moment matching."""
function read_oa(sh)
    curve_z=Float64[sh[r,1] for r in 4:37]
    curve_mean=Float64[sh[r,2] for r in 4:37]
    curve_sd=Float64[sh[r,3] for r in 4:37]
    # Column D is explicitly labelled "std err". Column C is the instantaneous SD.
    curve_sem=Float64[sh[r,4] for r in 4:37]
    order=sortperm(curve_z)
    curve=CurveData(curve_z[order],curve_mean[order],curve_sd[order],curve_sem[order])

    n=Float64[sh[r,7] for r in 4:13]
    populated=Int[]
    raw_headers=String[]
    probability_columns=Vector{Vector{Float64}}()
    for c in 8:39
        values=[sh[r,c] for r in 4:13]
        if all(!ismissing,values) && sum(Float64.(values))>0
            push!(populated,c)
            push!(raw_headers,string(sh[3,c]))
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
    data,curve,raw_headers,matched_mean
end

cb7=nothing; oa=nothing; oa_curve=nothing; oa_raw_headers=nothing; oa_matched_mean=nothing; cd=nothing
XLSX.openxlsx(INPUT) do xf
    global cb7=read_host(xf["CB7"];header_row=2,first_col=8,last_col=33,n_rows=10)
    global oa,oa_curve,oa_raw_headers,oa_matched_mean=read_oa(xf["OA"])
    global cd=read_host(xf["CD_Restrained"];header_row=3,first_col=9,last_col=59,n_rows=10,explicit_n_col=8)
end

# CB7 was sampled on the non-negative half-axis and is symmetric about z = 0.
# Mirror every positive-z P(N|z) column onto negative z, retaining z = 0 once.
function mirror_about_zero(data)
    positive=findall(>(0.0),data.z)
    order=reverse(positive)
    z=vcat(-data.z[order],data.z)
    p=hcat(data.p[:,order],data.p)
    OccupancyData(z,data.n,p,vcat(data.mean[order],data.mean),
        vcat(data.variance[order],data.variance),vcat(data.sd[order],data.sd))
end
cb7=mirror_about_zero(cb7)

function read_sem(path)
    raw=readdlm(path,'\t',Float64;skipstart=1)
    Dict(raw[i,1]=>raw[i,2] for i in axes(raw,1))
end

cb7_sem_positive=read_sem(joinpath(DERIVED_DIR,"cb7_cluster_bootstrap_sem.tsv"))
cd_sem_by_z=read_sem(joinpath(DERIVED_DIR,"betacd_cluster_bootstrap_sem.tsv"))
cb7_sem_by_z=Dict(z=>cb7_sem_positive[abs(z)] for z in cb7.z)

function aligned_sem(data,lookup)
    [get(lookup,z,NaN) for z in data.z]
end
cb7_sem=aligned_sem(cb7,cb7_sem_by_z)
cd_sem=aligned_sem(cd,cd_sem_by_z)
any(isnan,cb7_sem) && error("Missing CB7 SEM value")
any(isnan,cd_sem) && error("Missing beta-CD SEM value")

function read_transition_bands(path)
    bands=Dict{Tuple{String,String},Tuple{Float64,Float64}}()
    XLSX.openxlsx(path) do xf
        sh=xf["Summary"]
        for row in 2:5
            host=string(sh[row,1]); branch=string(sh[row,2])
            z25=Float64(sh[row,5]); z75=Float64(sh[row,7])
            bands[(host,branch)]=(min(z25,z75),max(z25,z75))
        end
    end
    bands
end

transition_bands=read_transition_bands(joinpath(DERIVED_DIR,"drying_transition_analysis.xlsx"))
cb7_positive=transition_bands[("CB7","symmetric portal")]
cb7_bands=[(-cb7_positive[2],-cb7_positive[1]),cb7_positive]
oa_bands=[transition_bands[("OA","open portal")]]
cd_bands=[
    transition_bands[("restrained β-CD","negative-z narrower primary portal")],
    transition_bands[("restrained β-CD","positive-z wider secondary portal")],
]

const CB7_SELECTED=[2.8,3.2,3.6,4.0,4.4,4.8]
const OA_SELECTED=[9.6,4.0,3.6,3.2,2.8,0.0]
const OA_DISPLAY_Z=Dict(9.6=>10.0)
const CD_SELECTED=[-4.4,-3.6,-2.8,-2.0,2.4,2.8,3.6,4.4]

default(fontfamily="DejaVu Sans",foreground_color_subplot=:black,
    background_color=:white,grid=false,linewidth=1.5,
    guidefontsize=14,tickfontsize=10,titlefontsize=11)

function selected_index(data,z)
    i=argmin(abs.(data.z .- z))
    abs(data.z[i]-z)<1e-8 || error("No z=$z")
    i
end

function nz_panel(z,mean,sem,title_text,letter,xlims,xticks,bands)
    p=plot(;xlabel=L"z\ (\mathrm{\AA})",ylabel=L"\langle N\rangle",
        xlims=xlims,xticks=xticks,ylims=(0,7),yticks=0:1:7,
        title=title_text,framestyle=:box,margin=6Plots.mm)
    for (lower,upper) in bands
        vspan!(p,[lower,upper];color=RGB(.72,.72,.72),alpha=.32,label="")
        vline!(p,[lower,upper];color=RGB(.42,.42,.42),linestyle=:dot,
            linewidth=.8,label="")
    end
    plot!(p,z,mean;yerror=sem,color=:black,marker=:circle,
        markercolor=:white,markerstrokecolor=:black,markersize=4,
        errorbar_color=RGB(.30,.30,.30),errorbar_linewidth=1.0,
        label=L"\langle N\rangle\pm\mathrm{SEM}",legend=:bottomright,
    )
    annotate!(p,xlims[1]+0.035*(xlims[2]-xlims[1]),6.68,text(letter,17,:bold,:left))
    p
end

p1=nz_panel(cb7.z,cb7.mean,cb7_sem,"CB7","a",(-10,10),-10:5:10,cb7_bands)
p2=nz_panel(oa_curve.z,oa_curve.mean,oa_curve.sem,"OA","b",(-3,11),-2:2:10,oa_bands)
p3=nz_panel(cd.z,cd.mean,cd_sem,"restrained β-CD","c",(-10,10),-10:5:10,cd_bands)
nz=plot(p1,p2,p3;layout=(1,3),size=(1800,560),dpi=600,
    left_margin=6Plots.mm,right_margin=5Plots.mm,bottom_margin=7Plots.mm,top_margin=5Plots.mm)
for ext in ("png","pdf","svg"); savefig(nz,joinpath(OUTPUT_DIR,"mean_occupancy_vs_z.$ext")); end

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
    records=(("CB7",cb7.z,cb7.mean,cb7.sd,cb7_sem,"cluster bootstrap"),
             ("OA",oa_curve.z,oa_curve.mean,oa_curve.sd,oa_curve.sem,"source std err column"),
             ("restrained beta-CD",cd.z,cd.mean,cd.sd,cd_sem,"cluster bootstrap"))
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
