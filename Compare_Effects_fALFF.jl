using JLD,Plots,DifferentialEquations,Parameters,Graphs,SimpleWeightedGraphs,DelimitedFiles,StatsBase
WORKDIR = "$(homedir())/NetworkDiffusionFUS"
include("NetworkSetup.jl")
include("RHS.jl")
c=13
u = load("SIR_thresh_sol.jld","SIR_thresh_sol")
SC,dist,lags,N,FC,missingROIs = networksetup(c;digits=3,nSC=2,nFC=1,N=140,normalise=false)
mutable struct regions
    name
    ROIs
    sizeROI
    start_time
    end_time
    duration
    T
    sp_from_lhc_dist
    sp_from_lhc_struct
    sp_from_lhc_diststruct
    simulated_start_time
end

ROIsize = readdlm("Paul.5z1_ROI_size.txt")[2,4:2:end]

name = ["ITC", "ITC", "ITC", "PONS", "PONS", "PONS","ACC","IPL","IPL","IPL",
        "SPL","SPL","SPL","PMC","PMC","PMC","PMC"]

ROIs = [94,95,96,132,133,134,1,16,17,
        18,13,14,15,19,20,89,90,74,75]

start_time = [46,46,46, 49, 49, 49,55,
        57,57,57,61,61,61,67,67,67,67]

end_time = [52,52,52,69,69,69,62,100,
        100,100,69,69,69,69,69,69,69]

T = [4.968, 4.968, 4.968,5.842,5.842,5.842,-5.09,
        -5.38,-5.38,-5.38,-4.62,-4.62,-4.62,-4.88,-4.88,-4.88,-4.88]



dataAffect = Array{regions}(undef,length(name))


BinGraph = zeros(N,N)
BinGraph .= dist
BinGraph[BinGraph .> 0.0 ] .= 1.
B = SimpleWeightedGraph(BinGraph)
S = SimpleWeightedGraph( 1 ./ SC)
D = SimpleWeightedGraph(dist)
V = SimpleWeightedGraph(dist ./13 )
sp_D = johnson_shortest_paths(D).dists
sp_V = 
sp_S = johnson_shortest_paths(S).dists
sp_B = johnson_shortest_paths(B).dists

sd = dist ./ SC
sd[isnan.(sd)] .= Inf

SD = SimpleWeightedGraph(sd)
sp_SD = johnson_shortest_paths(SD).dists

spD_lhc = sp_D[:,39]
spS_lhc = sp_S[:,39]
spSD_lhc = sp_SD[:,39]

for i = 1:length(name)
    dataAffect[i] = regions(name[i],ROIs[i],ROIsize[ROIs[i]],start_time[i],end_time[i],
    start_time[i]-end_time[i],T[i],spD_lhc[ROIs[i]],spS_lhc[ROIs[i]],spSD_lhc[ROIs[i]],findfirst(u[ROIs[i],:] .== 1))
end

effect = zeros(length(name))
dist_lhc = zeros(length(name))

for i = 1:length(name)
    effect[i] = dataAffect[i].start_time
    dist_lhc[i] = dataAffect[i].sp_from_lhc_dist
end

println("corr = ",corspearman(effect,dist_lhc))
scatter(dist_lhc,effect)
    