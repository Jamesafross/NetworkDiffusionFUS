using JLD,Plots,DifferentialEquations,Parameters,Graphs,SimpleWeightedGraphs,StatsBase,DelimitedFiles
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
        sp_from_lhc_delay
        simulated_start_time
end



ROIsize = readdlm("Paul.5z1_ROI_size.txt")[2,4:2:end]


name = ["ITC", "ITC", "ITC", "PONS", "PONS", "PONS","IPL","IPL","IPL",
        "PMC","PMC","PMC","PMC","THAL","THAL","THAL","THAL","THAL","THAL","HC","V4","M1","V1","OFC","OFC"]

ROIs = [94,95,96,132,133,134,16,17,
        18,19,20,89,90,120,121,122,123,124,125,39,104,9,106,74,75]

     
start_time = [60,60,60, 70, 70, 70,
        47,47,47,41,41,41,41,64,64,64,64,64,64,44,50,60,50,65,65]


end_time = [71,71,71, 86,86,86,70,
        70,70,55,55,55,55,70,70,70,70,70,70,56,63,67,63,100,100]

T = [5.678,5.678,5.678, 6.512,6.512,6.512,
        -5.924,-5.924,-5.924,-6.482,-6.482,-6.482,-6.482,5.948,5.948,5.948,5.948,5.948,5.948,
        -6.366,-5.959,-4.579,-6.310,5.929,5.929]

K = [75,75,75, 36,36,36,
    98,98,98,1396,1396,1396,1396,43,43,43,43,43,43,
    471,42,34,26,26,26]



dataAffect = Array{regions}(undef,length(name))


BinGraph = zeros(N,N)
BinGraph .= dist
BinGraph[BinGraph .> 0.0 ] .= 1.
B = SimpleWeightedGraph(BinGraph)
V = SimpleWeightedGraph(dist ./13.)
S = SimpleWeightedGraph( 1 ./ SC)
D = SimpleWeightedGraph(dist)
sd = dist ./ SC
sd[isnan.(sd)] .= Inf
SD = SimpleWeightedGraph(sd)

sp_V = johnson_shortest_paths(V).dists
sp_B = johnson_shortest_paths(B).dists
sp_D = johnson_shortest_paths(D).dists
sp_S = johnson_shortest_paths(S).dists
sp_SD = johnson_shortest_paths(SD).dists


delay_const = zeros(N)
delay_const = sp_B[39,:]
delay_const[delay_const .== 1] .= 0 
delay_const[delay_const .== 2] .= 1


spD_lhc = sp_D[:,39]
spS_lhc = sp_S[:,39]
spSD_lhc = sp_SD[:,39]
spV_lhc = sp_V[:,39] .+ delay_const*0.5



for i = 1:length(name)
    dataAffect[i] = regions(name[i],ROIs[i],ROIsize[ROIs[i]],start_time[i],end_time[i],
    start_time[i]- end_time[i],T[i],spD_lhc[ROIs[i]],spS_lhc[ROIs[i]],spSD_lhc[ROIs[i]],spV_lhc[ROIs[i]],findfirst(u[ROIs[i],:] .== 1))
end


effect = zeros(length(name))
dist_lhc = zeros(length(name))

for i = 1:length(name)
    effect[i] = dataAffect[i].start_time
    dist_lhc[i] = dataAffect[i].sp_from_lhc_delay
end

println("corr = ",corspearman(effect,dist_lhc))
scatter(dist_lhc,effect)

    