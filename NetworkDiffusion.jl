using JLD,Plots,DifferentialEquations,Parameters,DelimitedFiles,Parameters,Statistics,StatsBase,Graphs,SimpleWeightedGraphs
WORKDIR = "$(homedir())/NetworkDiffusionFUS"
include("NetworkSetup.jl")
include("RHS.jl")
include("NetworkDiffusionFunctions.jl")
V = readdlm("$WORKDIR/Paul.5z1_ROI_size.txt")[2,4:2:end]
c = 13
SC,dist,lags,N = networksetup(c;digits=3,nSC=2,nFC=1,N=140,normalise=false)
lags[lags .> 0.0] = lags[lags .> 0.0] .+ 10.
ROIsize1 = readdlm("$WORKDIR/Paul.5z1_ROI_size.txt")[2,4:2:end]
const ROIsize = 1 ./(ROIsize1/(maximum(ROIsize1)))




# load real data

ROIs = [94,95,96,132,133,134,16,17,
        18,19,20,89,90,120,121,122,123,124,125,39,104,9,106,74,75]

     
start_time = [60,60,60, 70, 70, 70,
        47,47,47,41,41,41,41,64,64,64,64,64,64,44,50,60,50,65,65]


δ = zeros(N)
for i = 1:N
    δ[i] = sum(SC[i,:])
end


S = SimpleWeightedGraph( 1 ./ SC)
D = SimpleWeightedGraph(dist)


sp_D = johnson_shortest_paths(D).dists
sp_S = johnson_shortest_paths(S).dists
#size_weighted_SC!(SC,ROIsize1,N)
SCo = make_struct_context(SC,39,dist)









nodes = collect(1:1:140)

spD_lhc = sp_D[39,:]
spS_lhc = sp_S[39,:]
spSD_lhc = spD_lhc.*spS_lhc

SCo = make_struct_context(SC,39,dist)


const ND = network_data(SC./maximum(SC),V,N,δ)



    h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 0.0 : zeros(N)


    β = 10.0
    α = 1.2
    th = 0.05
    τ = 0.0008

   
    tspan = (0.,80*100)
    x0 = zeros(3N)
    x0[N+1:2*N] .= 1.0
    NS = 39

    println("stim node = ", NS)
    x0[NS] = 1.0
    x0[N+NS] = 0.0 + (1-x0[NS])
    x0[2N+NS] = x0[NS]

    x0[13] = 1.0
    x0[N+13] = 0.0 + (1-x0[NS])
    x0[2N+13] = x0[NS]

    x0[15] = 1.0
    x0[N+15] = 0.0 + (1-x0[NS])
    x0[2N+15] = x0[NS]

    x0[19] = 1.0
    x0[N+19] = 0.0 + (1-x0[NS])
    x0[2N+19] = x0[NS]

    x0[36] = 1.0
    x0[N+36] = 0.0 + (1-x0[NS])
    x0[2N+36] = x0[NS]






    sol,u = run_diffusion_SIR(β,α,th,h,τ,tspan,x0)

    dists_roi = zeros(length(ROIs))
    str_roi = zeros(length(ROIs))
    simulated_start_time = zeros(length(ROIs))
    for i = 1:length(ROIs)
        try 
            dists_roi[i]=spD_lhc[ROIs[i]]
            str_roi[i]=spSD_lhc[ROIs[i]]
            simulated_start_time[i] = findfirst(u[ROIs[i],:] .== 1)
        catch
            simulated_start_time[i]= 100000
        end
    end
    
    corr2 = corspearman(simulated_start_time,dists_roi)
    corr1 = corspearman(simulated_start_time,start_time)
    corr3 = corspearman(simulated_start_time,str_roi)
    println("corr = ",corr1)
    println("corr with dist = ", corr2)
    println("corr with str = ", corr3)

p1 = heatmap(sol[1:N,1:10:end],xlabel="time",ylabel="ROI")
p2 = heatmap(u[1:N,1:10:end],xlabel="time",ylabel="ROI")

plot(p1,p2)
   