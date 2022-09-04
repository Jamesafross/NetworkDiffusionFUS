using JLD,Plots,DifferentialEquations,Parameters,DelimitedFiles,Parameters,Statistics,StatsBase
WORKDIR = "$(homedir())/NetworkDiffusionFUS"
include("NetworkSetup.jl")
include("RHS.jl")
include("NetworkDiffusionFunctions.jl")
V = readdlm("$WORKDIR/Paul.5z1_ROI_size.txt")[2,4:2:end]
c = 1
SC,dist,lags,N,FC,missingROIs = networksetup(c;digits=3,nSC=2,nFC=1,N=140,normalise=false)


# load real data

ROIs = [94,95,96,132,133,134,16,17,
        18,19,20,89,90,120,121,122,123,124,125,39,104,9,106,74,75]

     
start_time = [60,60,60, 70, 70, 70,
        47,47,47,41,41,41,41,64,64,64,64,64,64,44,50,60,50,65,65]


δ = zeros(N)
for i = 1:N
    δ[i] = sum(SC[i,:])
end




WSC = zeros(N,N)

    
for i = 1:N
    for j = 1:N
        WSC[i,j] = (1/V[i])*SC[i,j]*(1/(δ[j]))*V[j]
    end
end


n1 = 10
n2 = 10
n3 = 10
DfP = Array{DiffParams}(undef,n1,n2,n3)
 


const ND = network_data(SC./maximum(SC),V,N,δ)


β_vec = LinRange(1,4,n1)
α_vec = LinRange(0.01,1,n2)
th_vec = LinRange(0.01,0.1,n3)
h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 0.0 : zeros(N)
for ii = 1:n1;for jj = 1:n2; for kk = 1:n3;

    β = β_vec[ii]
    α = α_vec[jj]
    th = th_vec[kk]

   
    tspan = (0.,300.0)
    x0 = zeros(3N)
    x0[N+1:2*N] .= 1
    x0[39] = 1
    x0[N+39] = 0.0
    sol,u = run_diffusion_SIR(β,α,th,tspan,x0)

    simulated_start_time = zeros(length(ROIs))
    for i = 1:length(ROIs)
        try 
            simulated_start_time[i] = findfirst(u[ROIs[i],:] .== 1)
        catch
            simulated_start_time[i]= 10000000
        end
    end

    corr1 = corspearman(simulated_start_time,start_time)
    println("corr = ",corr1)

    DfP[ii,jj,kk] = DiffParams(β,α,th,corr1)
end;end;end



