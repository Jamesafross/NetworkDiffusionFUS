using JLD,Plots,DifferentialEquations,Parameters,Graphs,SimpleWeightedGraphs,StatsBase,DelimitedFiles
WORKDIR = "$(homedir())/NetworkDiffusionFUS"
include("NetworkSetup.jl")

c=13

SC,dist,lags,N = networksetup(c;digits=3,nSC=2,nFC=1,N=140,normalise=false)



D = SimpleWeightedGraph(dist)


sp_D = johnson_shortest_paths(D).dists



