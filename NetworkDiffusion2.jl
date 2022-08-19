using CircularArrayBuffers

CircularArrayBuffer(zeros(N,100))
c=13
SC,dist,lags,N,FC,missingROIs = networksetup(c;digits=3,nSC=2,nFC=1,N=140,normalise=false)
controlR2 = load("$WORKDIR/controlR2.jld","controlR2")
function maxk(a, k)
    b = partialsortperm(a, 1:k, rev=true)
    return collect(zip(b, a[b]))
end

mutable struct PARS
    SC::Matrix{Float64}
    lags::Matrix{Float64}
    N::Real
end



Δt = 0.1

#1t = 1ms
timespan = collect(0:0.1:10)


x0 = zeros(N)

for i = 1:length(timespan)
    
end


