function dNetwork(dx,x,h,p,t)
    β,α,τ = p
    @unpack W,V,N,δ = ND

    for i = 1:N
        
        d = 0.0
        for j = 1:N
            if W[i,j] .> 0 
            d += W[i,j]*h(p,t-lags[i,j];idxs=j)
            end
        end
        dx[i] = (β*d*x[i+N] - α*x[i])*τ
        dx[i+N] = (-β*d*x[i+N])*τ
        dx[i+2N] = (α*x[i])*τ
    end
  
end


function s(x)
    if x > 0.
        return x
    else 
        return 0.
    end
end

function run_diffusion_SIR(β,α,th,h,τ,tspan,x0)
    p = β,α,τ
    alg= MethodOfSteps(BS3())
    prob = DDEProblem(dNetwork,x0,h,tspan,p)
    sol = solve(prob,alg,saveat=0.1)

    u = sol[1:N,:]
    
    m = mean(sol[1:N,:])
    st = std(sol[1:N,:])
    #th = m +2*st
    u[u .< th] .= 0
    u[u .> th] .= 1

    return sol,u

end

mutable struct network_data
    W::Matrix{Real}
    V::Vector{Real}
    N::Real
    δ::Vector{Real}
end


mutable struct DiffParams
    β::Real
    α::Real
    th::Real
    cor::Real
end

function make_struct_context(SC,k,dist)
    N=size(dist,1)
    SCo = zeros(N,N)
    for i = 1:N
        for j = 1:N
            SCo[i,j] = (SC[k,i]*SC[k,j])/(dist[k,i]*SC[k,j] + dist[k,j]*SC[k,i])
        end
    end
    SCo[isnan.(SCo)]
    return SCo


end