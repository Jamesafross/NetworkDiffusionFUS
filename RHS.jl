function NetworkFlow(du,u,h,p,t)
    @unpack SC,lags,N = p
    for i = 1:N
        Pr = 0.0
        for j = 1:N
            if SC[j,i] > 0.0
                Pr += proba(SC[j,i],N,h(p,t-lags[j,i],idxs=j))
            end
        end

        du[i] = (-1*u[i] + Pr + stim(10,30,39,i,t))

        if u[i] > 1
            u[i] = 1
        end
       
        
    end

end

function proba(W,N,u)
    return 8*W*u/N
end

function stim(tstart,tend,node,i,t)
    if  tstart < t < tend && i == node
        return 1.
    else
        return 0.
    end
end



