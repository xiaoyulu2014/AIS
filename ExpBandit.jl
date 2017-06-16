
using PyPlot
using JLD
cc=logspace(-2,4,100)

regret=zeros(100,1600-1)

for i=1:100
    tmp=load("/data/ziz/not-backed-up/xlu/results/bandit_326215/$i.jld")["dict"]
    index = find(cc.==tmp["hyperparameter"])[1]
    regret[index,:] = regret[index,:] + tmp["out"][2]'
end

regret=regret./10


    
    
    
    
    
