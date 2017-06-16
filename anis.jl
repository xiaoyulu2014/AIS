
ENV["JULIA_PKGDIR"] = "/data/ziz/xlu/.julia"  #set package directory

using JLD    #run Pkg.add("JLD") if you do not have relevant packages.
using Iterators # package that can be used to generates hyperparameter grid

# hyperparameter values in sweep
cc = repeat(logspace(-3,3,100),inner=[10],outer=[1])
hyperparams = Iterators.product(cc)
numhyperparams = length(hyperparams)


job_id = parse(Int32,ENV["SLURM_ARRAY_JOB_ID"])    #parse environment variables in Julia
task_id = parse(Int32,ENV["SLURM_ARRAY_TASK_ID"])

if task_id > numhyperparams
  error("Number of jobs in array is more than number of hyperparameters")
end

hyperparam = Iterators.nth(hyperparams,task_id)
@show job_id
@show task_id
@show hyperparam



using StatsBase
using Distributions
using PyPlot


function bandit(C)
	function boost_func(i,nn,C) 
		return sqrt(C*log(i+1)/(nn))   ##choose C1
	end
	
	NN = 100
	num_iter = 16*NN
	n = 100
	sigma = boost_func(1,1,c)*ones(n)
	Y0 = 0 #0.1
	Y = Y0*ones(n)
	R = Float64[]; 
	q = ones(n)*(1/n);
	nn = ones(n);
	Z = 50.50 
	yy=collect((1:(n))/n)


	for i=1:num_iter
		k = Distributions.sample(1:n,WeightVec(q))
		sigma[k] = boost_func(i,nn[k],c)
		w = rand(Bernoulli(0.01))*k*n/n
		Y[k] = (Y[k]*(nn[k]-1) + w)/nn[k] 
		q = (Y+sigma)./sum(Y+sigma)	
		zz = collect(1:n)/Z
		push!(R,sum(yy.*log(yy./(Z*q))))
		nn[k] += 1
	end
	return(q,R[2:end])
end


out=bandit_func(hyperparam)
dict=Dict("out"=>out,"hyperparameter"=> hyperparam)
save("/data/localhost/not-backed-up/xlu/jobname_$(job_id)_$(task_id).jld","dict",dict)   




















