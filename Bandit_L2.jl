@everywhere begin
	using StatsBase
	using Distributions
	using PyPlot
end
@everywhere begin
	n = 100
	shape=(1:n)/10
	scale = 10*(1:n)

	pii(x) = x * (0<=x<=1) 
	xx=(1/n):(1/n):1
	yy=[pii(xx[i]) for i=1:length(xx)]
	Z = 5050
	q = ones(n)*(1/n)

	#initial proposal
	nn = zeros(n)
	Y0= 0.1
	Y = Y0*ones(n)
	sigma = 0.01*ones(n)
end


@everywhere NN = 100000
@everywhere num_iter = 16*NN
@everywhere boost = true



@everywhere function boost_func(i,k) 
	return sqrt(log(i+1))/nn[k]*2*sqrt(2)*(n^0.5)
end

@everywhere function bandit(num_iter::Int64, boost::Bool,boost_func,q,nn,Y,sigma)
	R = Float64[] 
	for i=1:num_iter
		k = sample(1:n,WeightVec(q))
		nn[k] += 1
		boost ? sigma[k] = boost_func(i,k) : sigma = zeros(n)
		w = rand(Bernoulli(0.01))*k[1]
		Y[k] = (Y[k]*(nn[k]-1) + w^2)/nn[k] 
		q = sqrt((Y+sigma))./sum(sqrt(Y+sigma))
		#zz = (collect(1:n)).^2*0.01*n
	   # w == 0? push!(R,0) :  push!(R,1/2*(sum(zz./sum(sqrt(zz))^2./q) - 1))
	    zz=(collect(1:n)).^2
		push!(R, sum(zz./q)/sum(sqrt(zz)).^2 - 1)	
  end
	return(q, R[2:end])
end

@everywhere function Bandit(x)
	return bandit(num_iter,true,boost_func,q,nn,Y,sigma)
end

res1 = pmap(Bandit, 1:10)



R_mat = Array(Float64, num_iter-1,10)
for i=1:10
	R_mat[:,i] = res1[i][2]
end


R6 = mean(R_mat,2)


plot(cumsum(R1),label="sqrt(log(t)/n)")
plot(cumsum(R5),label="sqrt(log(t))/n^0.75")
plot(cumsum(R6),label="sqrt(log(t))/n^0.25")
plot(cumsum(R2),label="sqrt(t/n)")
plot(cumsum(R3),label="log(t)/sqrt(n)")
plot(cumsum(R4),label="log(t)/log(n+1)")

legend(loc="upper left")
title("accummulated regrest, average over 10 runs")
xlabel("number of iteration")
ylabel("regret")

plot(log(1:159999),log(cumsum(R1)))





















