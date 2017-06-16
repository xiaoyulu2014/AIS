#target 
#number of arms is n
@everywhere begin
	using StatsBase
	using Distributions
	using PyPlot
end
@everywhere begin
	n = 100
	pii(x) = x * (0<=x<=1) 
	xx=(1/n):(1/n):1
	yy=collect((1:(n))/n)  ##yy= Y_a
	Z = 50.50    ##Z = (1+n)/2
	q = ones(n)*(1/n)

	#initial proposal
	nn = zeros(n) #ones(n)
	Y0 = 0 #0.1
	Y = Y0*ones(n)
	sigma = 0.1*ones(n) #zeros(n)
end


@everywhere NN = 1000000
@everywhere num_iter = 16*NN
@everywhere boost = true



@everywhere function boost_func(i,k) 
	return sqrt(0.1*log(i+1)/(nn[k]))   ##choose C1
end

@everywhere function bandit(num_iter::Int64, boost::Bool,boost_func,q,nn,Y,sigma)
	R = Float64[]; 
	for i=1:num_iter
		k = Distributions.sample(1:n,WeightVec(q))
		nn[k] += 1
		boost ? sigma[k] = boost_func(i,k) : sigma = zeros(n)
		w = rand(Bernoulli(0.01))*k*n/n
		Y[k] = (Y[k]*(nn[k]-1) + w)/nn[k] 
		q = (Y+sigma)./sum(Y+sigma)	
		zz = collect(1:n)/Z
		push!(R,sum(yy.*log(yy./(Z*q))))
	end
	return(q, R[2:end])
end

@everywhere function Bandit(x)
	return bandit(num_iter,true,boost_func,q,nn,Y,sigma)
end


























###########################################


#=
res = pmap(Bandit, 1:10)

R_mat = Array(Float64, num_iter-1,10)
for i=1:10
	R_mat[:,i] = res[i][2]
end


R7 = mean(R_mat,2)

plot(cumsum(R1),label="sqrt(log(t))/sqrt(n)")
plot(cumsum(R2),label="sqrt(log(t)/n")
plot(cumsum(R3),label="sqrt(log(t)/n^0.25)")
plot(cumsum(R4),label="log(t)/sqrt(n)")
plot(cumsum(R5),label="(t)^0.25/sqrt(n)")
plot(cumsum(R6),label="log(t)/n")
plot(cumsum(R7),label="no boost")
legend(loc="upper left")
title("cumulative KL regret",fontsize=20)
xlabel("number of iteration",fontsize=20)
ylabel("regret",fontsize=20)




regret=["R1" => R1, "R2" => R2 , "R3" => R3,"R4" => R4, "R5" => R5, "R6" => R6]
save("Bandit_KL_regret.jld","regret",regret)




q_mat = Array(Float64, 100,10)
for i=1:10
	q_mat[:,i] = res[i][1]
end
qq = mean(q_mat,2)


for i=1:(n-1)
	plot([(i-1)/n,i/n],[yy[i]/0.5,yy[i]/0.5],color="blue")
	plot([(i-1)/n,i/n],[n*qq[i],n*qq[i]],color="red")
end 
i=n
plot([(i-1)/n,i/n],[yy[i]/0.5,yy[i]/0.5],color="blue",label="target")
plot([(i-1)/n,i/n],[n*qq[i],n*qq[i]],color="red",label="proposal")
xlabel("number of iteration",fontsize=20)
ylabel("regret",fontsize=20)
legend(loc="upper left",fontsize=20)
title("target and proposal densities",fontsize=20)










n = 100
	shape=(1:n)/10
	scale = 10*(1:n)

	pii(x) = x * (0<=x<=1) 
	xx=(1/n):(1/n):1
	yy=[pii(xx[i]) for i=1:length(xx)]
	Z = 5050

	#initial proposal
	q = ones(n)/n
	nn = zeros(n)
	Y0= 0.1
	Y = Y0*ones(n)
	sigma = 0.001*ones(n)

fig, axs = plt[:subplots](4,4, figsize=(15, 6)); suptitle("Every $(NN) iterations, Y0=$(Y0), boost= sqrt(log(t)/n)");


Nat = Array(Float64,num_iter,n)
R = Float64[] #regret
for i=1:num_iter
	k = sample(1:n,WeightVec(q))
	nn[k] += 1
	Nat[i,:] = nn
	boost ? sigma[k] = sqrt(log(i+1)/nn[k]) : sigma = zeros(n)
	w = rand(Bernoulli(0.01))*k[1]
	Y[k] = (Y[k]*(nn[k]-1) + w)/nn[k] 
	q = (Y+sigma)./sum(Y+sigma)	
	zz = collect(1:n)/Z
    w == 0? push!(R,0) : push!(R,sum(zz.*log(zz./q)))
    # w == 0? push!(R,0) : push!(R,(w/(Z*q[k])*log(w/(q[k]*Z)) - w/(q[k]*Z)) +1)
	
	if (rem(i,NN) == 1) 
		ii=convert(Int64,(i-1)/NN+1)

		if (i < num_iter-NN)
			axs[ii][:plot]((1/n):(1/n):1,n*q,label="proposal")
			axs[ii][:plot](xx,yy/0.5,label="target")
			axs[ii][:set_xlabel]("x")
			axs[ii][:legend](loc="upper left")
		end
	end		
end


axs[16][:plot](1:(num_iter-1),cumsum(R[2:end]),label="loss")
axs[16][:set_xlabel]("iteration")
title("KL loss")
fig[:canvas][:draw]() 


figure()

plot((1/n):(1/n):1,n*q,label="proposal")
plot(xx,yy/0.5,label="target")
xlabel("x")
legend(loc="upper left")

#### see how Nat grows with t:
figure();
for i=1:5
	ii = i
	plot(Nat[:,ii],label = "arm $ii")
end
legend(loc="upper left")
title("boost =  sqrt(log(t)/n)")




=#
















