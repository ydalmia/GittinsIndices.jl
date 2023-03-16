"""
Calculate the gittins index for a given gaussian belief. Takes the following required parameters:
	mu: mean of gaussian
	tau: precision (inverse of variance)
	n: number of times the arm has been pulled
	gamma: discount factor.
Also optionally parameters concerning the accuracy of the returned approximation:
	xi: controls the range of the discretization, i.e. many standard deviations away from the mean to consider discretizing
	delta: the step size in the discretization of the gaussian
	N: horizon size, i.e. how many states you hallucinate into the future
	tol: the fineness of the calculation, given the above params 
"""
function calculate_gaussian_gittins(;mu, tau, n, gamma, tol, N=50, xi=3, delta=0.02)
    return nmab_gi(
        Sigma = mu * n, 
        n = n, 
        gamma = gamma, 
        tau = tau, 
        tol = tol, 
        N = N, 
        xi = xi, 
        delta = delta,
    )
end


function calibrate_arm(f::Function, lb::Union{Int64, Float64}, ub::Float64, tol::Float64, other_args...)
	while (ub - lb) > tol
		lambda = lb + (ub - lb) / 2
		if f(lambda, other_args...) > 0
			lb = lambda
		else
			ub = lambda
		end
	end
	return [lb, ub]
end

function nmab_giplus_value(lambda, mu, n, gamma::Float64)
	sd = sqrt(1 / n)
	return mu + gamma * pdf(Normal(), lambda / sd) * sd / (1 - gamma * cdf(Normal(), lambda / sd)) - lambda
end

function nmab_giplus(;Sigma::Union{Float64, Int64}, n::Int64, gamma::Float64, tol::Float64, ub = nothing, upper::Bool = false)
	if isnothing(ub)
		ub = gamma / (1 - gamma) / sqrt(n)
	end

	lb = 0
	interval = Sigma / n .+ calibrate_arm(
		nmab_giplus_value, lb, ub, tol, 0, n, gamma
	)
	
	if upper
		return interval[2]
	end
	
  	return mean(interval)
end

function nmab_kgi_value(lambda, mu, n::Int64, gamma::Float64, tau::Int64)
	sigt = sqrt(1 / n - 1 / (n + tau))
	z = (mu - lambda) / sigt
	v = sigt * (z * cdf(Normal(), z) + pdf(Normal(), z))
  	return v * gamma / (1 - gamma) + mu - lambda
end

function nmab_kgi(;Sigma::Union{Float64, Int64}, n::Int64, gamma::Float64, tau::Int64, tol::Float64, ub = nothing, lower=false)
	if isnothing(ub)
		ub = gamma / (1 - gamma) / sqrt(n)
	end

	lb = 0
	interval = Sigma / n .+ calibrate_arm(
		nmab_kgi_value, lb, ub, tol, 0, n, gamma, tau,
	)
	
	if lower
		return(interval[1])
	end
	
	return mean(interval)
end

function nmab_risky_reward(;
	mu::Float64, 
	y_lo_scaled, 
	y_hi_scaled, 
	tn_scaled, 
	tau::Int64, 
	s, 
	value_vec::Vector{Float64}, 
	discount::Union{Float64, Int64},
)
  yhi = y_hi_scaled .- mu .* tn_scaled
  ylo = y_lo_scaled .- mu .* tn_scaled
  p = cdf.(Normal(mu, s), yhi) - cdf.(Normal(mu, s), ylo)
	
  return discount * mu .+ (
	  cdf.(Normal(mu, s), ylo[1]) * value_vec[1] .+ (1 - cdf.(Normal(mu, s), yhi[end])) * value_vec[end] .+ sum(p .* value_vec)
  )
end

function nmab_gi_value(lambda::Float64, n::Int64, gamma::Float64, tau::Int64, N::Int64, xi::Union{Int64, Float64}, delta::Float64)
	extra = 1 # number of extra xi used for new states
	h = N + 1
	delta = delta / sqrt(n) # adjust delta so number of states is constant with n
	mu_range = 0 : delta : ((xi + extra) * sqrt(1 / n))
	mu_range = collect(mu_range)
	lr = length(mu_range)
	lr2 = 0 : delta : (xi * sqrt(1 / n))
	lr2 = length(collect(lr2))
	
	value = zeros(Float64, h, lr)
	
	# Value of end states (at stage N)
	value[h, :] = max.(mu_range, lambda) .* (gamma ^ N / (1 - gamma))
	rr = gamma .^ collect((0 : N)) ./ (1 - gamma)
	value[:, lr2 : lr] = rr * transpose(mu_range[lr2 : lr])
  
	lo = mu_range .- delta / 2 # lower end of discrete approximation to mu
  	hi = mu_range .+ delta / 2
	for j in N : -1 : 2
		t = j - 1
		tn = n + tau * (j - 1)
		# the next 3 variables are used for speed-up only
		y_hi_scaled = hi * (tn + tau) / tau
		y_lo_scaled = lo * (tn + tau) / tau
		tn_scaled = tn / tau
		s = sqrt(1 / tn + 1 / tau) #sd of y
		discount = gamma ^ t
		safe_reward = lambda * discount / ( 1 - gamma)
		value_vec = value[j + 1, :]
		for i in (lr2 - 1): -1 : 1
			risky_reward = nmab_risky_reward(
				mu = mu_range[i], 
				y_lo_scaled = y_lo_scaled, 
				y_hi_scaled = y_hi_scaled, 
				tn_scaled = tn_scaled, 
				tau = tau, 
				s = s, 
				value_vec = value_vec, 
				discount = discount,
			)
			if risky_reward > safe_reward
				value[j, i] = risky_reward
			else
				value[j, 1 : i] .= safe_reward
				break
			end
		end
		# If risky arm preferred in [j, 1] then it will be preferred in starting state
		if value[j, 1] > safe_reward
			return (value[j, 1] - safe_reward)
		end
	end
	
	# Value of risky arm in starting state at time 0
	s = sqrt(1 / n + 1 / tau)
	value_vec = value[2, :]
	risky_reward = nmab_risky_reward(
		mu = mu_range[1], 
		y_lo_scaled = y_lo_scaled = lo * (n + tau) / tau, 
		y_hi_scaled = y_hi_scaled = hi * (n + tau) / tau,
		tn_scaled = tn_scaled = n / tau, 
		tau = tau, 
		s = s, 
		value_vec = value_vec, 
		discount = 1,
	)
	return (risky_reward - lambda / (1 - gamma))
end

function nmab_gi(;Sigma::Union{Float64, Int64}, n::Int64, gamma::Float64, tau::Int64, tol::Float64, N::Int64, xi::Union{Float64, Int64}, delta::Float64, lb=nothing, ub=nothing)
	if isnothing(lb)
		lb = nmab_kgi(
			Sigma = 0, 
			n = n,
			gamma = gamma, 
			tau = tau, 
			tol = tol, 
			ub = ub, 
			lower=true,
		)
	end
	
	if isnothing(ub)
		ub = nmab_giplus(
			Sigma = 0, 
			n = n, 
			gamma = gamma, 
			tol = tol, 
			ub = ub, 
			upper=true,
		)
	end
	
	interval = Sigma / n .+ calibrate_arm(
		nmab_gi_value, lb, ub, tol, n, gamma, tau, N, xi, delta,
	)
	
	return mean(interval)
end