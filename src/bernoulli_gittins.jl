## Adapted from:
## James Edwards 19: 
## Practical Calculation of Gittins Indices for Multi-armed Bandits

"""
Calculates bernoulli gittins indices for a given arm prior,
	parameterized as Beta(alpha, beta). To calculate all the
	possible states from a given starting prior and a number of pulls,
	optionally pass the pulls argument. the function will then return 
	a triangular matrix where gi_matrix[i, j] = gittins index of 
	beta(alpha + i - 1, beta + j - 1) for i, j >= 1 and i + j <= pulls.

	alpha: current successes
	beta: current failures
	gamma: discount factor
	pulls: optional arg that tells the function how many priors to calculate

Optional args controlling accuracy of approximation:
	N: how far to hallucinate into the future
	tol: fineness of approximation given all of the above.
"""
function calculate_bernoulli_gittins(;alpha, beta, gamma, pulls=0, tol=1e-3, N=200)
    if pulls == 0
        return bmab_gi_ab(
            alpha=alpha, 
            beta=beta, 
            gamma=gamma, 
            tol=tol, 
            N=N,
        )
    else
        return bmab_gi_multiple_ab(
            alpha_start=alpha,
            beta_start=beta,
            gamma=gamma,
            tol=tol,
            N=N,
            num_actions=pulls,
        )
    end
end



function bmab_gi_multiple_ab(;
	alpha_start::Int64, 
	beta_start::Int64, 
	gamma::Float64, 
	N::Int64, 
	num_actions::Int64, 
	tol::Float64,
)	
	GI = Array{Union{Float64, Missing}, 2}(missing, num_actions, num_actions)
  	alpha_range = collect(alpha_start:(alpha_start + num_actions - 1))
	beta_range = collect(beta_start:(beta_start + num_actions - 1))
  	mu = alpha_start ./ (alpha_start .+ beta_range)
  	lb_vec = bmab_kgi(alpha_start, alpha_start .+ beta_range, gamma)
	
  	@progress for a in 1:num_actions
		ub = 1
		for b in 1:(num_actions - a + 1)
			GI[b, a] = bmab_gi_ab(
				alpha = alpha_range[a], 
				beta = beta_range[b], 
				gamma = gamma,
				tol = tol, 
				N = N, 
				lb = lb_vec[b], 
				ub = ub,
			)
		  	ub = GI[b, a]
		end
		lb_vec = GI[:, a]
	end
	return GI
end

function bmab_gi_ab(;
	alpha::Int64, 
	beta::Int64, 
	gamma::Float64,
	tol::Float64, 
	N::Int64, 
	lb::Union{Float64, Int64, Nothing} = nothing, 
	ub::Union{Float64, Int64, Nothing} = nothing,
)
	return bmab_gi(
		Sigma = alpha, 
		n = alpha + beta, 
		gamma = gamma, 
		tol = tol,
		N = N,
		lb = lb, 
		ub = ub,
	)
end

function calibrate_arm(f, lb, ub, tol, other_args...)
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

function bmab_giplus_value(lambda, Sigma::Int64, n::Int64, gamma::Float64)
	mu = Sigma / n
	mu_success = (Sigma + 1) / (n + 1)
	H = gamma / (1 - gamma)
	
	function f_continue(x) 		
		return pdf(Beta(Sigma + 1, n - Sigma), x) * x
	end	

	integral, err = quadgk(f_continue, lambda, 1)
	value_success = H * (integral + lambda * cdf(Beta(Sigma + 1, n - Sigma), lambda))
	
	value_fail = H * lambda
	
	return mu + mu * value_success + (1 - mu) * value_fail - lambda / (1 - gamma)
end

function bmab_giplus(
	Sigma::Int64, 
	n::Int64, 
	gamma::Float64, 
	tol::Float64, 
	upper::Bool = false,
)
	interval = calibrate_arm(
		bmab_giplus_value, 
		Sigma / n, 
		1, 
		tol, 
		Sigma, 
		n, 
		gamma,
	)
	
  	if upper
		return interval[2]
	end
	
	return mean(interval)
end

function bmab_kgi(
	Sigma::Union{Array{Int64}, Int64}, 
	n::Union{Array{Int64}, Int64},
	gamma::Float64,
)
  	mu = Sigma ./ n
  	H = gamma / (1 - gamma)
	
	return (mu .+ H .* mu .* (Sigma .+ 1) ./ (n .+ 1)) ./ (1 .+ H .* mu)
  	
end

function bmab_gi_value(lambda, Sigma, n, gamma, N)
	h = N + 1
	
	n_vec = collect(n : (n + N))
	s_vec = collect(Sigma : (Sigma + N))
	
	mu = s_vec ./ transpose(n_vec)
	
	value_mat = Array{Union{Missing, Float64}, 2}(missing, h, h)
	
	# Values of end states
	value_mat[:, h] = max.(mu[:, h], lambda) .* gamma .^ N / (1 - gamma)
	safe_reward = lambda * gamma .^ ((1:N) .- 1) ./ (1 - gamma)
	
	# Run DP to get values of other states
	for i in N:-1:1
		j = i + 1
		  
		risky_reward = mu[1 : i, i] .* (gamma ^ (i - 1) .+ value_mat[2 : j, j]) .+ (1 .- mu[1 : i, i]) .* value_mat[1 : i, j]
		  
		value_mat[1 : i, i] = max.(risky_reward, safe_reward[i])
  	end
	
	return (value_mat[1, 1] - lambda / (1 - gamma))
end

function bmab_gi(;
	Sigma::Int64, 
	n::Int64, 
	gamma::Float64, 
	N::Int64, 
	tol::Float64,
	lb::Union{Float64, Int64, Nothing} = nothing, 
	ub::Union{Float64, Int64, Nothing} = nothing,
)
	if isnothing(lb)
		lb = bmab_kgi(Sigma, n, gamma)
	end
	
	if isnothing(ub)
    	ub = bmab_giplus(Sigma, n, gamma, tol, true)
	end
	
  	return mean(calibrate_arm(bmab_gi_value, lb, ub, tol, Sigma, n, gamma, N))
end

