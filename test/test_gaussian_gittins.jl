function test_gaussian_gittins_single()
    n = 2
    gamma = 0.8

    gi = calculate_gaussian_gittins(
        mu = 0, 
        tau = 1, 
        n = n, 
        gamma = gamma, 
        tol = 1e-4, 
        N = 100,
        xi = 3,
        delta = 0.02,
    )
    
    table_val = n * sqrt(1 - gamma) * gi 
    return abs(table_val - 0.27584) < 0.01
end

@test test_gaussian_gittins_single() == true