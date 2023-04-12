function test_calculate_bernoulli_gittins_multiple()
    gi = calculate_bernoulli_gittins(
        alpha=Int64(1),
        beta=Int64(1),
        gamma=0.9,
        pulls=Int64(3),
        tol=1e-4,
        N=Int64(200),
    )

    expected_gi = [
        [0.7029 0.8001 0.8452]; 
        [0.5001 0.6346 missing]; 
        [0.3796 missing missing]
    ]

    for i in 1:3
        for j in 1:3
            if i + j <= 3
                diff = abs(expected_gi[i, j] - gi[i, j])
                if diff > 0.01
                    return false
                end
            end
        end
    end

    return true
end

function test_calculate_bernoulli_gittins_single()
    gi = calculate_bernoulli_gittins(
        alpha=Int64(1),
        beta=Int64(1),
        gamma=0.9,
        pulls=Int64(0),
        tol=1e-4,
        N=Int64(200),
    )

    expected_gi = 0.7029

    return abs(expected_gi - gi) < 0.01
end

@test test_calculate_bernoulli_gittins_multiple() == true
@test test_calculate_bernoulli_gittins_single() == true