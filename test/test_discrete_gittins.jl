function test_discrete_gittins_against_each_other()
    bp = GittinsIndices.BanditProcess(
        Int64(4),
        [
            0.1 0 0.8 0.1; 
            0.5 0 0.1 0.4; 
            0.2 0.6 0 0.2; 
            0 0.8 0 0.2
        ],
        [16.0, 19.0, 30.0, 4.0],
        Int64(2),
        0.75
    )
    sol_chen_katehakis = GittinsIndices.solve(
        GittinsIndices.ChenKatehakisLinearProgramming(bp)
    )
    sol_katehakis_veinott = GittinsIndices.solve(
        GittinsIndices.KatehakisVeinottRestartFormulation(bp)
    )
    sol_katehakis_veinott ≈ sol_chen_katehakis
end

function test_calculate_discrete_gittins()
    gi = calculate_discrete_gittins(
        Int64(4),
        [
            0.1 0 0.8 0.1; 
            0.5 0 0.1 0.4; 
            0.2 0.6 0 0.2; 
            0 0.8 0 0.2
        ],
        [16.0, 19.0, 30.0, 4.0],
        Int64(2),
        0.75,
    )

    gi2 = GittinsIndices.solve(
        GittinsIndices.ChenKatehakisLinearProgramming(
            GittinsIndices.BanditProcess(
                Int64(4),
                [
                    0.1 0 0.8 0.1; 
                    0.5 0 0.1 0.4; 
                    0.2 0.6 0 0.2; 
                    0 0.8 0 0.2
                ],
                [16.0, 19.0, 30.0, 4.0],
                Int64(2),
                0.75
            )
        )
    )
    
    return gi == gi2
end


@test test_discrete_gittins_against_each_other() == true
@test test_calculate_discrete_gittins() == true