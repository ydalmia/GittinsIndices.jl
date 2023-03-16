using GittinsIndices
using Test


@testset "all_gittins" begin
    @testset "bernoulli_gittins" begin
        include("test_bernoulli_gittins.jl")
    end

    @testset "gaussian_gittins" begin
        include("test_gaussian_gittins.jl")
    end
    
    @testset "discrete_gittins" begin
        include("test_discrete_gittins.jl")
    end 
end