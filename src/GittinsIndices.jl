module GittinsIndices

using LinearAlgebra
using Distributions
using JuMP
using HiGHS
using NLsolve
using Statistics
using ProgressLogging
using QuadGK

include("bernoulli_gittins.jl")
include("gaussian_gittins.jl")
include("discrete_gittins.jl")

export calculate_bernoulli_gittins
export calculate_gaussian_gittins
export calculate_discrete_gittins

end # module GittinsIndices
