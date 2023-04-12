var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = GittinsIndices","category":"page"},{"location":"#GittinsIndices","page":"Home","title":"GittinsIndices","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for GittinsIndices.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [GittinsIndices]","category":"page"},{"location":"#GittinsIndices.calculate_bernoulli_gittins-Tuple{}","page":"Home","title":"GittinsIndices.calculate_bernoulli_gittins","text":"Calculates bernoulli gittins indices for a given arm prior, \tparameterized as Beta(alpha, beta). To calculate all the \tpossible states from a given starting prior and a number of pulls, \toptionally pass the pulls argument. the function will then return  \ta triangular matrix where gi_matrix[i, j] = gittins index of  \tbeta(alpha + i - 1, beta + j - 1) for i, j >= 1 and i + j <= pulls.\n\nalpha: current successes\nbeta: current failures\ngamma: discount factor\npulls: optional arg that tells the function how many priors to calculate\n\nOptional args controlling accuracy of approximation: \tN: how far to hallucinate into the future \ttol: fineness of approximation given all of the above.\n\n\n\n\n\n","category":"method"},{"location":"#GittinsIndices.calculate_discrete_gittins-Tuple{Int64, Matrix{Float64}, Vector{Float64}, Int64, Float64}","page":"Home","title":"GittinsIndices.calculate_discrete_gittins","text":"Given the specified discrete stochastic process, compute the gittins  index of the starting state. \n\nm: number of states in stochastic process P: transition matrix between states r: reward vector with each entry specifying reward for that state initial_state: starting state, whose gittins index will be returned gamma: discount factor\n\n\n\n\n\n","category":"method"},{"location":"#GittinsIndices.calculate_gaussian_gittins-Tuple{}","page":"Home","title":"GittinsIndices.calculate_gaussian_gittins","text":"Calculate the gittins index for a given gaussian belief. Takes the following required parameters: \tmu: mean of gaussian \ttau: precision (inverse of variance) \tn: number of times the arm has been pulled \tgamma: discount factor. Also optionally parameters concerning the accuracy of the returned approximation: \txi: controls the range of the discretization, i.e. many standard deviations away from the mean to consider discretizing \tdelta: the step size in the discretization of the gaussian \tN: horizon size, i.e. how many states you hallucinate into the future \ttol: the fineness of the calculation, given the above params \n\n\n\n\n\n","category":"method"}]
}
