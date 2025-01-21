using MethylationReduction
using Test

@testset "MethylationReduction.jl" begin
    X = rand(100,100)
    sort!(X,dims=1) #to give some structure to X
    k = 10
    nmf = NMFCache(X, k; α=0.2);
    solveNMF(nmf, verbose=false, maxiter=200)

    @test residual(nmf) < 0.05
end
