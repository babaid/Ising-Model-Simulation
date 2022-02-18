using Revise
using Test

include("../src/models.jl")
using Ising
@testset "Hamiltonians" begin
    
    l = Lattice((100, 100), false)
    @test isreal(H_B(l, 1, 1))
    @test isreal(H_0(l, 1, 1))
            
    end
        
    
    