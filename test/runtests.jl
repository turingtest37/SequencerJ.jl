using SequencerJulia
using StatsBase
using Test

@testset "SequencerJulia.jl" begin
    # Write your tests here.

    @testset "Squared Euclidean Distance" begin
        objs = [[1,1,1], [2,3,4]]
        expected = [0. 14. ; 14. 0.]
        @test l2(objs) == expected
    end

    @testset "KL Divergence" begin
        objs = [[1/2, 1/2], [9/10, 1/10]]
        expected = [0. -0.2938933324510595; 0.8047189562170501 0.]
        @test kl(objs) == expected
    end

    @testset "EMD No Fail" begin
        objs = [[1,1,1], [2,3,4]]
        @test EMD(objs[1], objs[2]) != nothing
    end


    @testset "EMD with grid" begin
        objs = [[3.4, 3.9, 7.5, 7.8], [4.5, 1.4],
        [1.4, 0.9, 3.1, 7.2], [3.2, 3.5]]
        expected = 4.0781331438047861
        @test EMD(objs...) == expected
    end


end
