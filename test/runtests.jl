using SequencerJulia
using StatsBase
using Test

@testset "SequencerJulia.jl" begin
    # Write your tests here.

    # L2
    @testset "Squared Euclidean Distance" begin
        objs = [[1,1,1], [2,3,4]]
        expected = [0. 14. ; 14. 0.]
        @test l2dm(objs) == expected
    end

    @testset "l2 simple" begin
        objs = [[1,1,1], [2,3,4]]
        expected = 14.
        @test l2(objs...) == expected
    end

    # KL
    @testset "KL Divergence DM" begin
        a,b = ([1/2, 1/2], [9/10, 1/10])
# pretty slick !!
        expected = reshape(collect(kl(x,y) for x in (a,b), y in (a,b)),2,2)
        @test kldm(a,b) == expected
    end

    @testset "KL Divergence equivalent methods" begin
        objs = [[1/2, 1/2], [9/10, 1/10]]
        @test kl(objs) == kl(objs[1],objs[2])
    end

    @testset "Same as scipy.stats.entropy" begin
        u,v = [1/2, 1/2], [9/10, 1/10]
        expected = 0.5108256237659907
        @test kl(u,v) == expected
    end

    # EMD or 1-d Wasserstein
    @testset "EMD No Fail" begin
        objs = [[1,1,1], [2,3,4]]
        @test EMD(objs[1], objs[2]) != nothing
    end


    # EMD or 1-d Wasserstein
    @testset "EMD with grid" begin
        objs = [[3.4, 3.9, 7.5, 7.8], [4.5, 1.4],
        [1.4, 0.9, 3.1, 7.2], [3.2, 3.5]]
        expected = 4.0781331438047861
        @test EMD(objs...) == expected
    end

    @testset "Energy not null" begin
        objs = [[3.4, 3.9, 7.5, 7.8], [4.5, 1.4],
        [1.4, 0.9, 3.1, 7.2], [3.2, 3.5]]
        expected = 4.0781331438047861
        @test energy(objs...) != nothing
    end


    @testset "Energy no fail" begin
        objs = [[1,1,1], [2,3,4]]
        @test Energy(objs[1], objs[2]) != nothing
    end

end
