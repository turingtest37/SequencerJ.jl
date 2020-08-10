using SequencerJulia
using StatsBase
using Test

@testset "SequencerJulia.jl" begin
    # Write your tests here.
    using Images
    imgsmall = load(joinpath(@__DIR__,"..","resources","colony.png"))
    SMALL = convert(Matrix{Float32}, float32.(Gray.(imgsmall)))

    imgbig = load(joinpath(@__DIR__,"..","resources","h+v.jpg"))
    BIG = convert(Matrix{Float32}, float32.(Gray.(imgbig)))


    include("algotests.jl")

    @testset "real image does not fail." begin
        @test sequence(SMALL; scales=[2]) !== nothing
    end

    @testset "small image does not fail." begin
        @test sequence(SMALL; scales=[2]) !== nothing
    end

    @testset "small image, multiple scales." begin
        @test sequence(SMALL; scales=[1,2,8]) != nothing
    end

    @testset "large image does not fail." begin
        @test sequence(BIG; scales=[2]) !== nothing
    end



end
