using SequencerJulia
using StatsBase
using Test

@testset "SequencerJulia.jl" begin
    # Write your tests here.
    using Images
    imgsmall = load(joinpath(@__DIR__,"..","resources","colony.png"))
    SMALL = convert(Matrix{Float32}, float32.(Gray.(imgsmall)))

    imgbig = load(joinpath(@__DIR__,"..","resources","Hummingbird!.jpeg"))
    BIG = convert(Matrix{Float32}, float32.(Gray.(imgbig)))


    include("algotests.jl")

    @testset "real image does not fail." begin
        @test sequence(SMALL; scales=[2]) !== nothing
    end

    @testset "small image does not fail." begin
        @test sequence(SMALL; scales=[2]) !== nothing
    end

    @testset "small image, multiple scales." begin
        @test sequence(SMALL; scales=[1,2,8]) !== nothing
    end

    @testset "large image does not fail." begin
        @test sequence(BIG; scales=[2], metrics=(KLD, WASS1D)) !== nothing
    end

    @testset "reorder small image" begin
        using Random, LightGraphs

        Idx = shuffle(axes(SMALL,2))
        imgshuff = SMALL[:,Idx]
        MSTD, ηD, BFSD = sequence(imgshuff; scales=(1), metrics=(KLD, WASS1D))
        res = SMALL[:, collect(vertices(BFSD))]
        @test SMALL == res
    end


end
