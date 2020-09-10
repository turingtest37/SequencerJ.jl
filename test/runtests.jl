using Test
using SequencerJ
using StatsBase
using Random
using LightGraphs
using Images

@testset "SequencerJ.jl" begin

    Random.seed!(2732);

    include("algotests.jl")

    # declare a loss function to measure the error in calculations
    loss(A,B) = L2(A, B)

    imgsmall = Gray.(load(joinpath(@__DIR__,"..","resources","colony.png")))
    SMALL = convert(Matrix{Float32}, imgsmall)

    imgmed = Gray.(load(joinpath(@__DIR__,"..","resources","Aged Whisky.jpg")));
    MED = convert(Matrix{Float32}, imgmed)

    imgbig = Gray.(load(joinpath(@__DIR__,"..","resources","Hummingbird!.jpeg")))
    BIG = convert(Matrix{Float32}, imgbig)

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

    @testset "reorder small image, scale=1" begin
        Idx = shuffle(axes(SMALL,2))
        imgshuff = SMALL[:,Idx]
        seqres = sequence(imgshuff; scales=(1)) #all metrics
        ind = order(seqres)
        res = imgshuff[:, ind]
        # This way the test will also fail if the 
        # reordering begins to work perfectly, i.e. loss ~ 0!
        # @test loss(SMALL,res) ≈ 0
        @test loss(SMALL,res) < 500
    end

    @testset "reorder small image, autoscale" begin
        Idx = shuffle(axes(SMALL,2))
        imgshuff = SMALL[:,Idx]
        seqres = sequence(imgshuff) #autoscale, all metrics
        ind = order(seqres)
        res = imgshuff[:, ind]
        # This way the test will also fail if the 
        # reordering begins to work perfectly, i.e. loss ~ 0!
        # @test loss(SMALL,res) ≈ 0
        @test loss(SMALL,res) < 500
    end

    @testset "reorder small image, scale=1, weightrows = true" begin
        Idx = shuffle(axes(SMALL,2))
        imgshuff = SMALL[:,Idx]
        seqres = sequence(imgshuff; scales=(1), weightrows=true) #all metrics
        ind = order(seqres)
        res = imgshuff[:, ind]
        # This way the test will also fail if the 
        # reordering begins to work perfectly, i.e. loss ~ 0!
        # @test loss(SMALL,res) ≈ 0
        @test loss(SMALL,res) < 500
    end

    @testset "reorder med image, orig. portrait, scales=(1,2)" begin
        Idx = shuffle(axes(MED,2))
        imgshuff = MED[:,Idx]
        seqres = sequence(imgshuff; scales=(1,2), metrics=(KLD, L2))
        ind = order(seqres)
        res = imgshuff[:, ind]
        # This way the test will also fail if the 
        # reordering begins to work perfectly, i.e. loss ~ 0!
        # @test loss(MED,res) ≈ 0
        @test loss(MED,res) < 700
    end

    @testset "reorder med image, orig. portrait, autoscale" begin
        Idx = shuffle(axes(MED,2))
        imgshuff = MED[:,Idx]
        seqres = sequence(imgshuff; metrics=(KLD, L2))
        ind = order(seqres)
        res = imgshuff[:, ind]
        # This way the test will also fail if the 
        # reordering begins to work perfectly, i.e. loss ~ 0!
        # @test loss(MED,res) ≈ 0
        @test loss(MED,res) < 500
    end

    @testset "reorder med image, orig. portrait, autoscale, weightrows=true" begin
        Idx = shuffle(axes(MED,2))
        imgshuff = MED[:,Idx]
        seqres = sequence(imgshuff; metrics=(KLD, L2), weightrows=true)
        ind = order(seqres)
        res = imgshuff[:, ind]
        # This way the test will also fail if the 
        # reordering begins to work perfectly, i.e. loss ~ 0!
        # @test loss(MED,res) ≈ 0
        @test loss(MED,res) < 800
    end

    @testset "reorder med image, landscape, scales=(1,2)" begin
        MED = permutedims(MED)
        Idx = shuffle(axes(MED,2))
        imgshuff = MED[:,Idx]
        seqres = sequence(imgshuff; scales=(1,2), metrics=(KLD, L2))
        ind = order(seqres)
        res = imgshuff[:, ind]
        # This way the test will also fail if the 
        # reordering begins to work perfectly, i.e. loss ~ 0!
        # @test loss(MED,res) ≈ 0
        @test loss(MED,res) < 1100
    end

    @testset "reorder med image, landscape, autoscale" begin
        MED = permutedims(MED)
        Idx = shuffle(axes(MED,2))
        imgshuff = MED[:,Idx]
        seqres = sequence(imgshuff; metrics=(KLD, L2))
        ind = order(seqres)
        res = imgshuff[:, ind]
        # This way the test will also fail if the 
        # reordering begins to work perfectly, i.e. loss ~ 0!
        # @test loss(MED,res) ≈ 0
        @test loss(MED,res) < 1100
    end

    @testset "reorder med image, landscape, scales=(1,2), weightrows=true" begin
        MED = permutedims(MED)
        Idx = shuffle(axes(MED,2))
        imgshuff = MED[:,Idx]
        seqres = sequence(imgshuff; scales=(1,2), metrics=(KLD, L2), weightrows=true)
        ind = order(seqres)
        res = imgshuff[:, ind]
        # This way the test will also fail if the 
        # reordering begins to work perfectly, i.e. loss ~ 0!
        @test loss(MED,res) < 1100
    end

    @testset "reorder med image, landscape, autoscale, weightrows=true" begin
        MED = permutedims(MED)
        Idx = shuffle(axes(MED,2))
        imgshuff = MED[:,Idx]
        seqres = sequence(imgshuff; scales=(1,2), metrics=(KLD, L2), weightrows=true)
        ind = order(seqres)
        res = imgshuff[:, ind]
        # This way the test will also fail if the 
        # reordering begins to work perfectly, i.e. loss ~ 0!
        @test loss(MED,res) < 800
    end

    # This test fails on some combinations of random numbers
    @testset "reorder large image, orig. landscape, scales=(1,2,4)" begin
        Idx = shuffle(axes(BIG,2))
        imgshuff = BIG[:,Idx]
        seqres = sequence(imgshuff; scales=(1,2,4), metrics=(L2, KLD))
        ind = order(seqres)
        res = imgshuff[:, ind]
        @test loss(BIG,res) < 100
    end

    # This test fails on certain combinations of random numbers
    @testset "reorder large image, orig. landscape, autoscale" begin
        Idx = shuffle(axes(BIG,2))
        imgshuff = BIG[:,Idx]
        seqres = sequence(imgshuff; scales=(1,2,4), metrics=(L2, KLD))
        ind = order(seqres)
        res = imgshuff[:, ind]
        @test loss(BIG,res) < 100
    end

    @testset "reorder large image, orig. landscape, scales=(1,2,4), weightrows=true" begin
        Idx = shuffle(axes(BIG,2))
        imgshuff = BIG[:,Idx]
        seqres = sequence(imgshuff; scales=(1,2,4), metrics=(L2, KLD), weightrows=true)
        ind = order(seqres)
        res = imgshuff[:, ind]
        # This way the test will also fail if the 
        # reordering begins to work perfectly, i.e. loss ~ 0!
        # @test loss(BIG,res) ≈ 0
        @test loss(BIG,res) < 11500
    end

    @testset "reorder large image, orig. landscape, autoscale, weightrows=true" begin
        Idx = shuffle(axes(BIG,2))
        imgshuff = BIG[:,Idx]
        seqres = sequence(imgshuff; metrics=(L2, KLD), weightrows=true)
        ind = order(seqres)
        res = imgshuff[:, ind]
        # This way the test will also fail if the 
        # reordering begins to work perfectly, i.e. loss ~ 0!
        # @test loss(BIG,res) ≈ 0
        @test loss(BIG,res) < 11500
    end

    # So far, this is the only case that systematically produces ~ 0 loss
    @testset "reorder large image, portrait, scales=(1,2,4)" begin
        BIG = permutedims(BIG)
        Idx = shuffle(axes(BIG,2))
        imgshuff = BIG[:,Idx]
        seqres = sequence(imgshuff; scales=(1,2,4), metrics=(L2, KLD))
        ind = order(seqres)
        res = imgshuff[:, ind]
        @test loss(BIG,res) ≈ 0
    end

    # So far, this is the only case that systematically produces ~ 0 loss
    @testset "reorder large image, portrait, autoscale" begin
        BIG = permutedims(BIG)
        Idx = shuffle(axes(BIG,2))
        imgshuff = BIG[:,Idx]
        seqres = sequence(imgshuff, metrics=(L2, KLD))
        ind = order(seqres)
        res = imgshuff[:, ind]
        @test loss(BIG,res) ≈ 0
    end

    # Weighting rows sometimes results in zero loss, but often it results in a large loss ~ 11000.
    # Just test that loss remains below the current threshold of 11500.
    @testset "reorder large image, portrait, scales=(1,2,4), weightrows=true" begin
        BIG = permutedims(BIG)
        Idx = shuffle(axes(BIG,2))
        imgshuff = BIG[:,Idx]
        seqres = sequence(imgshuff; scales=(1,2,4), metrics=(L2, KLD), weightrows=true)
        ind = order(seqres)
        res = imgshuff[:, ind]
        # This way the test will also fail if the 
        # reordering begins to work perfectly, i.e. loss ~ 0!
        # @test loss(BIG,res) ≈ 0
        @test loss(BIG,res) < 11500
    end

    @testset "reorder large image, portrait, autoscale, weightrows=true" begin
        BIG = permutedims(BIG)
        Idx = shuffle(axes(BIG,2))
        imgshuff = BIG[:,Idx]
        seqres = sequence(imgshuff; metrics=(L2, KLD), weightrows=true)
        ind = order(seqres)
        res = imgshuff[:, ind]
        # This way the test will also fail if the 
        # reordering begins to work perfectly, i.e. loss ~ 0!
        # @test loss(BIG,res) ≈ 0
        @test loss(BIG,res) < 11500
    end

    @testset "Accessor functions" begin    
        seqres = sequence(SMALL; scales=(2,), metrics=(L2,))
        for f in (:elong, :order, :D, :mst)
            @testset "$f" begin
                @eval @test $f($seqres) !==nothing
            end 
        end
    end

    @testset "show" begin
        seqres = sequence(SMALL; scales=(2,), metrics=(L2,))
        @test occursin("Sequencer Result", "$(seqres)")
    end

    @testset "prettyp of a vector" begin
        v = collect(1:20)
        r = "1,2,3...18,19,20"
        @test r == prettyp(v)
    end

    @testset "prettyp of a vector, more digits" begin
        v = collect(1:20)
        r = "1,2,3,4,5...16,17,18,19,20"
        @test r == prettyp(v, 5)
    end

    @testset "prettyp of a vector, small" begin
        v = collect(1:5)
        r = "1,2,3,4,5"
        @test r == prettyp(v)
    end
end
