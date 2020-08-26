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
    SMALL = convert(Matrix{Float32}, float32.(imgsmall))

    imgmed = Gray.(load(joinpath(@__DIR__,"..","resources","Aged Whisky.jpg")));
    MED = convert(Matrix{Float32}, float32.(imgmed))

    imgbig = Gray.(load(joinpath(@__DIR__,"..","resources","Hummingbird!.jpeg")))
    BIG = convert(Matrix{Float32}, float32.(imgbig))

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
        Idx = shuffle(axes(SMALL,2))
        imgshuff = SMALL[:,Idx]
        seqres = sequence(imgshuff; scales=(1)) #all metrics
        ind = order(seqres)
        res = imgshuff[:, ind]
        # This way the test will also fail if the 
        # reordering begins to work perfectly, i.e. loss ~ 0!
        @test 0.1 < loss(SMALL,res) < 500
    end

    @testset "reorder small image, weightrows = true" begin
        Idx = shuffle(axes(SMALL,2))
        imgshuff = SMALL[:,Idx]
        seqres = sequence(imgshuff; scales=(1), weightrows=true) #all metrics
        ind = order(seqres)
        res = imgshuff[:, ind]
        # This way the test will also fail if the 
        # reordering begins to work perfectly, i.e. loss ~ 0!
        @test 0.1 < loss(SMALL,res) < 500
    end

    @testset "reorder med image, orig. portrait" begin
        Idx = shuffle(axes(MED,2))
        imgshuff = MED[:,Idx]
        seqres = sequence(imgshuff; scales=(1,2), metrics=(KLD, L2))
        ind = order(seqres)
        res = imgshuff[:, ind]
        # This way the test will also fail if the 
        # reordering begins to work perfectly, i.e. loss ~ 0!
        @test 0.1 < loss(MED,res) < 800
    end

    @testset "reorder med image, orig. portrait, weightrows=true" begin
        Idx = shuffle(axes(MED,2))
        imgshuff = MED[:,Idx]
        seqres = sequence(imgshuff; scales=(1,2), metrics=(KLD, L2), weightrows=true)
        ind = order(seqres)
        res = imgshuff[:, ind]
        # This way the test will also fail if the 
        # reordering begins to work perfectly, i.e. loss ~ 0!
        @test 0.1 < loss(MED,res) < 800
    end

    @testset "reorder med image, landscape" begin
        MED = permutedims(MED)
        Idx = shuffle(axes(MED,2))
        imgshuff = MED[:,Idx]
        seqres = sequence(imgshuff; scales=(1,2), metrics=(KLD, L2))
        ind = order(seqres)
        res = imgshuff[:, ind]
        # This way the test will also fail if the 
        # reordering begins to work perfectly, i.e. loss ~ 0!
        @test 0.1 < loss(MED,res) < 1100
    end

    @testset "reorder med image, landscape, weightrows=true" begin
        MED = permutedims(MED)
        Idx = shuffle(axes(MED,2))
        imgshuff = MED[:,Idx]
        seqres = sequence(imgshuff; scales=(1,2), metrics=(KLD, L2), weightrows=true)
        ind = order(seqres)
        res = imgshuff[:, ind]
        # This way the test will also fail if the 
        # reordering begins to work perfectly, i.e. loss ~ 0!
        @test 0.1 < loss(MED,res) < 800
    end

    @testset "reorder large image, orig. landscape" begin
        Idx = shuffle(axes(BIG,2))
        imgshuff = BIG[:,Idx]
        seqres = sequence(imgshuff; scales=(1,2,4), metrics=(L2, KLD))
        ind = order(seqres)
        res = imgshuff[:, ind]
        @test 0.1 < loss(BIG,res) < 83
    end

    @testset "reorder large image, orig. landscape, weightrows=true" begin
        Idx = shuffle(axes(BIG,2))
        imgshuff = BIG[:,Idx]
        seqres = sequence(imgshuff; scales=(1,2,4), metrics=(L2, KLD), weightrows=true)
        ind = order(seqres)
        res = imgshuff[:, ind]
        # This way the test will also fail if the 
        # reordering begins to work perfectly, i.e. loss ~ 0!
        @test 0.1 < loss(BIG,res) < 11500
    end

    # So far, this is the only case that systematically produces ~ 0 loss
    @testset "reorder large image, portrait" begin
        BIG = permutedims(BIG)
        Idx = shuffle(axes(BIG,2))
        imgshuff = BIG[:,Idx]
        seqres = sequence(imgshuff; scales=(1,2,4), metrics=(L2, KLD))
        ind = order(seqres)
        res = imgshuff[:, ind]
        @test loss(BIG,res) ≈ 0
    end

    # Weighting rows performs WORSE in this case!
    @testset "reorder large image, portrait, weightrows=true" begin
        BIG = permutedims(BIG)
        Idx = shuffle(axes(BIG,2))
        imgshuff = BIG[:,Idx]
        seqres = sequence(imgshuff; scales=(1,2,4), metrics=(L2, KLD), weightrows=true)
        ind = order(seqres)
        res = imgshuff[:, ind]
        # This way the test will also fail if the 
        # reordering begins to work perfectly, i.e. loss ~ 0!
        @test 0.1 < loss(BIG,res) < 11500
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

    # TODO #19 figure out how to test for zero log output....
    @testset "silent mode does not fail" begin
        @test sequence(SMALL; scales=(2,), metrics=(L2,), silent=true) !== nothing
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
