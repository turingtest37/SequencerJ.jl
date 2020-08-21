using Test
using SequencerJulia
using StatsBase
using Random
using LightGraphs
using Images

@testset "SequencerJulia.jl" begin

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
        @show loss(SMALL,res)
        @test SMALL == res
    end

    @testset "reorder med image, orig. portrait" begin
        Idx = shuffle(axes(MED,2))
        imgshuff = MED[:,Idx]
        seqres = sequence(imgshuff; scales=(1,2), metrics=(KLD, L2))
        ind = order(seqres)
        res = imgshuff[:, ind]
        @show loss(MED,res)
        @test MED == res
    end

    @testset "reorder med image, landscape" begin
        MED = permutedims(MED)
        Idx = shuffle(axes(MED,2))
        imgshuff = MED[:,Idx]
        seqres = sequence(imgshuff; scales=(1,2), metrics=(KLD, L2))
        ind = order(seqres)
        res = imgshuff[:, ind]
        @show loss(MED,res)
        @test MED == res
    end

    @testset "reorder large image, orig. landscape" begin
        Idx = shuffle(axes(BIG,2))
        imgshuff = BIG[:,Idx]
        seqres = sequence(imgshuff; scales=(1,2,4), metrics=(L2, KLD))
        ind = order(seqres)
        res = imgshuff[:, ind]
        @show loss(BIG,res)
        @test BIG == res
    end

    @testset "reorder large image, portrait" begin
        BIG = permutedims(BIG)
        Idx = shuffle(axes(BIG,2))
        imgshuff = BIG[:,Idx]
        seqres = sequence(imgshuff; scales=(1,2,4), metrics=(L2, KLD))
        ind = order(seqres)
        res = imgshuff[:, ind]
        @show loss(BIG,res)
        @test BIG == res
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

end
