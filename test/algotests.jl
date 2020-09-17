

@testset "Scaling and Distance Matrices" begin

    #### L2 = Euclidean
        @testset "Euclidean Distance Matrix" begin
            objs = [[1,1,1], [2,3,4]]
            expected = [0. sqrt(14) ; sqrt(14) 0.]
            A = permutedims(hcat(objs...)) # makes 2x3 array
            @test pairwise(Euclidean(), A, dims = 1) ≈ expected
        end

        @testset "Euclidean simple" begin
            objs = [[1,1,1], [2,3,4]]
            expected = float(sqrt(14))
            @test Euclidean()(objs...) ≈ expected
        end


    #### KL Divergence
        @testset "Same as scipy.stats.entropy" begin
            u,v = [1/2, 1/2], [9/10, 1/10]
            expected = 0.5108256237659907
            @test KLDivergence()(u,v) ≈ expected
        end

    ##### EMD or 1-d Wasserstein

        # EMD or 1-d Wasserstein
        @testset "EMD No Fail" begin
            objs = [[1,1,1], [2,3,4]]
            @test EMD(objs[1], objs[2]) !== nothing
        end


        @testset "EMD with grid provided, same result as python" begin
            objs = [[3.4, 3.9, 7.5, 7.8], [4.5, 1.4],
            [1.4, 0.9, 3.1, 7.2], [3.2, 3.5]]
            expected = 4.0781331438047861
            @test EMD(objs...) ≈ expected
        end

        @testset "EMD grid provided, different constructor, no fail" begin
            g = [0.5, 1.0, 1.5, 2.0]
            u = [3.4, 3.9, 7.5, 7.8]
            v = [1.4, 0.9, 3.1, 7.2]
            @test EMD((g,g))(u,v) !== nothing
        end

        @testset "EMD with grid, different constructor, same result as python" begin
            g1 = [3.4, 3.9, 7.5, 7.8]
            g2 = [4.5, 1.4]
            u = [1.4, 0.9, 3.1, 7.2]
            v = [3.2, 3.5]
            expected = 4.0781331438047861
            @test EMD((g1,g2))(u,v) ≈ expected
        end


    ##### Energy

        @testset "Energy() no fail" begin
            objs = [[1,1,1], [2,3,4]]
            @test Energy(objs[1], objs[2]) !== nothing
        end

        @testset "energy() not null" begin
            objs = [[3.4, 3.9, 7.5, 7.8], [4.5, 1.4],
            [1.4, 0.9, 3.1, 7.2], [3.2, 3.5]]
            expected = 4.0781331438047861
            @test energy(objs...) !== nothing
        end

        @testset "Energy() produces same result as python, simple." begin
            objs = ([0, 8], [0, 8], [3, 1], [2, 2])
            expected = 1.0
            @test Energy()(objs...) ≈ expected
        end

        @testset "Energy() produces same result as python, harder." begin
            objs = [[3.4, 3.9, 7.5, 7.8], [4.5, 1.4],
            [1.4, 0.9, 3.1, 7.2], [3.2, 3.5]]
            expected = 2.3674181638546394
            @test Energy()(objs...) ≈ expected
        end


        @testset "energy() w/ different u and v. Compare to python result." begin
            objs = ([0.7, 7.4, 2.4, 6.8], [1.4, 8. ],[2.1, 4.2, 7.4, 8. ], [7.6, 8.8])
            expected = Float32(0.8800334097615822) # from scipy.stats.energy_distance
            @test energy(objs...) ≈ expected
        end

        @testset "energy() w/ different u and v. Compare to python result. Float32 math." begin
            objs = (Float32[0.7, 7.4, 2.4, 6.8], Float32[1.4, 8. ],Float32[2.1, 4.2, 7.4, 8. ], Float32[7.6, 8.8])
            expected = 0.8800334097615822 # from scipy.stats.energy_distance
            @test Float32(energy(objs...)) ≈ expected
        end

        @testset "energy((g1,g2)) w/ different u and v. " begin
            g1 = [0.7, 7.4, 2.4, 6.8]
            g2 = [1.4, 8. ]
            u = [2.1, 4.2, 7.4, 8. ]
            v = [7.6, 8.8]
            expected = 0.8800334097615822 # from scipy.stats.energy_distance
            @test Energy((g1,g2))(u,v) ≈ expected
        end


end
