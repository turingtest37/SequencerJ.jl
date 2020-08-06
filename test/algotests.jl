

@testset "Scaling and Distance Matrices" begin

    #### L2 = Squared Euclidean
        @testset "Squared Euclidean Distance Matrix" begin
            objs = [[1,1,1], [2,3,4]]
            expected = [0. 14. ; 14. 0.]
            A = hcat(objs...)' # makes 2x3 array
            @test pairwise(SqEuclidean(), A, dims = 1) == expected
        end

        @testset "SquaredEuclidean simple" begin
            objs = [[1,1,1], [2,3,4]]
            expected = float(14)
            @test SqEuclidean()(objs...) == expected
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


            @testset "EMD with grid" begin
                objs = [[3.4, 3.9, 7.5, 7.8], [4.5, 1.4],
                [1.4, 0.9, 3.1, 7.2], [3.2, 3.5]]
                expected = 4.0781331438047861
                @test EMD(objs...) ≈ expected
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


            @testset "Energy() produces same result as python" begin
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



end
