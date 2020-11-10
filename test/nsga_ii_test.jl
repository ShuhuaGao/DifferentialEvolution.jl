@testset "NSGA-II" begin
    @testset "nondominated sorting" begin
        ws=[-1, -1]
        let fA = [20, 2.2], fB = [60, 4.4], fD = [15, 4.4]
            @test !dominate(fA, fD, ws)
            @test !dominate(fD, fA, ws)
            @test !dominate(fA, fA, ws)
            @test dominate(fA, fB, ws)
            @test dominate(fD, fB, ws)
            @test !dominate(fB, fD, ws)
        end
        fitness = [20 2.2; 60 4.4; 65 3.5; 15 4.4; 55 4.5; 50 1.8; 80 4.0; 25 4.6]'
        fronts = nondominated_sort(fitness, ws)
        # the order may differ, though we don't care about the order within each front
        @test all(issetequal.(fronts, [[1, 4, 6], [2, 3, 5, 8], [7]]))
    end
    
    @testset "crowding distance assignment" begin
        fitness = [20 2.2; 60 4.4; 65 3.5; 15 4.4; 55 4.5; 50 1.8; 80 4.0; 25 4.6]'
        ws=[-1, -1]
        fronts = nondominated_sort(fitness, ws)
        cds = assign_crowding_distance(fitness, fronts)
        @test cds[1] == 2.0
        @test cds[4] == Inf
        @test cds[6] == Inf
        @test cds[8] == Inf
        @test cds[3] == Inf
        @test cds[5] ≈ 1.0568181818181819
        @test cds[2] ≈ 1.1590909090909092
        @test cds[7] == Inf
    end
    @testset "crowded comparison" begin
        @test crowded_compare(1, 0.1, 4, Inf) == -1
        @test crowded_compare(5, 0.1, 4, Inf) == 1
        @test crowded_compare(1, 0.1, 1, Inf) == 1
        @test crowded_compare(1, Inf, 1, Inf) == 0
    end
end