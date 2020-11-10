@testset "Constraint" begin
    @testset "bound constraints" begin
        let bounds = [-10.0 10.0; 2.4 8.2; 0  3.4], r = [-11, 10.0, 1.23]
            bounce_back!(r, [-5.3, 0, 0], bounds)
            @test -10 <= r[1] <= -5.3 && 0 <= r[2] <= 8.2 && r[3] == 1.23
        end
    end
end