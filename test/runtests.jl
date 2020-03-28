using PeriDyn, Test

@testset "mysum" begin
    x, y  = 5, 7
    @test mysum(x, y) == 12
end
