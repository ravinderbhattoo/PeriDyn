@testset "materials/bond_based.jl" begin
    @safetestset "block1 block2" begin
        using PeriDyn
        PD = PeriDyn
        env = PD.virgin_state2(1.1)

        blockid = 1
        mask = env.type.==env.material_blocks[blockid].type
        acc = PD.force_density_T(env.y[:,mask],env.material_blocks[blockid])/env.material_blocks[blockid].general.density
        ans1 = [PD._magnitude(acc[:,i]) for i in 1:8]
        ans2 = [0.1+0.1*sqrt(3)+0.1*sqrt(6) for i in 1:8]
        @test(isapprox(ans1,ans2))

        blockid = 2
        mask = env.type.==env.material_blocks[blockid].type
        acc = PD.force_density_T(env.y[:,mask],env.material_blocks[blockid])/env.material_blocks[blockid].general.density
        ans1 = [PD._magnitude(acc[:,i]) for i in 1:8]
        ans2 = [0.01+0.01*sqrt(3)+0.01*sqrt(6) for i in 1:8]
        @test(isapprox(ans1,ans2))
    end
end
