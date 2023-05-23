@testset "materials/bond_based.jl" begin
    @safetestset "block1 block2" begin
        using PeriDyn
        PD = PeriDyn

        env = PD.virgin_state2(1.1)

        blockid = 1
        mat = env.material_blocks[blockid]
        mask = env.type .== mat.type
        acc = PD.force_density_T(env.y[:,mask], mat) / first(mat.specific.density)
        print(acc)
        ans1 = [PD.get_magnitude(acc[:,i]) for i in 1:8]
        ans2 = [0.1+0.1*sqrt(3)+0.1*sqrt(6) for i in 1:8] / first(mat.specific.density) / first(mat.general.volume)

        @test(isapprox(ans1,ans2))

        blockid = 2
        mat = env.material_blocks[blockid]
        mask = env.type .== mat.type
        acc = PD.force_density_T(env.y[:,mask], mat) / first(mat.specific.density)
        print(acc)
        ans1 = [PD.get_magnitude(acc[:,i]) for i in 1:8]
        ans2 = [0.01+0.01*sqrt(3)+0.01*sqrt(6) for i in 1:8] / first(mat.specific.density) / first(mat.general.volume)

        @test(isapprox(ans1,ans2))
    end
end
