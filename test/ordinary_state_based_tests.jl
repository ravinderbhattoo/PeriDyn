@testset "materials/ordinary_state_based.jl" begin
    @safetestset "block1 block2" begin
        using PeriDyn
        PD = PeriDyn
        env = PD.virgin_state(1.1)

        blockid = 1
        mask = env.type.==env.material_blocks[blockid].type
        Gen = env.material_blocks[blockid].general
        Spe = env.material_blocks[blockid].specific
        theta = PD.dilatation(env.y[:,mask],Gen,Spe.weighted_volume)
        m = Spe.weighted_volume
        K = Spe.bulk_modulus
        G = Spe.shear_modulus
        s = [0,0,0]

        for i in 2:8
            j = 1
            X = Gen.x[:,i] - Gen.x[:,j]
            Y = Gen.y[:,i] - Gen.y[:,j]
            xij = PD.magnitude(X)
            wij = PD.influence_function(X)
            wji = PD.influence_function(-X)
            e = (PD.magnitude(Y) - PD.magnitude(X))
            t = ( (((3*K-5*G)*(theta[i]*xij*wij/m[i]+theta[j]*xij*wji/m[j]) + 15*G*(e*wji/m[i]+e*wji/m[j]))) )*Gen.volume[j]/PD.magnitude(Y)
            s += t*Y
        end

        t = PD.force_density_T(env.y[:,mask],env.material_blocks[blockid])
        @test(isapprox(s,t[:,1]))
    end
end
