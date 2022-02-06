@testset "materials/ordinary_state_based.jl" begin
    @safetestset "block1 block2" begin
        using PeriDyn
        PD = PeriDyn
        env = PD.virgin_state(1.1)
        
        blockid = 1
        type_ = env.material_blocks[blockid].type
        mask = env.type .== type_
        Gen = env.material_blocks[blockid].general
        Spe = env.material_blocks[blockid].specific
        theta = PD.dilatation(env.y[:,mask],Gen,Gen.weighted_volume)
        m = Gen.weighted_volume
        K = first(Spe.bulk_modulus)
        G = first(Spe.shear_modulus)
        s = [0,0,0]
        
        for i in 2:8
            j = 1
            X = Gen.x[:,i] - Gen.x[:,j]
            Y = Gen.y[:,i] - Gen.y[:,j]
            xij = PD.get_magnitude(X)
            wij = PD.influence_function(X)
            wji = PD.influence_function(-X)
            e = (PD.get_magnitude(Y) - PD.get_magnitude(X))
            t1 = (3*K-5*G) * ( theta[i]*xij*wij/m[i] + theta[j]*xij*wji/m[j] )
            t2 = 15*G * ( e*wij/m[i] + e*wji/m[j] )  
            t = ( t1 + t2 ) * Gen.volume[j] / PD.get_magnitude(Y)
            s += t * Y
        end
        
        t = PD.force_density_T(env.y[:,mask],env.material_blocks[blockid])
        @test(isapprox(s, t[:,1]))
    end
end
