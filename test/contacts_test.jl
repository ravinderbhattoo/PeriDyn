@testset "constacts/simple.jl" begin
    @safetestset "b1" begin
        using PeriDyn
        PD = PeriDyn
        env = PD.virgin_state(1.1)

        epsilon = 1000.0
        alpha = 2.0

        equi_dist = env.short_range_repulsion[1].equi_dist
        val = PD.repulsion_acc([0.75,0,0],1,env.short_range_repulsion[1])
        val_ = epsilon*(0.75-equi_dist)*(alpha-1)/env.short_range_repulsion[1].densities[1]
        @test isapprox(PD.magnitude(val),abs(val_))

        val = PD.repulsion_acc([0.75,0,0],2,env.short_range_repulsion[1])
        val_ = epsilon*(0.75-equi_dist)*(alpha-1)/env.short_range_repulsion[1].densities[2]
        @test isapprox(PD.magnitude(val),abs(val_))

        equi_dist = env.short_range_repulsion[2].material.particle_size
        val = PD.repulsion_acc([0,0.75,0],env.short_range_repulsion[2])
        val_ = epsilon*(0.75-equi_dist)*(alpha-1)/env.short_range_repulsion[2].material.density
        @test isapprox(PD.magnitude(val),abs(val_))

        equi_dist = env.short_range_repulsion[3].material.particle_size
        val = PD.repulsion_acc([0,0,0.75],env.short_range_repulsion[3])
        val_ = epsilon*(0.75-equi_dist)*(alpha-1)/env.short_range_repulsion[3].material.density
        @test isapprox(PD.magnitude(val),abs(val_))
    end
end


@testset "contacts/contacts.jl" begin
    @safetestset "b1" begin
        using PeriDyn
        PD = PeriDyn

        env = PD.virgin_state(0.9)

        epsilon = 1000.0
        alpha = 2.0

        val1 = -PD.magnitude(PD.repulsion_acc([0.6,0,0],1,env.short_range_repulsion[1]))
        val2 = PD.magnitude(PD.repulsion_acc([0.6,0,0],2,env.short_range_repulsion[1]))
        PD.update_repulsive_neighs!(env.y,env.type,env.short_range_repulsion[1])
        PD.short_range_repulsion!(env.y,env.f,env.type,env.short_range_repulsion[1])

        @test isapprox(env.f[1,:],[0,0,0,0,val1,val1,val1,val1,val2,val2,val2,val2,0,0,0,0])
    end

    @safetestset "b2" begin
        using PeriDyn
        PD = PeriDyn
        env = PD.virgin_state(0.9)

        env.material_blocks[1].general.intact .*= 0
        env.material_blocks[2].general.intact .*= 0

        epsilon = 1000.0
        alpha = 2.0

        val1 = PD.magnitude(PD.repulsion_acc([0.9,0,0],env.short_range_repulsion[2]))
        val2 = PD.magnitude(PD.repulsion_acc([0.9,0,0],env.short_range_repulsion[3]))
        PD.update_repulsive_neighs!(env.y,env.type,env.short_range_repulsion[3])
        PD.short_range_repulsion!(env.y,env.f,env.type,env.short_range_repulsion[3])
        PD.update_repulsive_neighs!(env.y,env.type,env.short_range_repulsion[2])
        PD.short_range_repulsion!(env.y,env.f,env.type,env.short_range_repulsion[2])

        val_ = zeros(16)
        val_[1:4] .= -val1
        val_[5:8] .= +val1
        val_[9:12] .= -val2
        val_[13:16] .= +val2

        @test isapprox(env.f[1,:],val_)
    end

end
