function update_acc!(env)
    env.f .*= 0
    for i in 1:size(env.material_blocks,1)
        env.f[env.type.==env.material_blocks[i].type,:] .+= s_force_density_T(env.y[env.type.==env.material_blocks[i].type,:],env.material_blocks[i])/env.material_blocks[i].general.density
    end
    for i in 1:size(env.short_range_repulsion,1)
        short_range_repulsion!(env.y,env.f,env.type,env.short_range_repulsion[i])
    end
end


function velocity_verlet_step(env)
    update_acc!(env)
    c = env.dt/2
    env.v .+= c.*env.f
    env.y .+= env.dt.*env.v
    update_acc!(env)
    env.v .+= c.*env.f
end


function velocity_verlet(env::AbstractEnv;freq1=10,freq2=50)
    mkpath("./output")
    inc = 0
    for i in 1:env.N
        velocity_verlet_step(env)
        if i%freq1==0.0
            write_data(string("./output/datafile",inc),1,env.type,env.y)
            inc += 1
            per=i/env.N*100
            print(per,"% over\n")
        end
        if i%freq2==0.0
            for i in 1:size(env.short_range_repulsion,3)
            update_repulsive_neighs!(env.y,env.type,env.short_range_repulsion[i])
            end
        end
    end
    return pos
end
