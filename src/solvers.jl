function update_acc!(env::GeneralEnv)
    fill!(env.f,0.0)
    for i in 1:size(env.material_blocks,1)
        mask = env.type.==env.material_blocks[i].type
        env.f[:,mask] .+= s_force_density_T(env.y[:,mask],env.material_blocks[i])/env.material_blocks[i].general.density
    end
    for i in 1:size(env.short_range_repulsion,1)
        short_range_repulsion!(env.y,env.f,env.type,env.short_range_repulsion[i])
    end
end


function velocity_verlet_step!(env::GeneralEnv)
    c = env.dt/2
    env.v .+= c.*env.f
    env.y .+= env.dt.*env.v
    for bc in env.boundary_conditions
        apply_bc!(env,bc)
    end
    update_acc!(env)
    env.v .+= c.*env.f
    env.time_step += 1
end


function velocity_verlet!(envs::Any,N::Int64;freq1=10,freq2=50,file_prefix="datafile",start_at::Int64=0)
    mkpath("./output")
    print("\nUpdating neighbors for collision..................")
    for id in 1:size(envs,1)
        env = envs[id]
        for rm in 1:size(env.short_range_repulsion,1)
            update_repulsive_neighs!(env.y,env.type,env.short_range_repulsion[rm])
        end
    end
    print("Done\n")
    N = N + start_at
    for i in (1+start_at):N
        for id in 1:size(envs,1)
            if envs[id].state==2
                velocity_verlet_step!(envs[id])
            end
        end

        if i%freq1==0.0 || i==1
            print("\nWritting data file.......................")
            for id in 1:size(envs,1)
                env = envs[id]
                write_data(string("./output/",file_prefix,"_env_",env.id,"_step_",i,".data"),env.type,env.y,env.v,env.f)
            end
            print("Done\n")
            print(i/N*100,"% over\n")
        end

        if i%freq2==0.0
            print("\nUpdating neighbors for collision..................")
            for env in envs
                if env.state[1]==2
                    for rm in 1:size(env.short_range_repulsion,1)
                        update_repulsive_neighs!(env.y,env.type,env.short_range_repulsion[rm])
                    end
                end
            end
            print("Done\n")
        end
    end
end
