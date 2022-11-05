export velocity_verlet_step!, velocity_verlet!

"""
velocity_verlet_step!(env::GeneralEnv)

Implement a single step of velocity verlet algorithm.
"""
function velocity_verlet_step!(env::GeneralEnv)
    c = env.dt/2
    env.v .+= c.*env.f
    env.y .+= env.dt.*env.v
    
    for bc in env.boundary_conditions
        apply_bc!(env, bc, :position)
    end
    
    update_acc!(env)
    env.v .+= c.*env.f
    
    for bc in env.boundary_conditions
        apply_bc!(env, bc, :velocity)
    end
    
    for bc in env.boundary_conditions
        check!(bc, env)
    end

    env.time_step += 1
end


"""
velocity_verlet!(envs::Any,N::Int64;filewrite_freq=10,neigh_update_freq=50,file_prefix="datafile",start_at::Int64=0,write_from::Int=0)

Velocity verlet :).
"""
function velocity_verlet!(envs::Any, N::Int64; filewrite_freq=10,average_prop_freq=10,neigh_update_freq=50, out_dir="datafile",start_at::Int64=0,write_from::Int=0, ext::Symbol=:jld)
    mkpath(out_dir)
    _print_data_file!(step) = print_data_file!(envs, out_dir, step; ext=ext)

    # Updating neighbors for contact
    println("\nUpdating neighbors for collision..................")
    for id in 1:size(envs,1)
        env = envs[id]
        for rm in 1:size(env.short_range_repulsion,1)
            update_repulsive_neighs!(env.y, env.type, env.short_range_repulsion[rm])
        end
    end
    println("Done")

    # Starting simulation at start point 
    N = N + start_at

    # save initial state
    _print_data_file!(0+write_from)


    # apply boundary conditions at t = t0
    for id in 1:size(envs,1)
        if envs[id].state==2
            apply_bc_at0(envs[id], start_at)
        end
    end


    for i in (1+start_at):N
        
        # apply integrator step for all envs
        for id in 1:size(envs,1)
            if envs[id].state==2
                velocity_verlet_step!(envs[id])
                if envs[id].Collect! !== nothing
                    envs[id].Collect!(envs[id],i)
                end
            else
                println("Env $id is not active. step = $i")
            end
        end
        
        # calculate average property
        if i%average_prop_freq==0.0 || (i==1 && write_from==0)
            for env in envs
                println("Momentum: $i", sum(env.p, dims=2))
            end
        end
        
        # write file to disk
        if i%filewrite_freq==0.0 || i==1
            _print_data_file!(i+write_from)
        end
        
        # update neighbors
        if i%neigh_update_freq==0.0
            update_neighs!(envs)
        end
        println(round(i/N*100, digits=3),"%")
    end
end

SOLVERS[:vv] = velocity_verlet!

