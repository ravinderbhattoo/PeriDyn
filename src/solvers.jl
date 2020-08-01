
export update_acc!, velocity_verlet_step!, velocity_verlet!, minimize!, quasi_static!, update_neighs!

"""
    update_acc!(env::GeneralEnv)

Updates acceleration of all the material points in a simulation environment.
"""
function update_acc!(env::GeneralEnv)
    fill!(env.f,0.0)
    for i in 1:size(env.material_blocks,1)
        mask = env.type.==env.material_blocks[i].type
        env.f[:,mask] .+= force_density_T(env.y[:,mask],env.material_blocks[i])/env.material_blocks[i].general.density
    end
    for i in 1:size(env.short_range_repulsion,1)
        short_range_repulsion!(env.y,env.f,env.type,env.short_range_repulsion[i])
    end
end

"""
    velocity_verlet_step!(env::GeneralEnv)

Implement a single step of velocity verlet algorithm.
"""
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

"""
    velocity_verlet!(envs::Any,N::Int64;freq1=10,freq2=50,file_prefix="datafile",start_at::Int64=0)

Velocity verlet :).
"""
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

"""
    minimize!(env::GeneralEnv,step_size::Float64; max_iter::Int64=500,min_step_tol_per::Float64=5.0)

Minimize potential energy of simulation environment.
"""
function minimize!(env::GeneralEnv,step_size::Float64; max_iter::Int64=500,min_step_tol_per::Float64=5.0)
    ref_particle_size = env.material_blocks[1].general.particle_size
    orig_step_size = step_size
    for i in 1:max_iter
        update_acc!(env)
        if maximum(abs.(env.f))*step_size < min_step_tol_per/100*ref_particle_size
            step_size *= 2
        elseif maximum(abs.(env.f))*step_size > ref_particle_size/2
            step_size /= 2
        else
        end
        env.y .+= step_size*env.f
        for bc in env.boundary_conditions
            apply_bc!(env,bc)
        end
        if step_size>1000*orig_step_size
            break
        end
    end
    env.time_step += 1
end

"""
    quasi_static!(envs::Any,N::Int64,step_size::Float64; max_iter::Int64=100, min_step_tol_per::Float64=0.5, freq1::Int64=10, freq2::Int64=50, file_prefix::String="datafile",start_at::Int64=0)

Implements quasi static simulation using minimize for each step.
"""
function quasi_static!(envs::Any,N::Int64,step_size::Float64; max_iter::Int64=100, min_step_tol_per::Float64=0.5, freq1::Int64=10, freq2::Int64=50, file_prefix::String="datafile",start_at::Int64=0)
    mkpath("./output")
    update_neighs!(envs)
    N = N + start_at
    for i in (1+start_at):N
        for id in 1:size(envs,1)
            if envs[id].state==2
                fill!(envs[id].v,0.0)
                minimize!(envs[id],step_size)
            end
            if envs[id].Collect! != nothing
                envs[id].Collect!(envs[id].Params,envs[id].Out,i)
            end
        end


        if i%freq1==0.0 || i==1
            print_data_file!(envs,file_prefix,i)
            print(i/N*100,"% over\n")
        end

        if i%freq2==0.0
            update_neighs!(envs)
        end
    end
end


"""
    function update_neighs!(envs)

Updates neighbors of each material point for a list of simulation environments.
"""
function update_neighs!(envs)
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

"""
    print_data_file!(envs::Array{GeneralEnv}, file_prefix::String,i::Int64)

Writes data file to disk.
"""
function print_data_file!(envs::Array{GeneralEnv}, file_prefix::String,i::Int64)
    print("\nWritting data file................................")
    for id in 1:size(envs,1)
        env = envs[id]
        write_data(string("./output/",file_prefix,"_env_",env.id,"_step_",i,".data"),env.type,env.y,env.v,env.f)
    end
    print("Done\n")
end
