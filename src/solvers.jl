
export save_state!, update_acc!, velocity_verlet_step!, velocity_verlet!, minimize!, quasi_static!, update_neighs!
using Flux
using Dates

function filepath_(file_prefix)
    dtime = replace(string(ceil(Dates.now(), Dates.Second(1))), ":"=>"-")
    return mkpath("./output/"*file_prefix*"_"*dtime)*"/"
end


"""
"""
function save_state!(filename, env)
    update_acc!(env)
    write_data(filename, env.type, env.y, env.v, env.f)
    print("Initial state saved.")
end

function apply_bc_at0(env, start_at)
    if start_at==0
        for bc in env.boundary_conditions
            apply_bc_at0!(env, bc)
        end
    end
end

"""
    update_acc!(env::GeneralEnv)

Updates acceleration of all the material points in a simulation environment.
"""
function update_acc!(env::GeneralEnv)
    fill!(env.f, 0.0)
    for i in 1:size(env.material_blocks,1)
        mask = false
        for j in env.material_blocks[i].type
            m = (env.type .== j)
            mask = mask .| m
        end
        env.f[:,mask] .+= force_density_T(env.y[:,mask], env.material_blocks[i]) 
    end
    for i in 1:size(env.short_range_repulsion,1)
        short_range_repulsion!(env.y,env.f,env.type,env.short_range_repulsion[i])
    end

    env.p[1, :] = 1*env.volume
    env.p[2, :] = 1*env.volume
    env.p[3, :] = 1*env.volume
    for i in 1:size(env.material_blocks,1)
        for j in env.material_blocks[i].type
            m = (env.type .== j)
            t = j - env.material_blocks[i].type.start + 1
            env.f[:, m] ./= env.material_blocks[i].specific.density[t]
            env.p[:, m] .*= env.v[:, m] .* env.material_blocks[i].specific.density[t]
        end
    end

    if any(isnan,env.f)
        error("NaN in env.f")
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
        apply_bc!(env, bc)
    end
    update_acc!(env)
    env.v .+= c.*env.f
    for bc in env.boundary_conditions
        apply_bc!(env, bc)
    end
    env.time_step += 1
end

"""
    velocity_verlet!(envs::Any,N::Int64;filewrite_freq=10,neigh_update_freq=50,file_prefix="datafile",start_at::Int64=0)

Velocity verlet :).
"""
function velocity_verlet!(envs::Any, N::Int64; filewrite_freq=10,average_prop_freq=10,neigh_update_freq=50,file_prefix="datafile",start_at::Int64=0)
    foldername = filepath_(file_prefix)
    print("\nUpdating neighbors for collision..................")
    for id in 1:size(envs,1)
        env = envs[id]
        for rm in 1:size(env.short_range_repulsion,1)
            update_repulsive_neighs!(env.y, env.type, env.short_range_repulsion[rm])
        end
    end
    print("Done\n")
    N = N + start_at
    for id in 1:size(envs,1)
        env = envs[id]
        filename = string(foldername,"env_",env.id,"_step_",0,".data")
        save_state!(filename, env)
    end
    for id in 1:size(envs,1)
        if envs[id].state==2
            apply_bc_at0(envs[id], start_at)
        end
    end
    for i in (1+start_at):N
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

        if i%average_prop_freq==0.0 || i==1
            for env in envs
                println("Momentum: $i", sum(env.p, dims=2))
            end
        end

        if i%filewrite_freq==0.0 || i==1
            print_data_file!(envs, foldername, i)
        end

        if i%neigh_update_freq==0.0
            update_neighs!(envs)
        end
        println(round(i/N*100, digits=3),"%")
    end
end


function minimize!(env::GeneralEnv, step_size::Float64; max_iter::Int64=50, x_tol::Float64=1.0e-6, f_tol::Float64=1.0e-6)
    for bc in env.boundary_conditions
        apply_bc!(env, bc)
    end
    update_acc!(env)

    opt = Flux.ADAM(step_size)
    mask = false
    for bc in env.boundary_conditions
        if ~(bc.onlyatstart)
            mask = mask .| bc.bool
        end
    end
    mask = .!(mask)
    println("Minimizing...")
    for i in 1:max_iter
        f_tol_ = maximum(abs.(env.f[:, mask]))
        Δy = Flux.Optimise.apply!(opt, env.y, env.f)
        env.y .+= Δy
        x_tol_ = maximum(abs.(Δy[:, mask]))

        for bc in env.boundary_conditions
            if ~(bc.onlyatstart)
                apply_bc!(env, bc)
            end
        end
        update_acc!(env)

        if ( (x_tol_ < x_tol) || (f_tol_ < f_tol) )
            println("Tolerance reached. $x_tol_ <= $x_tol or $f_tol_ <= $f_tol")
            break
        end
        if i==max_iter
            println("Maximum iteration ($max_iter) reached. $x_tol_ !<= $x_tol and $f_tol_ !<= $f_tol")
        end
    end
    env.time_step += 1
end


"""
    quasi_static!(envs::Any,N::Int64,step_size::Float64; max_iter::Int64=100, min_step_tol_per::Float64=0.5, filewrite_freq::Int64=10, neigh_update_freq::Int64=50, file_prefix::String="datafile",start_at::Int64=0)

Implements quasi static simulation using minimize for each step.
"""
function quasi_static!(envs::Any,N::Int64,step_size::Float64; max_iter::Int64=100, filewrite_freq::Int64=10, neigh_update_freq::Int64=1, file_prefix::String="datafile",start_at::Int64=0)
    foldername = filepath_(file_prefix)
    print_data_file!(envs, foldername, 0)
    update_neighs!(envs)
    for id in 1:size(envs,1)
        env = envs[id]
        filename = string(foldername,"env_",env.id,"_step_",0,".data")
        save_state!(filename,env)
    end

    for id in 1:size(envs,1)
        if envs[id].state==2
            apply_bc_at0(envs[id], start_at)
        end
    end
    N = N + start_at
    for i in (1+start_at):N
        for id in 1:size(envs,1)
            if envs[id].state==2
                fill!(envs[id].v,0.0)
                minimize!(envs[id],step_size; max_iter=max_iter)
            end
            if envs[id].Collect! !== nothing
                envs[id].Collect!(envs[id],i)
            end
        end
        if i%filewrite_freq==0.0 || i==1
            print_data_file!(envs, foldername, i)
        end

        if i%neigh_update_freq==0.0
            update_neighs!(envs)
        end
        println(round(i/N*100, digits=3),"%")
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
end

"""
    print_data_file!(envs::Array{GeneralEnv}, file_prefix::String, i::Int64)

Writes data file to disk.
"""
function print_data_file!(envs::Array{GeneralEnv}, file_prefix::String, i::Int64)
    print("\nWritting data file................................")
    for id in 1:size(envs,1)
        env = envs[id]
        write_data(string(file_prefix,"env_",env.id,"_step_",i,".data"),env.type,env.y,env.v,env.f)
    end
    print("Done\n")
end
