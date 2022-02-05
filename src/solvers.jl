export save_state!, update_acc!, velocity_verlet_step!, velocity_verlet!, minimize!, quasi_static!, update_neighs!

function filepath_(file_prefix; append_date=false)
    if append_date
        dtime = replace(string(ceil(Dates.now(), Dates.Second(1))), ":"=>"-")
        return mkpath("./output/"*file_prefix*"_"*dtime)*"/"
    else
        return mkpath("./output/"*file_prefix)*"/"
    end
end


"""
"""
function save_state!(filename, env)
    update_acc!(env)
    write_data(filename, 1:length(env.type), env.type, env.y, env.v, env.f, env.mass, env.volume)
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

̈u(xᵢ, t) = from material forces + from contact forces
 
from material force, ̈u(xᵢ, t) = [∑ᵏⱼ₌₁ {T[xᵢ, t]<xⱼ-xᵢ> - T[xⱼ, t]<xᵢ-xⱼ> }*Vⱼ + b(xᵢ, t)] / ρᵢ 

"""
function update_acc!(env::GeneralEnv)
    # fill force vector with zeros
    fill!(env.f, 0.0)
    
    ###################
    # MATERIAL FORCES #
    ###################
    # calculate force density * volume i.e
    # ̈ρu(xᵢ, t) = [∑ᵏⱼ₌₁ {T[xᵢ, t]<xⱼ-xᵢ> - T[xⱼ, t]<xᵢ-xⱼ> }*Vⱼ + b(xᵢ, t)]  
    # for each block (type)    
    for i in 1:size(env.material_blocks,1)
        mask = false
        for j in env.material_blocks[i].type
            m = (env.type .== j)
            mask = mask .| m
        end
        env.f[:,mask] .+= force_density(env.y[:,mask], env.material_blocks[i])
        
        # if NaN occurs, throw error
        @check_nan env.f "env.f from material forces ($(env.material_blocks[i].name))."
    end
    
    # if NaN occurs, throw error
    @check_nan env.f "env.f from material forces."
    
    ##################
    # CONTACT FORCES #
    ##################
    # calculate short range repulsion (should be force density * volume)
    for i in 1:size(env.short_range_repulsion,1)
        short_range_repulsion!(env.y, env.f, env.type, env.volume, env.short_range_repulsion[i])
    end
    
    # if NaN occurs, throw error
    @check_nan env.f "env.f from repulsive forces."
    
    # calculate aceeleration and momentum from force density and velocity respectively 
    env.p[1, :] = 1*env.volume
    env.p[2, :] = 1*env.volume
    env.p[3, :] = 1*env.volume
    for i in 1:size(env.material_blocks,1)
        for j in env.material_blocks[i].type
            m = (env.type .== j)
            t = j - env.material_blocks[i].type.start + 1
            
            # acceleration = (force density * volume) / density
            env.f[:, m] ./= env.material_blocks[i].specific.density[t]
            
            # momentum = velocity * volume * density
            env.p[:, m] .*= env.v[:, m] .* env.material_blocks[i].specific.density[t]
        end
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
    for bc in env.boundary_conditions
        check!(bc, env)
    end
    env.time_step += 1
end


"""
    velocity_verlet!(envs::Any,N::Int64;filewrite_freq=10,neigh_update_freq=50,file_prefix="datafile",start_at::Int64=0,write_from::Int=0)

Velocity verlet :).
"""
function velocity_verlet!(envs::Any, N::Int64; filewrite_freq=10,average_prop_freq=10,neigh_update_freq=50, out_dir="datafile",start_at::Int64=0,write_from::Int=0)
    mkpath(out_dir)
    
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
    print_data_file!(envs, out_dir, 0+write_from)
 
    
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
            print_data_file!(envs, out_dir, i+write_from)
        end
        
        # update neighbors
        if i%neigh_update_freq==0.0
            update_neighs!(envs)
        end
        println(round(i/N*100, digits=3),"%")
    end
end
SOLVERS[:vv] = velocity_verlet!


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
    for bc in env.boundary_conditions
        check!(bc, env)
    end
    env.time_step += 1
end


"""
    quasi_static!(envs::Any,N::Int64,step_size::Float64; max_iter::Int64=100, min_step_tol_per::Float64=0.5, filewrite_freq::Int64=10, neigh_update_freq::Int64=50, file_prefix::String="datafile",start_at::Int64=0,write_from::Int=0)

Implements quasi static simulation using minimize for each step.
"""
function quasi_static!(envs::Any,N::Int64,step_size::Float64; max_iter::Int64=100, filewrite_freq::Int64=10, neigh_update_freq::Int64=1, out_dir::String="datafile",start_at::Int64=0,write_from::Int=0)
    mkpath(out_dir)
    print_data_file!(envs, out_dir, 0)
    update_neighs!(envs)
    
    print_data_file!(envs, out_dir, 0+write_from)
 
    
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
        if i%filewrite_freq==0.0 || (i==1 && write_from==0)
            print_data_file!(envs, out_dir, i+write_from)
        end
        
        if i%neigh_update_freq==0.0
            update_neighs!(envs)
        end
        println(round(i/N*100, digits=3),"%")
    end
end
SOLVERS[:qs] = quasi_static!


"""
    function update_neighs!(envs)

Updates neighbors of each material point for a list of simulation environments.
"""
function update_neighs!(envs)
    println("\nUpdating neighbors for collision..................")
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
    println("\nWritting data file................................$i\n")
    for id in 1:size(envs,1)
        env = envs[id]
        filename = string(file_prefix,"/env_",env.id,"_step_",i,".data")
        save_state!(filename, env)
    end
    println("Done")
end
