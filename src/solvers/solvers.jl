export save_state!, update_acc!, update_neighs!

function filepath_(file_prefix; append_date=false)
    if append_date
        dtime=replace(string(ceil(Dates.now(), Dates.Second(1))), ":"=>"-")
        return mkpath("./output/"*file_prefix*"_"*dtime)*"/"
    else
        return mkpath("./output/"*file_prefix)*"/"
    end
end


"""
    save_state!(filename, env)

Save `env::GeneralEnv` to disk.
"""
function save_state!(filename, env)
    update_acc!(env)
    write_data(filename; id=1:length(env.type), type=env.type, position=env.y, velocity=env.v, acceleration=env.f, mass=env.mass, volume=env.volume)
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
            
            # acceleration = (force_density[force/vol/vol] * volume) / density
            env.f[:, m] ./= env.material_blocks[i].specific.density[t]
            
            # momentum = velocity * volume * density
            env.p[:, m] .*= env.v[:, m] .* env.material_blocks[i].specific.density[t]
        end
    end

end


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
function print_data_file!(envs::Array{GeneralEnv}, file_prefix::String, i::Int64; ext::Symbol=:jld)
    println("\nWritting data file................................$i\n")
    for id in 1:size(envs,1)
        env = envs[id]
        filename = string(file_prefix,"/env_",env.id,"_step_",i,".$(String(ext))")
        save_state!(filename, env)
    end
    println("Done")
end

include("./velocity_verlet.jl")
include("./minimize.jl")
