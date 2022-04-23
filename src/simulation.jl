export Env, @env, simulate

mutable struct GeneralEnv
    id::Int64
    type::Array{Int64,1}
    ghost_atoms::Array{Int64,1}
    state::Int64
    y::Array{Float64,2}
    v::Array{Float64,2}
    f::Array{Float64,2}
    p::Array{Float64,2}
    volume::Array{Float64,1}
    mass::Array{Float64,1}
    time_step::Int64
    dt::Float64
    neighs::Array{Int64,2}
    boundary_conditions::Any
    short_range_repulsion::Any
    material_blocks::Any
    Collect!::Any
    Params::Any
    Out::Any
end

"""
    Env(id::Int64,materials,short_range_repulsion,boundary_conds,dt;state=2)

Create a GeneralEnv for holding parameters for a simulation.
"""
function Env(id::Int64,materials,short_range_repulsion,boundary_conds,dt;state=2)
    mat = materials[1]
    type = mat.general.type
    y = mat.general.y
    v = mat.general.velocity
    volume = mat.general.volume
    mass = 1*volume
    for j in mat.type
        mask = (mat.general.type .== j)
        t = j - mat.type.start + 1
        mass[mask] .*= mat.specific.density[t]
    end

    for mat in materials[2:end]
        type = vcat(type, mat.general.type)
        y = hcat(y,mat.general.y)
        v = hcat(v,mat.general.velocity)
        volume = vcat(volume,mat.general.volume)

        mass_ = 1*mat.general.volume
        for j in mat.type
            mask = (mat.general.type .== j)
            t = j - mat.type.start + 1
            mass_[mask] .*= mat.specific.density[t]
        end

        mass = vcat(mass,mass_)
    end
    return GeneralEnv(id,type,0*type,state,y,v,0v,0v,volume,mass,0,dt,zeros(2,2),boundary_conds,short_range_repulsion,materials,nothing,nothing,nothing)
end

"""
    set_ghost_atoms!(env,ghost_atoms)

Set ghost atoms for an environment.
"""
function set_ghost_atoms!(env,ghost_atoms)
    env.ghost_atoms[1:end,1:end] = ghost_atoms
end

"""
    set_env_active!(env)

Set environment state active.
"""
function set_env_active!(env)
    env.state = 2
end

"""
    set_env_idel!(env)

Set environment state idel.
"""
function set_env_idel!(env)
    env.state = 1
end


"""
    set_env_inactive!(env)

Set environment state inactive.
"""
function set_env_inactive!(env)
    env.state = 0
end


"""
"""
function simulate(args...; solver=:vv, out_dir="datafile", append_date=true, kwargs...)
    sim = SOLVERS[solver]
    println("Using solver: $(sim)")
    foldername = filepath_(out_dir; append_date=append_date)
    println("Output folder: $foldername")
    return sim(args...; kwargs..., out_dir=foldername)
end



#
