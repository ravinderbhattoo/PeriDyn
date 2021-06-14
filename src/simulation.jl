export Env, @env

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
    type = materials[1].general.type
    for i in 2:size(materials,1)
        type = vcat(type, materials[i].general.type)
    end
    y = materials[1].general.y
    for mat in materials[2:end]
        y = hcat(y,mat.general.y)
    end
    v = materials[1].general.velocity
    for mat in materials[2:end]
        v = hcat(v,mat.general.velocity)
    end
    volume = materials[1].general.volume
    for mat in materials[2:end]
        volume = vcat(volume,mat.general.volume)
    end
    return GeneralEnv(id,type,0*type,state,y,v,0v,0v,volume,0,dt,zeros(2,2),boundary_conds,short_range_repulsion,materials,nothing,nothing,nothing)
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

#
