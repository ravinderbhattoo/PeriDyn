mutable struct GeneralEnv
    id::Int64
    type::Array{Int64,1}
    ghost_atoms::Array{Int64,1}
    state::Int64
    y::Array{Float64,2}
    v::Array{Float64,2}
    f::Array{Float64,2}
    time_step::Int64
    dt::Float64
    neighs::Array{Int64,2}
    boundary_conditions::Any
    short_range_repulsion::Any
    material_blocks::Any
end


function Env(id::Int64,materials,short_range_repulsion,boundary_conds,dt;state=2)
    type = materials[1].type*ones(Int64,size(materials[1].general.y,2))
    for i in 2:size(materials,1)
        type = vcat(type,materials[i].type*ones(Int64,size(materials[i].general.y,2)))
    end
    y = materials[1].general.y
    for mat in materials[2:end]
        y = hcat(y,mat.general.y)
    end
    v = materials[1].general.velocity
    for mat in materials[2:end]
        v = hcat(v,mat.general.velocity)
    end
    return GeneralEnv(id,type,0*type,state,y,v,0v,0,dt,zeros(2,2),boundary_conds,short_range_repulsion,materials)
end

function set_ghost_atoms!(env,ghost_atoms)
    env.ghost_atoms[1:end,1:end] = ghost_atoms
end

function set_env_active!(env)
    env.state = 2
end

function set_env_idel!(env)
    env.state = 1
end

function set_env_inactive!(env)
    env.state = 0
end

#
