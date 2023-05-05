export Env, @env, simulate

mutable struct GeneralEnv
    id::Int64
    type::AbstractArray{Int64,1}
    bid::AbstractArray{Int64,1}
    ghost_atoms::AbstractArray{Int64,1}
    state::Int64
    y::AbstractArray{Float64,2}
    v::AbstractArray{Float64,2}
    f::AbstractArray{Float64,2}
    p::AbstractArray{Float64,2}
    volume::AbstractArray{Float64,1}
    intact0::AbstractArray{Int64, 1}
    mass::AbstractArray{Float64,1}
    time_step::Int64
    dt::Float64
    neighs::AbstractArray{Int64,2}
    boundary_conditions::Any
    short_range_repulsion::Any
    material_blocks::Any
    boundaries::Tuple
    Collect!::Any
    Params::Any
    Out::Any
end

function _cudaconvert(x::Vector{GeneralEnv})
    _cudaconvert.(x)
end

function _cudaconvert(x::T) where T <: GeneralEnv
    function fn(x, k)
        if k!=:Out
            _cudaconvert(getfield(x, k))
        else
            getfield(x, k)
        end
    end
    T(( _cudaconvert(getfield(x, k)) for k in fieldnames(T))...)
end

"""
    Env(id::Int64,materials,short_range_repulsion,boundary_conds,dt;state=2)

Create a GeneralEnv for holding parameters for a simulation.
"""
function Env(id::Int64,materials,short_range_repulsion,boundary_conds,dt;state=2,bskin=0.5)
    mat = materials[1]
    type = mat.general.type
    intact = sum(mat.general.intact, dims=1)
    y = mat.general.y
    v = mat.general.velocity
    volume = mat.general.volume
    mass = 1*volume
    bid = fill(mat.blockid, length(mass))

    for j in mat.type
        mask = (mat.general.type .== j)
        t = j - mat.type.start + 1
        mass[mask] .*= mat.specific.density[t]
    end

    for mat in materials[2:end]
        type = vcat(type, mat.general.type)
        y = hcat(y, mat.general.y)
        v = hcat(v, mat.general.velocity)
        volume = vcat(volume,mat.general.volume)
        intact = hcat(intact, sum(mat.general.intact, dims=1))
        mass_ = 1*mat.general.volume
        for j in mat.type
            mask = (mat.general.type .== j)
            t = j - mat.type.start + 1
            mass_[mask] .*= mat.specific.density[t]
        end

        mass = vcat(mass,mass_)
        bid_ = fill(mat.blockid, length(mass_))
        bid = vcat(bid,bid_)

    end

    _min = minimum(y, dims=2)
    _max = maximum(y, dims=2)
    _cm = @. (_min + _max ) / 2
    _L = @. _max - _min

    boundaries = (_cm .- (0.5 + bskin)*_L, _cm .+ (0.5 + bskin)*_L)


    env = GeneralEnv(id,type,bid,0*type,state,y,v,0v,0v,volume, reshape(intact, :),
                mass,0,dt,zeros(2,2),
                boundary_conds,short_range_repulsion,materials,boundaries,
                nothing,nothing,nothing)

    return deviceconvert(env)
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
