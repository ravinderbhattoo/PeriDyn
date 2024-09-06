export GeneralEnvironment, Env, printdata

"""
    SimulationEnvironment

Abstract type for holding parameters for a simulation.
"""
abstract type SimulationEnvironment end


"""
    GeneralEnvironment <: SimulationEnvironment

Type for holding parameters for a simulation.

# Fields
- `id::Int64`: ID of the environment
- `type::AbstractArray{Int64,1}`: Type of each particle
- `bid::AbstractArray{Int64,1}`: Block ID of each particle
- `ghost_particles::AbstractArray{Int64,1}`: Ghost particles
- `state::Int64`: State of the environment
- `time_step::Int64`: Time step of the environment
- `dt::T where T<:QF`: Time step size
- `y::AbstractArray{T,2} where T<:QF`: Position of each particle
- `v::AbstractArray{T,2} where T<:QF`: Velocity of each particle
- `f::AbstractArray{T,2} where T<:QF`: Force of each particle
- `p::AbstractArray{T,2} where T<:QF`: Momentum of each particle
- `volume::AbstractArray{T,1} where T<:QF`: Volume of each particle
- `mass::AbstractArray{T,1} where T<:QF`: Mass of each particle
- `density::AbstractArray{T,1} where T<:QF`: Density of each particle
- `intact0::AbstractArray{Int64, 1}`: Intact particles information
- `neighs::AbstractArray{Int64,2}`: Neighbors of each particle
- `boundary_conditions::AbstractArray{T, 1} where T`: Boundary conditions
- `short_range_contact::AbstractArray{T, 1} where T`: Short range contact
- `material_blocks::AbstractArray{T, 1} where T`: Material blocks
- `boundaries::Tuple`: Boundaries
- `Collect!::Function`: Function for collecting data
- `Params::Dict{Symbol, Any}`: Parameters
- `Out::Dict{Symbol, Any}`: Output
- `cprint::Function`: Function for printing
"""
mutable struct GeneralEnvironment <: SimulationEnvironment
    id::Int64
    type::AbstractArray{Int64,1}
    bid::AbstractArray{Int64,1}
    ghost_particles::AbstractArray{Int64,1}
    state::Int64

    time_step::Int64
    dt::T where T<:QF

    y::AbstractArray{T,2} where T<:QF
    v::AbstractArray{T,2} where T<:QF
    f::AbstractArray{T,2} where T<:QF
    p::AbstractArray{T,2} where T<:QF
    volume::AbstractArray{T,1} where T<:QF
    mass::AbstractArray{T,1} where T<:QF
    density::AbstractArray{T,1} where T<:QF

    intact0::AbstractArray{Int64, 1}
    neighs::AbstractArray{Int64,2}

    boundary_conditions::AbstractArray{T, 1} where T
    short_range_contact::AbstractArray{T, 1} where T
    material_blocks::AbstractArray{T, 1} where T

    boundaries::Tuple
    Collect!::Function
    Params::Dict{Symbol, Any}
    Out::Dict{Symbol, Any}
    cprint::Function
end

function Base.show(io::IO, env::GeneralEnvironment)
    print(io, getPanel(env))
end

function getPanel(env::GeneralEnvironment, ptype=SPanel; width=Term.default_width())
    info2(x) = begin
        txt = @green string(x) * "\n"
        txt
    end
    info3(x) = begin
        txt = @green string(x) * "\n"
        txt
    end
    info1(x) = [Panel(
                info3(item),
                title= "$(item.name)",
                title_style="bold blue",
                width=60
                )
                for item in x]

    params = @bold(@yellow "Params key(s): ") * "$(keys(env.Params))"
    outt = @bold(@yellow "Out key(s): ") * "$(keys(env.Out))"
    name(j) = variable_color(j)
    txt = [
            Panel("$(name("EnvID")): $(env.id)",
            "$(name("State")): $(env.state)",
            "$(name("BID(s)")): $(unique(env.bid))",
            "$(name("Type(s)")): $(unique(env.type))",
            "$(name("Time Step")): $(env.time_step)",
            "$(name("dt")): $(env.dt)",
            "$(name("Momentum")): $(array_repr(sum(env.p, dims=2)))",
            "$(name("Force(acc)")): $(array_repr(sum(env.f, dims=2)))",
            params,
            outt,
            title="📜 Info",
            title_style="bold red",
            width=width-6),

            Panel([getPanel(mat; ptype=SPanel, width=width-12) for mat in env.material_blocks],
            title="🧱 Material Blocks", title_style="bold red", width=width-6),

            Panel([getPanel(cm; ptype=SPanel, width=width-12) for cm in env.short_range_contact],
            title="💥 Contact Models", title_style="bold red", width=width-6),

            Panel([getPanel(cm; ptype=SPanel, width=width-12) for cm in env.boundary_conditions],
            title="📍 Boundary Conditions", title_style="bold red", width=width-6),
        ]
    return ptype(txt,
        title="🗃  General Environment: $(env.id)",
        subtitle="End of General Environment: $(env.id)",
        subtitle_style="bold blue",
        subtitle_justify=:center,
        )
end



function _cudaconvert(x::Vector{GeneralEnvironment})
    _cudaconvert.(x)
end

function _cudaconvert(x::T) where T <: GeneralEnvironment
    function fn(x, k)
        if !(k in [:Out, :Params])
            _cudaconvert(getfield(x, k))
        else
            _cudaconvert(getfield(x, k))
        end
    end
    args = (fn(x, k) for k in fieldnames(T))
    T(args...)
end

"""
    Env(id::Int64, materials,short_range_contact, boundary_conds, dt; state=2, bskin=0.5, units=false)

Constructor for `GeneralEnvironment` type.

# Arguments
- `id::Int64`: Environment ID
- `materials::AbstractArray{Material, 1}`: Materials
- `short_range_contact::AbstractArray{ContactModel, 1}`: Short range contact models
- `boundary_conds::AbstractArray{BoundaryCondition, 1}`: Boundary conditions
- `dt::T where T<:QF`: Time step
- `state::Int64`: State of the environment
- `bskin::T where T<:QF`: Boundary skin
- `units::Bool`: Units

# Returns
- `env::GeneralEnvironment`: General environment
"""
function Env(id::Int64, materials, short_range_contact, boundary_conds, dt; state=2, bskin=0.5, units=false)
    mat = materials[1]
    type = mat.general.type
    intact = sum(mat.general.intact, dims=1)
    y = mat.general.y
    v = mat.general.velocity
    volume = mat.general.volume
    mass = volume * unit(eltype(mat.specific.density))
    bid = fill(mat.blockid, length(mass))

    for j in mat.type
        mask = (mat.general.type .== j)
        t = j - mat.type.start + 1
        mass[mask] .*= Unitful.ustrip(mat.specific.density[t])
    end

    for mat in materials[2:end]
        type = vcat(type, mat.general.type)
        y = hcat(y, mat.general.y)
        v = hcat(v, mat.general.velocity)
        volume = vcat(volume,mat.general.volume)
        intact = hcat(intact, sum(mat.general.intact, dims=1))
        mass_ = mat.general.volume * unit(eltype(mat.specific.density))
        for j in mat.type
            mask = (mat.general.type .== j)
            t = j - mat.type.start + 1
            mass_[mask] .*= Unitful.ustrip(mat.specific.density[t])
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

    ghost_particles = 0*type
    v = Unitful.ustrip(v) * unit(eltype(y)) / unit(dt)
    p = v * unit(eltype(mass))
    f = 0 * v / unit(dt)
    time_step = 0
    neighs = zeros(2,2)
    density = mass ./ volume

    env = GeneralEnvironment(id, type, bid, ghost_particles, state,
                time_step, dt,
                y, v, f, p,
                volume,
                mass,
                density,
                reshape(intact, :),
                neighs,
                boundary_conds,
                short_range_contact,
                materials,boundaries,
                (env, i)->(nothing),
                Dict(),
                Dict(),
                (x)->(nothing))
    if ~units
        env = ustrip(env)
    end
    return deviceconvert(env)
end

################################################################################
# unit conversion
################################################################################

function uconvert_to(DIMS, env::GeneralEnvironment)
    uconvert_ = (item) -> uconvert_to(DIMS, item)
    GeneralEnvironment(env.id, env.type, env.bid, env.ghost_particles, env.state,
           env.time_step, uconvert_(env.dt),

            uconvert_(env.y),
            uconvert_(env.v),
            uconvert_(env.f),
            uconvert_(env.p),
            uconvert_(env.volume),
            uconvert_(env.mass),
            uconvert_(env.density),
            env.intact0,
            env.neighs,

            uconvert_.(env.boundary_conditions),
            uconvert_.(env.short_range_contact),
            uconvert_.(env.material_blocks),

            env.boundaries,
            env.Collect!,
            env.Params,
            env.Out,
            env.cprint,
            )
end

function ustrip(env::GeneralEnvironment)
    GeneralEnvironment(env.id, env.type, env.bid, env.ghost_particles, env.state,
            env.time_step, ustrip(env.dt),

            ustrip(env.y), ustrip(env.v),
            ustrip(env.f), ustrip(env.p),
            ustrip(env.volume),
            ustrip(env.mass),
            ustrip(env.density),

            env.intact0,
            env.neighs,

            ustrip.(env.boundary_conditions),
            ustrip.(env.short_range_contact),
            ustrip.(env.material_blocks),

            env.boundaries,
            env.Collect!,
            ustrip(env.Params),
            ustrip(env.Out),
            env.cprint)
end

"""
    printdata(env)

Print (and save to log.txt) data of an environment during simulation.

If `PRINT_DEFAULT_DATA` is set to `true`, the following data will be printed:
- `envID`: Environment ID
- `time`: Simulation time
- `px`, `py`, `pz`: Total momentum in x, y, z directions
- `Fx`, `Fy`, `Fz`: Total force in x, y, z directions
- `damage`: Total damage
- and the data printed by `cprint` function of the environment

# Arguments
- `env`: Environment
"""
function printdata(env)
    if PRINT_DEFAULT_DATA[]
        log_data(envID=env.id)
        log_data(time=env.time_step*env.dt)
        p = sum(env.p, dims=2)
        log_data(px=p[1], py=p[2], pz=p[3])
        f = sum(env.f, dims=2)
        log_data(Fx=f[1], Fy=f[2], Fz=f[3])
        damage = cal_damage(env)
        log_data(damage=sum(damage)/length(damage))
    end
    env.cprint(env)
end

"""
    set_ghost_particles!(env,ghost_particles)

Set ghost particles for an environment.
"""
function set_ghost_particles!(env, ghost_particles)
    env.ghost_particles[1:end,1:end] = ghost_particles
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
