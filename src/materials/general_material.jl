"""
This module contain definition of Peridynamics material type.
"""

export GeneralMaterial, PeridynamicsMaterial, SpecificMaterial
export force_density!, force_density_T!
export @newmaterial

"""
    PeridynamicsMaterial

    Abstract Peridynamics material type.
"""
abstract type PeridynamicsMaterial end

"""
    SpecificMaterial

    Abstract specific material type.
"""
abstract type SpecificMaterial end

"""
    GeneralMaterial

General peridynamics material type.

# Fields
- `y::AbstractArray{Float64,2}`: Deformed position of the material.
- `velocity::AbstractArray{Float64,2}`: Velocity of the material.
- `x::AbstractArray{Float64,2}`: Position of the material.
- `volume::AbstractArray{Float64,1}`: Volume of the material.
- `type::AbstractArray{Int64,1}`: Type of the material.
- `particle_size::Float64`: Particle size of the material.
- `horizon::Float64`: Horizon of the material.
- `family::VeryBigArray{Int64,2}`: Family of the material.
- `intact::VeryBigArray{Bool, 2}`: Intact of the material.
- `weighted_volume::AbstractArray{Float64,1}`: Weighted volume of the material.
- `deformed::AbstractVector{Bool}`: Deformed of the material.
- `skip_bb::Bool`: Skip bond based material (no idea why it is there).

"""
struct GeneralMaterial
    y::AbstractArray{T,2} where T <: QF
    velocity::AbstractArray{T,2} where T <: QF
    x::AbstractArray{T,2} where T <: QF
    volume::AbstractArray{T,1} where T <: QF
    type::AbstractArray{Int64,1}
    particle_size::T where T <: QF
    horizon::T where T <: QF
    family::VeryBigArray{Int64,2}
    intact::VeryBigArray{Bool, 2}
    weighted_volume::AbstractArray{T,1} where T <: QF
    deformed::AbstractVector{Bool}
    skip_bb::Bool
end



"""
    GeneralMaterial(y0, v0, x, volume, type, horizon; max_neigh=100, particle_size=0, skip_bb=false)

Creates a `GeneralMaterial` type.

# Arguments
- `y0::AbstractArray{Float64,2}`: Initial Deformed position of the material.
- `v0::AbstractArray{Float64,2}`: Initial velocity of the material.
- `x::AbstractArray{Float64,2}`: Initial position of the material.
- `volume::AbstractArray{Float64,1}`: Volume of the material.
- `type::AbstractArray{Int64,1}`: Type of the material.
- `horizon::Float64`: Horizon of the material.

# Keyword Arguments
- `max_neigh::Int64` = 100: Maximum number of neighbors.
- `particle_size::Float64` = 0: Particle size of the material.
- `skip_bb::Bool` = false: Skip bond based material.

# Returns
- `mat::GeneralMaterial`: General material type.
"""
function GeneralMaterial(y0, v0, x, volume, type, horizon; max_neigh=100, particle_size=0, skip_bb=false)
    if particle_size == zero(particle_size)
        particle_size = volume[1]^(1/3)
    end
    family = cal_family(x, horizon, max_neigh)
    intact = family .> 0
    log_info("Average family members: $(sum(intact)/size(intact, 2))")
    k = (size(intact,1)-maximum(sum(intact,dims=1)))::Int64
    intact = intact[end:-1:k,:]
    family = family[end:-1:k,:]
    m = weighted_volume(x, volume, particle_size, family, horizon)
    if y0 == x
        deformed = false
    else
        deformed = true
    end
    return GeneralMaterial(y0, v0, x, volume, type, particle_size, horizon, family, intact, m, [deformed], skip_bb)
end


"""
    GeneralMaterial(items::Dict, args...; kwargs...)

Creates a `GeneralMaterial` type.

# Arguments
- `items::Dict`: Dictionary of the fields of the `GeneralMaterial` type.
- `args...`: Arguments of the `GeneralMaterial` type.
- `kwargs...`: Keyword arguments of the `GeneralMaterial` type.

# Returns
- `mat::GeneralMaterial`: General material type.

# See also
- [`GeneralMaterial`](@ref)
- [`GeneralMaterial(y0, v0, x, volume, type, horizon; max_neigh=100, particle_size=0, skip_bb=false)`](@ref)
"""
function GeneralMaterial(items::Dict, args...; kwargs...)
    GeneralMaterial(unpack(items)..., args...; kwargs...)
end

"""
    init(spc::T1, gen::T2) where {T1<:SpecificMaterial, T2<:GeneralMaterial}

Initialize specific material type. It uses information from `GeneralMaterial` to create `SpecificMaterial` type.
This is default implementation of `init` function which returns `spc` without any change.

# Arguments
- `spc::T1`: Specific material type.
- `gen::T2`: General material type.

# Returns
- `spc::T1`: Initialized specific material type.
"""
function init(spc::T1, gen::T2) where {T1<:SpecificMaterial, T2<:GeneralMaterial}
    spc
end

"""
    force_density_t_ij(mat, args...; kwargs...)

Force density function for peridynamics material. This is default implementation of `force_density_t_ij` function

"""
function force_density_t_ij(mat::PeridynamicsMaterial, args...; kwargs...)
    error("force_density_t_ij not specified for type **$(typeof(mat))**.
    It should be implemented with the definition of $(typeof(mat)).")
end

# Generate macro for adding common fields to PeridynamicsMaterial
@def PeridynamicsMaterial_gf begin
    name::String
    type::Union{UnitRange{Int64}, Vector{Int64}}
    blockid::Int64
    general::GeneralMaterial
end


"""
    PeridynamicsMaterial(name, type, bid, gen, spc)

Creates a peridynamics material type with given name, type, block id, general material type
and specific material type.

# Arguments
- `name::String`: Name of the material.
- `type::Union{UnitRange{Int64}, AbstractVector{Int64}}`: Type of the material.
- `bid::Int64`: Block id of the material.
- `gen::GeneralMaterial`: General material type.
- `spc::SpecificMaterial`: Specific material type.

# Returns
- `mat::PeridynamicsMaterial`: Peridynamics material type.
"""
function PeridynamicsMaterial(name, type, bid, gen, spc)
    error("Not specified for type **$(typeof(spc))**")
end


"""
    PeridynamicsMaterial(bid, gen, spc; name="PeridynamicsMaterial")

Creates a peridynamics material type with given block id, general and
specific material type. `name` is optional.

# Arguments
- `bid::Int64`: Block id of the material.
- `gen::GeneralMaterial`: General material type.
- `spc::SpecificMaterial`: Specific material type.
- `name::String`: Name of the material.

# Returns
- `mat::PeridynamicsMaterial`: Peridynamics material type.

"""
function PeridynamicsMaterial(bid, gen, spc; name="PeridynamicsMaterial")
    type = minimum(gen.type):maximum(gen.type)
    PeridynamicsMaterial(name, type, bid, gen, spc)
end


"""
    PeridynamicsMaterial(gen, spc; name="PeridynamicsMaterial")

Creates a peridynamics material type with given general and specific material
type. `name` is optional.

# Arguments
- `gen::GeneralMaterial`: General material type.
- `spc::SpecificMaterial`: Specific material type.
- `name::String`: Name of the material.

# Returns
- `mat::PeridynamicsMaterial`: Peridynamics material type.

"""
function PeridynamicsMaterial(gen, spc; name="PeridynamicsMaterial")
    bid = PDBlockID[]
    PDBlockID[] += 1
    type = minimum(gen.type):maximum(gen.type)
    PeridynamicsMaterial(name, type, bid, gen, spc)
end


"""
    force_density!(f, y, limits, mat::PeridynamicsMaterial)

Calculates the force density of the material.
Calls `force_density_T` with `device` if all the particles are deformed.

# Arguments
- `f::AbstractArray`: Force of the material.
- `y::AbstractArray`: Deformed position of the material.
- `limits::AbstractVector{Int64}`: Limits of the material.
- `mat::PeridynamicsMaterial`: Peridynamics material type.

# Returns
- `nothing`: In-place modification of `f`.
"""
function force_density!(f, y, limits, mat::PeridynamicsMaterial)
    device = DEVICE[]

    if all(mat.general.deformed)
        force_density_T!(f, y, limits, mat, device)
    else
        mask = (collect(1:size(y, 2)) .>= limits[1]) .& (collect(1:size(y, 2)) .<= limits[2])
        y_ = Array(y[:, mask])
        x = Array(mat.general.x)
        if isapprox(y_ .- y_[:, 1] .+ x[:, 1], x)
            log_info("Material force calculation passed ($(mat.name)). [Undeformed]")
        else
            mat.general.deformed .= true
            force_density_T!(f, y, limits, mat, device)
        end
    end
end


"""
    force_density_T!(f, y, limits, mat::PeridynamicsMaterial, device::Symbol)

Calculates force density of material pair particles. Converts `device` to `Type{Val{:device}}` and calls `force_density_T`.

# Arguments
- `f::AbstractArray`: Force of the material.
- `y::AbstractArray`: Deformed position of the material.
- `limits::AbstractVector{Int64}`: Limits of the material.
- `mat::PeridynamicsMaterial`: Peridynamics material type.
- `device::Symbol`: Device

# Returns
- `nothing`: In-place modification of `f`.
"""
function force_density_T!(f, y, limits, mat::PeridynamicsMaterial, device::Symbol)
    force_density_T!(f, y, limits, mat, Val{device})
end


"""
    force_density_T!(f, y, limits, mat::GeneralMaterial, device::Type{Val{:cpu}}; particles=nothing)

Calculates force density of material pair particles using CPU. Calls `force_density_t_ij` for each pair of
particles in the material. The `force_density_t_ij` function is defined in the specific material type.

# Arguments
- `f::AbstractArray`: Force of the material.
- `y::AbstractArray`: Deformed position of the material.
- `limits::AbstractVector{Int64}`: Limits of the material.
- `mat::GeneralMaterial`: General material type.
- `device::Type{Val{:cpu}}`: Device
- `particles::Union{Nothing, AbstractVector{Int64}}`: Particles to calculate the force.

# Returns
- `nothing`: In-place modification of `f`.
"""
function force_density_T!(f, y, limits, mat::PeridynamicsMaterial, device::Type{Val{:cpu}}; particles=nothing)

    y_ = @view y[:, limits[1]:limits[2]]
    f_ = @view f[:, limits[1]:limits[2]]
    numparticles = size(mat.general.family, 2)
    max_neighs = size(mat.general.family, 1)

    if isnothing(particles)
        Ns = 1:numparticles
        N_ = Ns .+ (limits[1] .- 1)
    else
        N_ = particles
        Ns = N_ .- (limits[1] .- 1)
    end
    M = max_neighs

    @timeit timings "Out of force loop" out_of_loop!(y_, mat, device)

    # ARGS = map((i) -> (i, 1:M), Ns)
    # inner_map(i, inds) = map_reduce((j)-> force_ij(i, j, mat, y_, f_), +, inds; init=[0.0, 0.0, 0.0])
    # # inner_map(i, inds) = mapreduce((j)-> force_ij(i,j), +, inds)
    # outer_map(ARGS) = map((x)->inner_map(x[1], x[2]), ARGS)
    # f_ .= hcat(outer_map(ARGS)...)

    # nested loops for each particle with its neighbors
    intact = mat.general.intact
    family = mat.general.family
    start = mat.type.start
    x = mat.general.x
    critical_stretch = mat.specific.critical_stretch
    type = mat.general.type .- (start - 1)
    volume = mat.general.volume
    density = mat.specific.density
    items = get_items(mat)
    @timeit timings "In force loop" Threads.@threads for i in Ns
        _den_i = density[type[i]]
        for k in 1:M
            if intact[k, i]
                j = family[k, i]
                _den_j = density[type[j]]
                X = get_ij(j, i, x)
                Y = get_ij(j, i, y_)
                yij = get_magnitude(Y)
                xij = get_magnitude(X)
                extension = yij - xij
                s = extension/xij
                type1 = type[i]
                type2 = type[j]
                if s < critical_stretch[type1, type2]
                    wij = influence_function(X)
                    wji = influence_function(.-X)
                    tij = force_density_t_ij(mat,
                                        i, j, k, type1, type2,
                                        xij, yij, extension, s,
                                        wij, wji, items)
                    xx = 0.5 * tij / yij 
                    # directly convert to acceleration from force density
                    # force_density  = force / vol_i / vol_j
                    # force_i = force_density * vol_j * vol_i 
                    # acceleration_i = force_i / mass_i = force_density * vol_j * vol_i / mass_i
                    # acceleration_i = force_density * vol_j * / den_i
                    f_[1, i] += xx * Y[1] * volume[j] / _den_i
                    f_[2, i] += xx * Y[2] * volume[j] / _den_i
                    f_[3, i] += xx * Y[3] * volume[j] / _den_i
                    f_[1, j] -= xx * Y[1] * volume[i] / _den_j
                    f_[2, j] -= xx * Y[2] * volume[i] / _den_j
                    f_[3, j] -= xx * Y[3] * volume[i] / _den_j
                else
                    intact[k,i] = false
                end
            end
        end
    end
end

"""
    force_density_T!(f, y, limits, mat::GeneralMaterial, device::Type{Val{:cuda}}; particles=nothing)

Calculates force density of material pair particles using CUDA. Calls `force_density_t_ij` for each pair of
particles in the material. The `force_density_t_ij` function is defined in the specific material type.

# Arguments
- `f::AbstractArray`: Force of the material.
- `y::AbstractArray`: Deformed position of the material.
- `limits::AbstractVector{Int64}`: Limits of the material.
- `mat::GeneralMaterial`: General material type.
- `device::Type{Val{:cuda}}`: Device
- `particles::Union{Nothing, AbstractVector{Int64}}`: Particles to calculate the force.

# Returns
- `nothing`: In-place modification of `f`.
"""
function force_density_T!(f, y, limits, mat::PeridynamicsMaterial, device::Type{Val{:cuda}}; particles=nothing)
    error("CUDA version is not implemented yet.")
end


"""
    newmaterial(name, fields...)

Creates a new material type. The macro creates two subtypes, one for the `PeridynamicsMaterial`
and one for the `SpecificMaterial` with the given name. The name of the `SpecificMaterial` is the
given name with the suffix `Specific`. The name of the `PeridynamicsMaterial` is the given name
with the suffix `Material`.

For example, if the name is `Elastic`, the macro creates two types `ElasticMaterial` and
`ElasticSpecific`. The `ElasticMaterial` is a subtype of `PeridynamicsMaterial` and the
`ElasticSpecific` is a subtype of `SpecificMaterial`. The `ElasticSpecific` type is used to
store the specific parameters of the material.

```julia
@newmaterial(Elastic,
                "young_modulus::Matrix",
                "poisson_ratio::Matrix",
                "density::Vector")
```

The user is expected to define the `force_density_t_ij` function for the `ElasticSpecific` type.

# Arguments
- `name::Symbol`: Name of the material.
- `fields::Symbol`: Fields of the specific material.

"""
macro newmaterial(name, schema...)
    name = string(name)
    items = map(string, schema)
    if_array_in(x) = any([occursin(i, x) for i in
                    ["Array", "Vector", "Matrix"]])
    skip_if(x) = (~occursin("Int", x) &&
                         ~occursin("{", x))
    cleanit(x) = if_array_in(x) && skip_if(x) ? x*"{<:QF}" : x
    fields = map(x->Meta.parse(cleanit(x)),
                    items)
    spc_name = Symbol(name*"Specific")
    mat_name = Symbol(name*"Material")
    return esc(quote
        """
            $($spc_name)

        Specific material type for the $($name) material.
        """
        struct $(spc_name) <: SpecificMaterial
                $(fields...)
        end

        """
            $($mat_name)

        Material type for the $name material.

        # Arguments
        - `name::String`: Name of the material.
        - `type::Union{UnitRange{Int64}, Vector{Int64}}`: Type of the material.
        - `bid::Int64`: Block id of the material.
        - `gen::GeneralMaterial`: General material type.
        - `spc::$($spc_name)`: Specific material type.
        """
        struct $mat_name <: PeridynamicsMaterial
            name::String
            type::Union{UnitRange{Int64}, Vector{Int64}}
            blockid::Int64
            general::GeneralMaterial
            specific::$spc_name
            function BondBasedMaterial(name, type, bid, gen, spc)
                new(name, type, bid, gen, init(spc, gen))
            end
        end

        """
            PeridynamicsMaterial(name, type, bid, gen, spc::$($spc_name))

        Creates a new material of the $($mat_name) type.

        # Arguments
        - `name::String`: Name of the material.
        - `type::Union{UnitRange{Int64}, Vector{Int64}}`: Type of the material.
        - `bid::Int64`: Block id of the material.
        - `gen::GeneralMaterial`: General material type.
        - `spc::$($spc_name)`: Specific material type.
        """
        function PeridynamicsMaterial(name, type, bid, gen, spc::$spc_name)
            $mat_name(name, type, bid, gen, spc)
        end

     end)
end

# """
#     force_ij(i, k, mat, y_, f_)

#     Calculates the Fᵢⱼ of the material.

# # Arguments
# - `i::Int64`: Particle index.
# - `k::Int64`: Neighbor index.
# - `mat::PeridynamicsMaterial`: Peridynamics material type.
# - `y_::AbstractArray`: Deformed position of the material.
# - `f_::AbstractArray`: Force of the material.

# # Returns (may be in-place)
# - `f_::AbstractArray`: Force of the material.
# """

# @inline function force_ij(i, k, mat, y_, f_)
#     j = mat.general.family[k, i]
#     X = get_ij(j, i, mat.general.x)
#     Y = get_ij(j, i, y_)
#     yij = get_magnitude(Y)
#     xij = get_magnitude(X)
#     extension = yij - xij
#     s = extension/xij
#     wij = influence_function(X)
#     wji = influence_function(-X)
#     type1 = mat.general.type[i] - mat.type.start + 1
#     type2 = mat.general.type[j] - mat.type.start + 1
#     if s < mat.specific.critical_stretch[type1, type2]
#         tij =  force_density_t_ij(mat, i, j, X, Y, xij, yij, extension, s, wij, wji, type1, type2)
#         # # return tij
#         f_[:, i] .+= tij
#     else
#         mat.general.intact[k,i] = false
#         # return [0.0, 0.0, 0.0]
#     end
# end


# CUDA convert functions
function _cudaconvert(x::Vector{T}) where T <: PeridynamicsMaterial
    _cudaconvert.(x)
end

function _cudaconvert(x::T) where T <: PeridynamicsMaterial
    T((_cudaconvert(getfield(x, k)) for k in fieldnames(T))...)
end

function _cudaconvert(x::Vector{GeneralMaterial})
    _cudaconvert.(x)
end

function _cudaconvert(x::T) where T <: GeneralMaterial
    T((_cudaconvert(getfield(x, k)) for k in fieldnames(T))...)
end

function _cudaconvert(x::Vector{SpecificMaterial})
    _cudaconvert.(x)
end

function _cudaconvert(x::T) where T <: SpecificMaterial
    T((_cudaconvert(getfield(x, k)) for k in fieldnames(T))...)
end

function Base.show(io::IO, i::GeneralMaterial)
    print(io, getPanel(i))
end

function getPanel(i::GeneralMaterial; ptype=SPanel, width=Term.default_width())
    # txt = "type: $(unique(i.type))" * "\n"
    # txt = txt * "size: $(size(i.x))" * "\n"
    # txt = txt * "horizon: $(i.horizon)" * "\n"
    # txt = txt * "particle size: $(unique(i.particle_size))" * "\n"

    txt = ""
    # for j in fieldnames(typeof(i))
    for j in [:type, :x, :y, :velocity, :volume,
                :particle_size, :horizon, :weighted_volume, :deformed]
        item = getproperty(i, j)
        # check if iterable
        txt = txt * "$(variable_color(j)): "
        if isa(item, AbstractArray)
            txt = txt * array_repr(item)
        else
            txt = txt * "$(item)\n"
        end
    end

    txt = txt[1:end-1]
    return ptype(txt,
        title="General Material",
        justify=:left,
        width=width,
        )
end

function Base.show(io::IO, i::SpecificMaterial)
    print(io, getPanel(i))
end

function getPanel(i::SpecificMaterial; ptype=SPanel, width=Term.default_width())
    txt = ""
    for j in fieldnames(typeof(i))
        item = getproperty(i, j)
        # check if iterable
        txt = txt * "$(variable_color(j)): "
        if isa(item, AbstractArray)
            txt = txt * array_repr(item)
        else
            txt = txt * "$(item)\n"
        end
    end
    txt = txt[1:end-1]
    return ptype(txt,
        title=string(typeof(i)),
        justify=:left,
        width=width)
end

function Base.show(io::IO, mat::PeridynamicsMaterial)
    print(io, getPanel(mat))
end

function getPanel(mat::PeridynamicsMaterial; ptype=SPanel, width=Term.default_width())
    txt = [
            "$(variable_color("BlockID(s)")): $(unique(mat.blockid))",
            "$(variable_color("Type(s)")): $(unique(mat.type))",
            getPanel(mat.general; ptype=DPanel, width=width-6),
            getPanel(mat.specific; ptype=DPanel, width=width-6),
        ]
    out = ptype(txt,
                title="$(mat.name)",
                width=width)
    return out
end

#################################################
# Unit conversion functions
# Unit stripping
###############################################

function uconvert_to(DIMS, x::T) where T <: Union{GeneralMaterial,
            SpecificMaterial, PeridynamicsMaterial}
    args = (uconvert_to(DIMS, getfield(x, k)) for k in fieldnames(T))
    return T(args...)
end

function ustrip(x::T) where T <: Union{GeneralMaterial,
            SpecificMaterial, PeridynamicsMaterial}
    args = (ustrip(getfield(x, k)) for k in fieldnames(T))
    return T(args...)
end

