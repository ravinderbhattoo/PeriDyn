"""
This module contain definition of Peridynamics material type.
"""

export GeneralMaterial, PeridynamicsMaterial, SpecificMaterial, force_density


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

function Base.show(io::IO, i::SpecificMaterial)
    println(io, typeof(i))
    for j in fieldnames(typeof(i))
        item = getproperty(i, j)
        # check if iterable
        if isa(item, AbstractArray) && mapreduce(*, +, size(item)) > 10
            mapreduce(*, +, size(item)) > 10
            println(io, j, ": ", typeof(item), size(item))
        else
            println(io, j, ": ", item)
        end
    end
end

"""
    GeneralMaterial

General peridynamics material type.

# Fields
- `y::AbstractArray{Float64,2}`: Displacement of the material.
- `velocity::AbstractArray{Float64,2}`: Velocity of the material.
- `x::AbstractArray{Float64,2}`: Position of the material.
- `volume::AbstractArray{Float64,1}`: Volume of the material.
- `type::AbstractArray{Int64,1}`: Type of the material.
- `particle_size::Float64`: Particle size of the material.
- `horizon::Float64`: Horizon of the material.
- `family::AbstractArray{Int64,2}`: Family of the material.
- `intact::AbstractArray{Bool, 2}`: Intact of the material.
- `weighted_volume::AbstractArray{Float64,1}`: Weighted volume of the material.
- `deformed::AbstractVector{Bool}`: Deformed of the material.
- `skip_bb::Bool`: Skip bond based material.

    `skip_bb` ? why do we need this?
"""
struct GeneralMaterial
    y::AbstractArray{Float64,2}
    velocity::AbstractArray{Float64,2}
    x::AbstractArray{Float64,2}
    volume::AbstractArray{Float64,1}
    type::AbstractArray{Int64,1}
    particle_size::Float64
    horizon::Float64
    family::AbstractArray{Int64,2}
    intact::AbstractArray{Bool, 2}
    weighted_volume::AbstractArray{Float64,1}
    deformed::AbstractVector{Bool}
    skip_bb::Bool
end


function Base.show(io::IO, i::GeneralMaterial)
    println(io, typeof(i))
    println(io, "type: ", unique(i.type))
    println(io, "size: ", size(i.x))
    println(io, "horizon: ", i.horizon)
    println(io, "particle size: ", unique(i.particle_size))
end

"""
    GeneralMaterial(y0, v0, x, volume, type, horizon; max_neigh=100, particle_size=0, skip_bb=false)

    Creates a `GeneralMaterial` type.

# Arguments
- `y0::AbstractArray{Float64,2}`: Initial displacement of the material.
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
    family = cal_family(x, horizon, max_neigh)
    intact = family .> 0
    log_info("Average family members: $(sum(intact)/size(intact, 2))")
    k = (size(intact,1)-maximum(sum(intact,dims=1)))::Int64
    intact = intact[k:end,:]
    family = family[k:end,:]
    if particle_size==0
        particle_size = volume[1]^(1/3)
    end
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

    Creates a `GeneralMaterial` type. `items` is a dictionary of the fields of the `GeneralMaterial` type
    typical comes from PDMaterialPoints package.

# Arguments
- `items::Dict`: Dictionary of the fields of the `GeneralMaterial` type.
- `args...`: Arguments of the `GeneralMaterial` type.
- `kwargs...`: Keyword arguments of the `GeneralMaterial` type.

# Returns
- `mat::GeneralMaterial`: General material type.

# see also
- `GeneralMaterial(y0, v0, x, volume, type, horizon; max_neigh=100, particle_size=0, skip_bb=false)`

"""
function GeneralMaterial(items::Dict, args...; kwargs...)
    GeneralMaterial(unpack(items)..., args...; kwargs...)
end

"""
    init(spc::T1, gen::T2) where {T1<:SpecificMaterial, T2<:GeneralMaterial}

Initialize specific material type. It uses information from `GeneralMaterial` to create `SpecificMaterial` type.

# Arguments
- `spc::T1`: Specific material type.
- `gen::T2`: General material type.

# Returns
- `spc::T1`: Initialized specific material type.
"""
function init(spc::T1, gen::T2) where {T1<:SpecificMaterial, T2<:GeneralMaterial}
    spc
end

# Generate macro for adding common fields to PeridynamicsMaterial
@def PeridynamicsMaterial_gf begin
    name::String
    type::Union{UnitRange{Int64}, AbstractVector{Int64}}
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
    PeridynamicsMaterial(bid, gen, spc; name="PeriMat")

    Creates a peridynamics material type with given block id, general and specific material type. `name` is optional.

# Arguments
- `bid::Int64`: Block id of the material.
- `gen::GeneralMaterial`: General material type.
- `spc::SpecificMaterial`: Specific material type.
- `name::String`: Name of the material.

# Returns
- `mat::PeridynamicsMaterial`: Peridynamics material type.

"""
function PeridynamicsMaterial(bid, gen, spc; name="PeriMat")
    type = minimum(gen.type):maximum(gen.type)
    PeridynamicsMaterial(name, type, bid, gen, spc)
end


"""
    PeridynamicsMaterial(gen, spc; name="PeriMat")

    Creates a peridynamics material type with given general and specific material type. `name` is optional.

# Arguments
- `gen::GeneralMaterial`: General material type.
- `spc::SpecificMaterial`: Specific material type.
- `name::String`: Name of the material.

# Returns
- `mat::PeridynamicsMaterial`: Peridynamics material type.

"""
function PeridynamicsMaterial(gen, spc; name="PeriMat")
    bid = PDBlockID[]
    PDBlockID[] += 1
    type = minimum(gen.type):maximum(gen.type)
    PeridynamicsMaterial(name, type, bid, gen, spc)
end

function force_density_T(f, y, limits, mat::PeridynamicsMaterial, ::Union{Type{Val{:cpu}}, Type{Val{:cuda}}}; kwargs...)
    force_density_T(f, y, limits, mat)
end

"""
    force_density_T(f, y, limits, mat::PeridynamicsMaterial, device::Symbol)

    Calculates the force of the material.

# Arguments
- `f::AbstractArray{Float64,2}`: Force of the material.
- `y::AbstractArray{Float64,2}`: Displacement of the material.
- `limits::AbstractVector{Int64}`: Limits of the material.
- `mat::PeridynamicsMaterial`: Peridynamics material type.
- `device::Symbol`: Device of the material.

# Returns (may be in-place)
- `f::AbstractArray{Float64,2}`: Force of the material.
"""
function force_density_T(f, y, limits, mat::PeridynamicsMaterial, device::Symbol)
    force_density_T(f, y, limits, mat, Val{device})
end


"""
    force_density_T(f, y, limits, mat::PeridynamicsMaterial)

    Calculates the force of the material.

# Arguments
- `f::AbstractArray{Float64,2}`: Force of the material.
- `y::AbstractArray{Float64,2}`: Displacement of the material.
- `limits::AbstractVector{Int64}`: Limits of the material.
- `mat::PeridynamicsMaterial`: Peridynamics material type.

# Returns (may be in-place)
- `f::AbstractArray{Float64,2}`: Force of the material.
"""
function force_density(f, y, limits, mat::PeridynamicsMaterial)
    device = DEVICE[]
    # force_density_T(f, y, limits, mat, device)

    if all(mat.general.deformed)
        force_density_T(f, y, limits, mat, device)
    else
        mask = (collect(1:size(y, 2)) .>= limits[1]) .& (collect(1:size(y, 2)) .<= limits[2])
        y_ = Array(y[:, mask])
        x = Array(mat.general.x)
        if isapprox(y_ .- y_[:, 1] .+ x[:, 1], x)
            log_info("Material force calculation passed ($(mat.name)). [Undeformed]")
        else
            mat.general.deformed .= true
            force_density_T(f, y, limits, mat, device)
        end
    end
end


"""
    force_density_T(f, y, limits, mat::GeneralMaterial, device::Type{Val{:cpu}}; particles=nothing)

    Calculates the force of the material.

# Arguments
- `f::AbstractArray{Float64,2}`: Force of the material.
- `y::AbstractArray{Float64,2}`: Displacement of the material.
- `limits::AbstractVector{Int64}`: Limits of the material.
- `mat::GeneralMaterial`: General material type.
- `device::Type{Val{:cpu}}`: Device of the material.
- `particles::Union{Nothing, AbstractVector{Int64}}`: Particles to calculate the force.

# Returns (may be in-place)
- `f::AbstractArray{Float64,2}`: Force of the material.
"""
function force_density_T(f, y, limits, mat::PeridynamicsMaterial, device::Type{Val{:cpu}}; particles=nothing)

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
    # ARGS = map((i) -> (i, 1:M), Ns)

    out_of_loop!(y_, mat, device)

    # inner_map(i, inds) = map_reduce((j)-> force_ij(i, j, mat, y_, f_), +, inds; init=[0.0, 0.0, 0.0])
    # # inner_map(i, inds) = mapreduce((j)-> force_ij(i,j), +, inds)
    # outer_map(ARGS) = map((x)->inner_map(x[1], x[2]), ARGS)
    # f_ .= hcat(outer_map(ARGS)...)

    # nested loops for each particle with its neighbors
    Threads.@threads for i in Ns
        for k in 1:M
            force_ij(i, k, mat, y_, f_)
        end
    end

end

"""
    force_ij(i, k, mat, y_, f_)

    Calculates the Fᵢⱼ of the material.

# Arguments
- `i::Int64`: Particle index.
- `k::Int64`: Neighbor index.
- `mat::PeridynamicsMaterial`: Peridynamics material type.
- `y_::AbstractArray`: Displacement of the material.
- `f_::AbstractArray`: Force of the material.

# Returns (may be in-place)
- `f_::AbstractArray`: Force of the material.
"""

@inline function force_ij(i, k, mat, y_, f_)
    if mat.general.intact[k,i]
        j = mat.general.family[k,i]
        X = get_ij(j,i, mat.general.x)
        Y = get_ij(j,i, y_)
        yij = get_magnitude(Y)
        xij = get_magnitude(X)
        extention = yij - xij
        s = extention/xij
        wij = influence_function(X)
        wji = influence_function(-X)
        type1 = mat.general.type[i] - mat.type.start + 1
        type2 = mat.general.type[j] - mat.type.start + 1
        if s < mat.specific.critical_stretch[type1, type2]
            tij =  force_density_t_ij(mat, i, j, X, Y, xij, yij, extention, s, wij, wji, type1, type2)
            # return tij
            f_[:, i] .+= tij
        else
            mat.general.intact[k,i] = false
            # return [0.0, 0.0, 0.0]
        end
    else
        # return [0.0, 0.0, 0.0]
    end
end


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

