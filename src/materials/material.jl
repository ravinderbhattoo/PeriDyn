"""
This module contain definition of Peridynamics material type.
"""

export GeneralMaterial, PeridynamicsMaterial, SpecificMaterial, force_density


"""
Abstract PeridynamicsMaterial type.
"""
abstract type PeridynamicsMaterial end

abstract type SpecificMaterial end

function init(spc::T, gen) where T <: SpecificMaterial
    spc
end

@def PeridynamicsMaterial_gf begin
    name::String
    type::Union{UnitRange{Int64}, AbstractVector{Int64}}
    blockid::Int64
    general::GeneralMaterial
end

"""
    PeridynamicsMaterial(gen, spc)

Create peridynamics material.
"""
function PeridynamicsMaterial(name, type, bid, gen, spc)
    error("Not specified for type **$(typeof(spc))**")
end


"""
    PeridynamicsMaterial(gen, spc::PeridynamicsMaterial; name="PM")
"""
function PeridynamicsMaterial(bid, gen, spc; name="PM")
    type = minimum(gen.type):maximum(gen.type)
    PeridynamicsMaterial(name, type, bid, gen, spc)
end


"""
    PeridynamicsMaterial(gen, spc::PeridynamicsMaterial; name="PM")
"""
function PeridynamicsMaterial(gen, spc; name="PM")
    bid = PDBlockID[]
    PDBlockID[] += 1
    type = minimum(gen.type):maximum(gen.type)
    PeridynamicsMaterial(name, type, bid, gen, init(spc, gen))
end

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
    GeneralMaterial(y0,v0,x,volume,type,horizon; particle_size=0,max_neigh=50)

General peridynamics material type.
"""
function GeneralMaterial(y0, v0, x, volume, type, horizon; max_neigh=100, particle_size=0, skip_bb=false)
    family = cal_family(x, horizon, max_neigh)
    intact = family .> 0
    println("Average family members: ", sum(intact)/size(intact, 2))
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
"""
function GeneralMaterial(out::Dict, args...; kwargs...)
    GeneralMaterial(unpack(out)..., args...; kwargs...)
end

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


function force_density_T(f, y, limits, mat::PeridynamicsMaterial, ::Union{Type{Val{:cpu}}, Type{Val{:cuda}}}; kwargs...)
    force_density_T(f, y, limits, mat)
end

function force_density_T(f, y, limits, mat::PeridynamicsMaterial, device::Symbol)
    force_density_T(f, y, limits, mat, Val{device})
end


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
            println("Material force calculation passed ($(mat.name)). [Undeformed]")
        else
            mat.general.deformed .= true
            force_density_T(f, y, limits, mat, device)
        end
    end
end



"""
    force_density_T(f, y, limits, mat::PeridynamicsMaterial, device::Type{Val{:cpu}}; particles=nothing)

Calculates force density (actually acceleration) for ordinary state based material type.
"""
function force_density_T(f_, y_, limits, mat::PeridynamicsMaterial, device::Type{Val{:cpu}}; particles=nothing)

    y = y_[:, limits[1]:limits[2]]
    numparticles = size(mat.general.family, 2)
    max_neighs = size(mat.general.family, 1)

    if isnothing(particles)
        N = 1:numparticles
        N_ = N .+ (limits[1] .- 1)
    else
        N_ = particles
        N = N_ .- (limits[1] .- 1)
    end
    M = max_neighs
    ARGS = map((i) -> (i, 1:M), N)

    out_of_loop!(y, mat, device)

    @inbounds function with_if_cal_force_ij(i, k)
        if mat.general.intact[k,i]
            j = mat.general.family[k,i]
            X = get_ij(j,i, mat.general.x)
            Y = get_ij(j,i, y)
            yij = get_magnitude(Y)
            xij = get_magnitude(X)
            extention = yij - xij
            s = extention/xij
            wij = influence_function(X)
            wji = influence_function(-X)
            type1 = mat.general.type[i] - mat.type.start + 1
            type2 = mat.general.type[j] - mat.type.start + 1
            if s < mat.specific.critical_stretch[type1, type2]
                return force_density_t_ij(mat, i, j, X, Y, xij, yij, extention, s, wij, wji, type1, type2)
            else
                mat.general.intact[k,i] = false
                return [0.0, 0.0, 0.0]
            end
        else
            return [0.0, 0.0, 0.0]
        end
    end

    inner_map(i, inds) = map_reduce((j)-> with_if_cal_force_ij(i,j), +, inds; init=[0.0, 0.0, 0.0])
    # inner_map(i, inds) = mapreduce((j)-> with_if_cal_force_ij(i,j), +, inds)
    outer_map(ARGS) = map((x)->inner_map(x[1], x[2]), ARGS)

    f_[:, N_] .+= hcat(outer_map(ARGS)...)
end


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

