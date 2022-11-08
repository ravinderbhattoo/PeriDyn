"""
This module contain definition of Peridynamics material type.
"""

export GeneralMaterial, force_density


"""
Abstract PeridynamicsMaterial type.
"""
abstract type PeridynamicsMaterial end

@def PeridynamicsMaterial_gf begin
    name::String
    type::UnitRange{Int64}
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
    PeridynamicsMaterial(name, type, bid, gen, spc)
end

struct GeneralMaterial
    y::Array{Float64,2}
    velocity::Array{Float64,2}
    x::Array{Float64,2}
    volume::Array{Float64,1}
    type::Array{Int64,1}
    particle_size::Float64
    horizon::Float64
    family::Array{Int64,2}
    intact::BitArray{2}
    weighted_volume::Array{Float64,1}
    deformed::Vector{Bool}
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


abstract type SpecificMaterial end

function Base.show(io::IO, i::SpecificMaterial)
    println(io, typeof(i))
    for j in fieldnames(typeof(i))
        println(io, j, ": ", getproperty(i, j))
    end
end



function force_density_T(y::Array{Float64,2}, mat::PeridynamicsMaterial, ::Union{Type{Val{:cpu}}, Type{Val{:cuda}}}; kwargs...)
    force_density_T(y, mat)
end

function force_density_T(y::Array{Float64,2}, mat::PeridynamicsMaterial, device::Symbol)
    force_density_T(y::Array{Float64,2}, mat::PeridynamicsMaterial, Val{device})
end


function force_density(y::Array{Float64,2}, mat::PeridynamicsMaterial)
    device = DEVICE[]
    if mat.general.deformed[1]
        return force_density_T(y, mat, device)
    else
        y_ = y .- y[:, 1] .+ mat.general.x[:, 1]
        if isapprox(y_, mat.general.x)
            println("Material force calculation passed ($(mat.name)). [Undeformed]")
            return 0*y
        else
            mat.general.deformed[1] = true
            return force_density_T(y, mat, device)
        end
    end
end