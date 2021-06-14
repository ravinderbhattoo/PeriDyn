"""
This module contain definition of Peridynamics material type.
"""

export GeneralMaterial

"""
Abstract PeridynamicsMaterial type.
"""
abstract type PeridynamicsMaterial end

"""
    PeridynamicsMaterial(gen, spc)

Create peridynamics material.
"""
function PeridynamicsMaterial(gen, spc)
    error("Not specified for type **$(typeof(spc))**")
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
end


function Base.show(io::IO, i::GeneralMaterial) 
    println(io, typeof(i))
    println(io, "type: ", unique(i.type))
    println(io, "horizon: ", i.horizon)
    println(io, "particle size: ", unique(i.particle_size))
end

"""
    GeneralMaterial(y0,v0,x,volume,type,horizon; particle_size=0,max_neigh=50)

General peridynamics material type.
"""
function GeneralMaterial(y0, v0, x, volume, type, horizon; max_neigh=100, particle_size=0)
    family = cal_family(x, horizon, max_neigh)
    intact = family .> 0.5
    println("Average family members: ", sum(intact)/size(intact, 2))
    k = (size(intact,1)-maximum(sum(intact,dims=1)))::Int64
    intact = intact[k:end,:]
    family = family[k:end,:]
    if particle_size==0
        particle_size = volume[1]^(1/3)
    end
    m = weighted_volume(x, volume, particle_size, family, horizon)
    return GeneralMaterial(y0, v0, x, volume, type, particle_size, horizon, family, intact, m)
end


abstract type SpecificMaterial end

function Base.show(io::IO, i::SpecificMaterial)
    println(io, typeof(i))
    for j in fieldnames(typeof(i))
        println(io, j, ": ", getproperty(i, j))
    end
end