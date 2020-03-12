struct AbstractMaterial
    y::Array{Float64,2}
    velocity::Array{Float64,2}
    x::Array{Float64,2}
    particle_size::Float64
    volume::Array{Float64,1}
    density::Float64
    horizon::Float64
    critical_stretch::Float64
    family::Array{Int64,2}
    intact::BitArray{2}
end

abstract type PeridynamicsMaterial end

function Material(y0,v0,x,volume,density,horizon,critical_stretch; particle_size=0,max_neigh=50)
    family = cal_family(x,horizon,max_neigh)
    intact = family.>0.5
    if particle_size==0
        particle_size = volume[1]^(1/3)
    end
    return deepcopy(AbstractMaterial(y0,v0,x,particle_size,volume,density, horizon,critical_stretch,family,intact))
end
