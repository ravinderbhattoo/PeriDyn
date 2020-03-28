abstract type PeridynamicsMaterial end

struct GeneralMaterial
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

function GeneralMaterial(y0,v0,x,volume,density,horizon,critical_stretch; particle_size=0,max_neigh=50)
    family = cal_family(x,horizon,max_neigh)
    intact = family.>0.5
    k = (size(intact,1)-maximum(sum(intact,dims=1)))::Int64
    intact = intact[k:end,:]
    family = family[k:end,:]
    if particle_size==0
        particle_size = volume[1]^(1/3)
    end
    return deepcopy(GeneralMaterial(y0,v0,x,particle_size,volume,density, horizon,critical_stretch,family,intact))
end
