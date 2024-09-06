export SkipMaterial, SkipSpecific, force_density

# @newmaterial Skip "density::Array"

# """
#     force_density!(f, y, limits, mat::SkipMaterial)

# Overload of the function `force_density` for the material `SkipMaterial` to skip
# froce calculations. Ideally force_density_t_ij should be defined for force density
# calculations for every `SpecificMaterial`.

# # Arguments
# - `f::Array`: Force vector.
# - `y::Array`: Displacement vector.
# - `limits::Array`: Limits of the domain.
# - `mat::SkipMaterial`: Material to calculate the force density.

# # Returns
# - `nothing`
# """
# function force_density!(f, y, limits, mat::SkipMaterial)
#     log_info("Material force calculation passed ($(mat.name)). [SkipSpecific]")
# end



"""
    SkipSpecific <: SpecificMaterial

Skip specific material type.

# Fields
- `density::AbstractArray{<:QF, 1}`: Density.
- `critical_stretch::AbstractArray{Float64, 2}`: Critical stretch.

# Returns
- `SkipSpecific`: Skip specific material.
"""
struct SkipSpecific <: SpecificMaterial
    density::Array{T, 1} where T
    critical_stretch::Array{T, 2} where T
end


function SkipSpecific(density)
    SkipSpecific(density, Inf.*(ones(eltype(density), length(density), length(density))))
end

SkipSpecific() =  SkipSpecific([1.0])

struct SkipMaterial <: PeridynamicsMaterial
    @PeridynamicsMaterial_gf
    specific::SkipSpecific
end


function PeridynamicsMaterial(name, type, bid, gen, spc::SkipSpecific)
    SkipMaterial(name, type, bid, gen, spc)
end


function force_density!(f, y, limits, mat::SkipMaterial)
    log_info("Material force calculation passed ($(mat.name)). [SkipSpecific]")
end


