export SkipMaterial, SkipSpecific, force_density_T, PeridynamicsMaterial

"""
Specific skip materail type.
"""
struct SkipSpecific <: SpecificMaterial
end

"""
    SkipMaterial(type::UnitRange{Int64}, general::GeneralMaterial, specific::SkipSpecific)
"""
struct SkipMaterial <: PeridynamicsMaterial
    type::UnitRange{Int64}
    general::GeneralMaterial
    specific::SkipSpecific
end

"""
    PeridynamicsMaterial(gen, spc::SkipSpecific)
"""
function PeridynamicsMaterial(gen, spc::SkipSpecific)
    type = minimum(gen.type):maximum(gen.type)
    SkipMaterial(type, gen, spc)
end


"""
    force_density_T(y::Array{Float64,2},mat::SkipMaterial)

Calculates force density (actually acceleration) for skip material type.
"""
function force_density_T(y::Array{Float64,2}, mat::SkipMaterial)
    return 0*y
end

#
