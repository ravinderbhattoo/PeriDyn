export SkipMaterial, SkipSpecific, force_density_T, PeridynamicsMaterial

"""
Specific skip materail type.
"""
struct SkipSpecific <: SpecificMaterial
    density::Array{Float64, 1}
end

"""
Specific skip materail type.
"""
function SkipSpecific()
    SkipSpecific([1.0])
end



"""
    SkipMaterial(type::UnitRange{Int64}, general::GeneralMaterial, specific::SkipSpecific)
"""
struct SkipMaterial <: PeridynamicsMaterial
    @PeridynamicsMaterial_gf 
    specific::SkipSpecific
end

"""
    PeridynamicsMaterial(gen, spc::SkipSpecific)
"""
function PeridynamicsMaterial(name, type, bid, gen, spc::SkipSpecific)
    SkipMaterial(name, type, bid, gen, spc)
end


"""
    force_density_T(y::Array{Float64,2},mat::SkipMaterial)

Calculates force density (actually acceleration) for skip material type.
"""
function force_density_T(y::Array{Float64,2}, mat::SkipMaterial)
    return 0*y
end

#
