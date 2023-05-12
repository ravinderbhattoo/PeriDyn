export SkipMaterial, SkipSpecific, force_density_T, PeridynamicsMaterial

"""
Specific skip materail type.
"""
struct SkipSpecific <: SpecificMaterial
    density::Array{AbstractFloat, 1}
    critical_stretch::Array{AbstractFloat, 2}
end

"""
Specific skip materail type.
"""
function SkipSpecific(density)
    SkipSpecific(density, Inf.*(ones(eltype(density), length(density), length(density))))
end

SkipSpecific() =  SkipSpecific([1.0])


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
    force_density_T(mat::SkipMaterial)

Calculates force density (actually acceleration) for bond based material type.
"""
function force_density_T(f, y::AbstractArray{AbstractFloat,2}, limits, mat::SkipMaterial; kwargs...)
    force_density_T(f, y, limits, mat, :cpu; kwargs)
end

function out_of_loop!(y_, mat::SkipMaterial, device)
    nothing
end

function force_density_t_ij(mat::SkipMaterial, i, j, X, Y, xij, yij, extention, s, wij, wji, type1, type2)
    return 0.0 .* Y
end

function force_density(f, y, limits, mat::SkipMaterial)
    println("Material force calculation passed ($(mat.name)). [SkipSpecific]")
end
