export force_density_t_ij, ElastoPlasticSolidMaterial,
        ElastoPlasticSolidSpecific
export PlasticFailureCriteria

"""
    PlasticFailureCriteria

Abstract type for plastic failure criteria.
"""
abstract type PlasticFailureCriteria end

"""
    ElasticPlasticSolidMaterial <: PeridynamicsMaterial

Elastic-plastic solid material.

# Fields
- `bulk_modulus::AbstractArray{<:QF,2}`: Bulk modulus
- `shear_modulus::AbstractArray{<:QF,2}`: Shear modulus
- `critical_stretch::AbstractArray{Float64,2}`: Critical stretch
- `density::AbstractArray{<:QF,1}`: Density
- `σγ::AbstractArray{<:QF,2}`: Yield stress
- `ψ::AbstractArray{<:QF,2}`: Yield stress value
- `edp::AbstractArray{<:QF,2}`: Plastic deviatoric strain
- `td_norm2::AbstractArray{<:QF,1}`: td_norm2
- `theta::AbstractArray{<:QF,1}`: Dilatation
- `criteria::PlasticFailureCriteria`: Plastic failure criteria
"""
mutable struct ElastoPlasticSolidSpecific <: SpecificMaterial
    bulk_modulus::AbstractArray{<:QF,2}
    shear_modulus::AbstractArray{<:QF,2}
    critical_stretch::AbstractArray{Float64,2}
    density::AbstractArray{<:QF,1}
    σγ::AbstractArray{<:QF,2}
    ψ::AbstractArray{<:QF,2}
    edp::AbstractArray{<:QF,2}
    td_norm2::AbstractArray{<:QF,1}
    theta::AbstractArray{<:QF,1}
    criteria::PlasticFailureCriteria
end


"""
    init(spc::ElastoPlasticSolidSpecific, gen::GeneralMaterial)

Initialize the specific material. It sets the yield stress value `ψ` equal to
`` 75/8π * \\sigma_e^2 / \\delta^5 `` where ``\\sigma_e`` is the effective yield stress
and ``\\delta`` is the horizon.
It initializes the dilatation `theta`, the plastic deviatoric strain `edp` and the
`td_norm2` to zero.

# Arguments
- `spc::ElastoPlasticSolidSpecific`: Specific material
- `gen::GeneralMaterial`: General material

# Returns
- `spc::ElastoPlasticSolidSpecific`: Initialized specific material
"""
function init(spc::ElastoPlasticSolidSpecific, gen::GeneralMaterial)
    δ = gen.horizon
    estress(x) = begin
        σ = [1.0, 0, 0] * x
        effective_stress2(σ..., spc.criteria)
    end
    σe2 = estress.(spc.σγ)
    spc.ψ = 75/8pi * σe2 / δ^5
    spc.edp = zeros(eltype(gen.x), size(gen.family)...)
    spc.td_norm2 = zeros(eltype(spc.ψ), size(gen.volume)...)
    spc.theta = zeros(typeof(1.0u"m"/u"m"), size(gen.volume)...)
    spc
end

"""
    ElastoPlasticSolidSpecific(
        bulk_modulus::AbstractArray{<:QF,1},
        shear_modulus::AbstractArray{<:QF,1},
        critical_stretch::Array{Float64,1},
        density::AbstractArray{<:QF,1},
        σγ::AbstractArray{<:QF,1};
        criteria = VonMises())

Construct an `ElastoPlasticSolidSpecific` material.

# Arguments
- `bulk_modulus::AbstractArray{<:QF,1}`: Bulk modulus
- `shear_modulus::AbstractArray{<:QF,1}`: Shear modulus
- `critical_stretch::Array{Float64,1}`: Critical stretch
- `density::AbstractArray{<:QF,1}`: Density
- `σγ::AbstractArray{<:QF,1}`: Yield stress
- `criteria::PlasticFailureCriteria`: Plastic failure criteria

# Returns
- `spc::ElastoPlasticSolidSpecific`: Specific material
"""
function ElastoPlasticSolidSpecific(
    bulk_modulus::AbstractArray{<:QF,1},
    shear_modulus::AbstractArray{<:QF,1},
    critical_stretch::Array{Float64,1},
    density::AbstractArray{<:QF,1},
    σγ::AbstractArray{<:QF,1};
    criteria = VonMises())
    return ElastoPlasticSolidSpecific(
        make_matrix(bulk_modulus), # bulk modulus
        make_matrix(shear_modulus), # shear modulus
        make_matrix(critical_stretch), # critical stretch
        density, # density
        make_matrix(σγ), # σγ
        make_matrix(σγ), # ψ
        make_matrix(zero(σγ)), # edp
        Float64[], # td_norm2
        Float64[], # theta
        criteria, # failure criteria
    )
end

"""
    ElastoPlasticSolidMaterial <: PeridynamicsMaterial

Elastic-plastic solid material.

# Fields
- `name::String`: Name of the material
- `type::Array`: Type of the material
- `bid::Int`: Block id
- `general::GeneralMaterial`: General material
- `specific::ElastoPlasticSolidSpecific`: Specific material
"""
struct ElastoPlasticSolidMaterial <: PeridynamicsMaterial
    @PeridynamicsMaterial_gf
    # macro will fill the following fields
    # name::String
    # type::Array
    # bid::Int
    # general::GeneralMaterial
    specific::ElastoPlasticSolidSpecific
end

"""
    PeridynamicsMaterial(name, type, bid, gen, spc::ElastoPlasticSolidSpecific)

Construct an `ElastoPlasticSolidMaterial`.

# Arguments
- `name::String`: Name of the material
- `type::Array`: Type of the material
- `bid::Int`: Block id
- `gen::GeneralMaterial`: General material
- `spc::ElastoPlasticSolidSpecific`: Specific material

# Returns
- `mat::ElastoPlasticSolidMaterial`: Material
"""
function PeridynamicsMaterial(name, type, bid, gen, spc::ElastoPlasticSolidSpecific)
    ElastoPlasticSolidMaterial(name, type, bid, gen, init(spc, gen))
end


"""
    out_of_loop!(y_, mat::ElastoPlasticSolidMaterial, device)

Calculate the dilatation and the td_norm2 out of the loop.

# Arguments
- `y_::AbstractArray{<:QF,2}`: Displacement
- `mat::ElastoPlasticSolidMaterial`: Material
- `device`: Device
"""
function out_of_loop!(y_, mat::ElastoPlasticSolidMaterial, device)
    @timeit timings "dilatation" dilatation!(y_, mat, device) # calculate dilatation
    fill!(mat.specific.td_norm2, zero(eltype(mat.specific.td_norm2)))
    @timeit timings "td_norm2" loop_over_neighs!(y_, mat, td_norm2_ij!, device) # calculate td_norm2
end


"""
    get_items(mat::ElastoPlasticSolidMaterial)

Get the items of the material.

# Arguments
- `mat::ElastoPlasticSolidMaterial`: Material

# Returns
- `items`: Tuple of items
(`weighted_volume`, `theta`, `td_norm2`, `bulk_modulus`, `shear_modulus`, `σγ`, `ψ`, `edp`, `horizon`, `volume`, `criteria`)
"""
function get_items(mat::ElastoPlasticSolidMaterial)
    return mat.general.weighted_volume,
                    mat.specific.theta,
                    mat.specific.td_norm2,
                    mat.specific.bulk_modulus,
                    mat.specific.shear_modulus,
                    mat.specific.σγ,
                    mat.specific.ψ,
                    mat.specific.edp,
                    mat.general.horizon,
                    mat.general.volume,
                    mat.specific.criteria
end


"""
    force_density_t_ij(mat::ElastoPlasticSolidMaterial,
        i, j, k, type1, type2,
        xij, yij, extension, s,
        wij, wji, items)

Calculate the force density.

# Arguments
- `mat::BondBasedMaterial`: Bond based material type.
- `i::Int`: Index of particle ``i``.
- `j::Int`: Index of particle ``j \\in N_{i}``.
- `k::Int`: Location index of particle ``j`` in `family` and `intact`.
- `type1::Int`: Type of particle ``i``.
- `type2::Int`: Type of particle ``j``.
- `xij::Real`: ``\\lVert X_{ij} \\rVert``
- `yij::Real`: ``\\lVert Y_{ij} \\rVert``
- `extension::Real`: Bond extension.
- `s::Real`: Bond stretch.
- `wij::Real`: Influence function (``\\omega_{ij}``) for pair ``(i,j)``.
- `wji::Real`: Influence function (``\\omega_{ji}``) for pair ``(j,i)``.
- `items::Tuple`: Items for force calculation.

# Returns
- `Real`: Force density.
"""
function force_density_t_ij(mat::ElastoPlasticSolidMaterial,
        i, j, k, type1, type2,
        xij, yij, extension, s,
        wij, wji, items)
    m, theta, td_norm2, K, G, σγ, ψ, edp, δ, vol, criteria = items
    Kij = K[type1, type2]
    Gij = G[type1, type2]
    ψij = ψ[type1, type2]

    C = sqrt(2 * ψij / td_norm2[i])

    a, b = wij/m[i], wji/m[j]
    θω_m = (theta[i] * a + theta[j] * b)
    theta_ij = θω_m / (a + b)
    ed = extension - theta_ij * xij / 3

    ede = ed - edp[k, i]
    td_trial = 15 * Gij * (a + b) * ede
    if C > 1
        td = td_trial
    else
        td = C * td_trial
        edp[k, i] = edp[k, i] + ede * (1 - C)
    end
    tn = (3 * Kij) * xij * θω_m
    t = tn + td
    return t
end

"""
    td_norm2_ij!(mat, i, j, k, type1, type2,
        xij, yij, extension, s,
        wij, wji, items)

Calculate the td_norm2.

# Arguments
- `mat::BondBasedMaterial`: Bond based material type.
- `i::Int`: Index of particle ``i``.
- `j::Int`: Index of particle ``j \\in N_{i}``.
- `k::Int`: Location index of particle ``j`` in `family` and `intact`.
- `type1::Int`: Type of particle ``i``.
- `type2::Int`: Type of particle ``j``.
- `xij::Real`: ``\\lVert X_{ij} \\rVert``
- `yij::Real`: ``\\lVert Y_{ij} \\rVert``
- `extension::Real`: Bond extension.
- `s::Real`: Bond stretch.
- `wij::Real`: Influence function (``\\omega_{ij}``) for pair ``(i,j)``.
- `wji::Real`: Influence function (``\\omega_{ji}``) for pair ``(j,i)``.
- `items::Tuple`: Items for force calculation.

# Returns
- 'nothing': Update `td_norm2` in-place.
"""
function td_norm2_ij!(mat, i, j, k, type1, type2,
                        xij, yij, extension, s,
                        wij, wji, items)
    m, theta, td_norm2, K, G, σγ, ψ, edp, δ, vol, criteria = items
    Gij = G[type1, type2]
    edpij = edp[k, i]
    a, b = wij/m[i], wji/m[j]
    ω_m = (a + b)
    θω_m = (theta[i] * a + theta[j] * b)
    theta_ij = θω_m / ω_m
    ed = extension - theta_ij * xij / 3
    td = 15Gij * (ed - edpij) * ω_m
    td_norm2[i] += td^2 * vol[j]
end


# =======================================
#       ustrip
# =======================================

function uconvert_to_default(spc::T) where T <: PlasticFailureCriteria
    uconvert_to(DEFAULT_UNITS, spc)
end

function uconvert_to_si(spc::T) where T <: PlasticFailureCriteria
    uconvert_to(SI_UNITS, spc)
end

function uconvert_to(DIMS, spc::T) where T <: PlasticFailureCriteria
    args = (uconvert_to(DIMS, getfield(spc, k)) for k in fieldnames(T))
    return T(args...)
end

function ustrip(x::T) where T<:PlasticFailureCriteria
    x
end

function ustrip_to_default(x::T) where T<:PlasticFailureCriteria
    ustrip(uconvert_to_default(DEFAULT_UNITS, x))
end

function ustrip_to_si(x::T) where T<:PlasticFailureCriteria
    ustrip(uconvert_to_si(SI_UNITS, x))
end

# =======================================
#       FAILURE CRITERIA
# =======================================

export VonMises, DruckerPrager, MohrCoulomb, Tresca
export effective_stress2

function effective_stress2(σ₁, σ₂, σ₃, criteria::T) where T <: PlasticFailureCriteria
    error("effective_stress not implemented for $(T)")
end

"""
    VonMises <: PlasticFailureCriteria

Von Mises failure criterion.
"""
struct VonMises <: PlasticFailureCriteria end

"""
    effective_stress2(σ₁, σ₂, σ₃, criteria::VonMises)

Calculate the effective stress for Von Mises failure criterion.

```math
σₑ^2 = ((σ₁ - σ₂)^2 + (σ₂ - σ₃)^2 + (σ₁ - σ₃)^2) / 2
```


# Arguments
- `σ₁::Real`: First principal stress.
- `σ₂::Real`: Second principal stress.
- `σ₃::Real`: Third principal stress.
- `criteria::VonMises`: Von Mises failure criterion.

# Returns
- `Real`: Effective stress.
"""
function effective_stress2(σ₁, σ₂, σ₃, criteria::VonMises)
    return @. ((σ₁ - σ₂)^2 + (σ₂ - σ₃)^2 + (σ₁ - σ₃)^2) / 2
end

"""
    DruckerPrager <: PlasticFailureCriteria

Drucker-Prager failure criterion.
"""
struct DruckerPrager <: PlasticFailureCriteria end

"""
    MohrCoulomb <: PlasticFailureCriteria

Mohr-Coulomb failure criterion.
"""
struct MohrCoulomb <: PlasticFailureCriteria end

"""
    Tresca <: PlasticFailureCriteria

Tresca failure criterion.
"""
struct Tresca <: PlasticFailureCriteria end


#



















# """
#     force_density_T!(y::Array{Float64,2},mat::ElastoPlasticSolidMaterial)

# Calculates force density (actually acceleration) for ordinary state based material type.
# """
# function force_density_T!(y::Array{Float64,2}, mat::ElastoPlasticSolidMaterial)
#     types = mat.general.type
#     force = zeros(size(y))
#     X = Float64[1.0, 0.0, 0.0]
#     Y = Float64[1.0, 0.0, 0.0]
#     S = mat.general
#     m = mat.general.weighted_volume
#     theta = dilatation(y, S, m)
#     td_norm2 = td_norm2_fn(y, theta, mat)
#     δ = S.horizon
#     edp = mat.specific.edp
#     N = size(S.x, 2)
#     _max = 0.0
#     for i = 1:N
#         for k = 1:size(S.family, 1)
#             j = S.family[k, i]::Int64
#             if !S.intact[k, i] || i == j
#                 # force[1,i] += 0.0
#                 # force[2,i] += 0.0
#                 # force[3,i] += 0.0
#                 S.intact[k, i] = false
#             else
#                 type1 = types[i] - mat.type.start + 1
#                 type2 = types[j] - mat.type.start + 1
#                 K = mat.specific.bulk_modulus[type1, type2]
#                 G = mat.specific.shear_modulus[type1, type2]
#                 σ = mat.specific.σγ[type1, type2]
#                 X[1], X[2], X[3] =
#                     S.x[1, j] - S.x[1, i], S.x[2, j] - S.x[2, i], S.x[3, j] - S.x[3, i]
#                 Y[1], Y[2], Y[3] = y[1, j] - y[1, i], y[2, j] - y[2, i], y[3, j] - y[3, i]
#                 _Y, _X = @_magnitude(Y), @_magnitude(X)
#                 xij = @_magnitude(X)::Float64
#                 wij = influence_function(X)::Float64
#                 wji = influence_function(-X)::Float64

#                 e = _Y - _X
#                 ed = e - theta[i] * _X / 3
#                 ee = ed - edp[k, i]

#                 td_trial = 15 * G * (ee * wij / m[i] + ee * wji / m[j])

#                 pyf, ψ = plastic_yield_function(td_norm2[i], σ, δ)

#                 if pyf < 0
#                     td = td_trial
#                 else
#                     td = sqrt(2 * ψ) * td_trial / td_norm2[i]
#                     td_extra = td_trial - td
#                     plastic_ed = td_extra / (15 * G * (wij / m[i] + wji / m[j]))
#                     edp[k, i] = edp[k, i] + plastic_ed
#                 end
#                 if (e / xij) < mat.specific.critical_stretch[type1, type2]
#                     t =
#                         (3 * K) *
#                         (theta[i] * xij * wij / m[i] + theta[j] * xij * wji / m[j])
#                     t = t + td
#                     t = t * S.volume[i] / @_magnitude(Y)
#                     force[1, i] += t * Y[1]
#                     force[2, i] += t * Y[2]
#                     force[3, i] += t * Y[3]
#                 else
#                     # force[1,i] += 0.0
#                     # force[2,i] += 0.0
#                     # force[3,i] += 0.0
#                     S.intact[k, i] = false
#                 end
#             end
#         end
#     end
#     return force
# end


# """
#     _tij(x,mi,thetai,eij,K,G)::Float64

# calculates tij force density magnitude (actually acceleration).
# """
# function _tij(x, mi, thetai, eij, K, G)::Float64
#     xij = @_magnitude(x)::Float64
#     wij = influence_function(x)::Float64
#     return ((3 * K - 5 * G) * (thetai * xij) + 15 * G * eij) * wij / mi
# end


# """
#     td_norm2_fn
# """
# function td_norm2_fn(y, theta, mat)
#     S = mat.general
#     types = mat.general.type
#     intact = S.intact
#     family = S.family
#     x = S.x
#     N = size(family, 2)
#     M = size(family, 1)
#     ARGS = map((i) -> (i, 1:M), 1:N)
#     m = mat.general.weighted_volume
#     function with_if_cal_td_2(i, k)
#         if intact[k, i]
#             j = family[k, i]
#             type1 = types[i] - mat.type.start + 1
#             type2 = types[j] - mat.type.start + 1
#             G = mat.specific.shear_modulus[type1, type2]
#             edp = mat.specific.edp[k, i]

#             X = [x[1, j] - x[1, i], x[2, j] - x[2, i], x[3, j] - x[3, i]]
#             Y = [y[1, j] - y[1, i], y[2, j] - y[2, i], y[3, j] - y[3, i]]

#             _X = @_magnitude(X)
#             _Y = @_magnitude(Y)

#             e = _Y - _X

#             ed = e - theta[i] * _X / 3

#             ee = ed - edp

#             s = e / _X

#             wij = influence_function(X)::Float64
#             wji = influence_function(-X)::Float64
#             type1 = types[i] - mat.type.start + 1
#             type2 = types[j] - mat.type.start + 1
#             if s < mat.specific.critical_stretch[type1, type2]
#                 K = mat.specific.bulk_modulus[type1, type2]
#                 G = mat.specific.shear_modulus[type1, type2]
#                 td = 15 * G * (ee * wij / m[i] + ee * wji / m[j])
#                 return td^2
#             else
#                 intact[k, i] = false
#                 return 0.0
#             end
#         else
#             return 0.0
#         end
#     end

#     inner_map(i, inds) = map_reduce((j)-> with_if_cal_force_ij(i,j), +, inds)

#     # inner_map(i, inds) = mapreduce((j)-> with_if_cal_force_ij(i,j), +, inds)
#     outer_map(ARGS) = map((x) -> inner_map(x[1], x[2]), ARGS)
#     return outer_map(ARGS)
# end

