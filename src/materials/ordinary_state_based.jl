export force_density_t_ij, OrdinaryStateBasedMaterial, OrdinaryStateBasedSpecific

"""
    OrdinaryStateBasedSpecific <: SpecificMaterial

Ordinary state based specific material type.

# Fields
- `bulk_modulus::AbstractArray{<:QF, 2}`: Bulk modulus.
- `shear_modulus::AbstractArray{<:QF, 2}`: Shear modulus.
- `critical_stretch::AbstractArray{Float64, 2}`: Critical stretch.
- `density::AbstractArray{<:QF, 1}`: Density.
- `theta::AbstractArray{<:QF, 1}`: Dilatation.

# Returns
- `OrdinaryStateBasedSpecific`: Ordinary state based specific material.
"""
mutable struct OrdinaryStateBasedSpecific <: SpecificMaterial
    bulk_modulus::AbstractArray{<:QF, 2}
    shear_modulus::AbstractArray{<:QF, 2}
    critical_stretch::AbstractArray{Float64, 2}
    density::AbstractArray{<:QF, 1}
    theta::AbstractArray{<:QF, 1}
end

"""
    init(spc::OrdinaryStateBasedSpecific, gen::GeneralMaterial)

Initializes ordinary state based specific material type. Set `spc.theta` to zero.

# Arguments
- `spc::OrdinaryStateBasedSpecific`: Specific material.
- `gen::GeneralMaterial`: General material.

# Returns
- `OrdinaryStateBasedSpecific`: Specific material.
"""
function init(spc::OrdinaryStateBasedSpecific, gen::GeneralMaterial)
    spc.theta = zero(gen.volume) / unit(eltype(gen.volume))
    return spc
end

"""
    OrdinaryStateBasedSpecific(bulk_modulus, shear_modulus, critical_stretch, density)

Creates ordinary state based specific material type.

# Arguments
- `bulk_modulus::AbstractVector{<:QF}`: Bulk modulus.
- `shear_modulus::AbstractVector{<:QF}`: Shear modulus.
- `critical_stretch::AbstractVector{Float64}`: Critical stretch.
- `density::AbstractVector{<:QF}`: Density.

# Returns
- `OrdinaryStateBasedSpecific`: Ordinary state based specific material.
"""
function OrdinaryStateBasedSpecific(
                            bulk_modulus::AbstractVector{<:QF},
                            shear_modulus::AbstractVector{<:QF},
                            critical_stretch::AbstractVector{Float64},
                            density::AbstractVector{<:QF})
    return OrdinaryStateBasedSpecific(make_matrix(bulk_modulus), make_matrix(shear_modulus), make_matrix(critical_stretch), density, 0*density)
end

"""
    OrdinaryStateBasedMaterial <: PeridynamicsMaterial

Ordinary state based material type.

# Fields
- `name::String`: Name of material.
- `type::Array`: Type of material particles.
- `bid::Int`: Material block id.
- `general::GeneralMaterial`: General material type.
- `specific::OrdinaryStateBasedSpecific`: Specific material type.

# Returns
- `OrdinaryStateBasedMaterial`: Ordinary state based material.
"""
struct OrdinaryStateBasedMaterial <: PeridynamicsMaterial
    @PeridynamicsMaterial_gf
    # macro will fill the following fields
    # name::String
    # type::Array
    # bid::Int
    # general::GeneralMaterial
    specific::OrdinaryStateBasedSpecific
end


"""
    PeridynamicsMaterial(name, type, bid, gen, spc::OrdinaryStateBasedSpecific)

Creates ordinary state based material type.

# Arguments
- `name::String`: Name of material.
- `type::Int`: Type of material.
- `bid::Int`: Bond type.
- `gen::GeneralMaterial`: General material.
- `spc::OrdinaryStateBasedSpecific`: Specific material.

# Returns
- `OrdinaryStateBasedMaterial`: Ordinary state based material.
"""
function PeridynamicsMaterial(name, type, bid, gen, spc::OrdinaryStateBasedSpecific)
    OrdinaryStateBasedMaterial(name, type, bid, gen, init(spc, gen))
end

"""
    out_of_loop!(y_, mat::OrdinaryStateBasedMaterial, device)

Calculates dilatation for ordinary state based material type.

# Arguments
- `y_::AbstractArray`: Deformed positions.
- `mat::OrdinaryStateBasedMaterial`: Ordinary state based material.
- `device::Type{Val{:device}}`: Device type.

# Returns
- `nothing`
"""
function out_of_loop!(y_, mat::OrdinaryStateBasedMaterial, device)
    @timeit timings "dilatation" dilatation!(y_, mat, device)
end

"""
    get_items(mat::OrdinaryStateBasedMaterial)

Returns items for force calculation.

# Arguments
- `mat::BondBasedMaterial`: Ordinary state based material.

# Returns
- `Tuple`: Items for force calculation.
(`m`, `theta`, `K`, `G`) where
`m` is weighted volume,
`theta` is dilatation,
`K` is bulk modulus,
`G` is shear modulus.
"""
function get_items(mat::OrdinaryStateBasedMaterial)
    return mat.general.weighted_volume,
                    mat.specific.theta,
                    mat.specific.bulk_modulus,
                    mat.specific.shear_modulus
end

"""
    force_density_t_ij(mat::OrdinaryStateBasedMaterial,
                            i, j, k, type1, type2,
                            xij, yij, extension, s,
                            wij, wji, items)

Calculates force density (actually acceleration) for ordinary state based material type.

# Arguments
- `mat::BondBasedMaterial`: Ordinary state based material.
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
function force_density_t_ij(mat::OrdinaryStateBasedMaterial,
                    i, j, k, type1, type2,
                    xij, yij, extension, s,
                    wij, wji, items)
    m, theta, K, G = items
    Kij = K[type1, type2]
    Gij = G[type1, type2]
    a, b = wij/m[i], wji/m[j]
    t = (3*Kij-5*Gij) * (theta[i]*xij*a + theta[j]*xij*b) +
            15*Gij*extension*(a + b)
    return t
end





# """
#     force_density_T!(force::AbstractArray{Float64,2}, y::AbstractArray{Float64,2}, limits, mat::OrdinaryStateBasedMaterial, device::Type{Val{:cuda}}; particles=nothing)

# Calculates force density (actually acceleration) for bond based material type.
# """
# function force_density_T!(force::AbstractArray{Float64,2}, y::AbstractArray{Float64,2}, limits, mat::OrdinaryStateBasedMaterial, device::Type{Val{:cuda}}; particles=nothing)

#     types = mat.general.type
#     tstart = mat.type[1]
#     x = mat.general.x
#     intact = mat.general.intact
#     family = mat.general.family
#     volume = mat.general.volume
#     m = mat.general.weighted_volume
#     critical_stretch = mat.specific.critical_stretch
#     bulk_modulus = mat.specific.bulk_modulus
#     shear_modulus = mat.specific.shear_modulus

#     if isnothing(particles)
#         __N = collect(1:size(family, 2))
#         _N = __N .+ (limits[1] .- 1)
#     else
#         _N = particles
#         __N = _N .- (limits[1] .- 1)
#     end

#     theta = mat.specific.theta
#     l1 = limits[1] - 1
#     dilatation!(theta, l1, y, x, intact, family, volume, m,
#                 mat.general.particle_size, mat.general.horizon, device)


#     __N = CUDA.CuArray(__N)
#     _N = CUDA.CuArray(_N)


#     if isa(mat.specific.kernel, Nothing)

#         function cal_force(force, y, x, family, intact, volume, types, _N, __N,
#                             critical_stretch, theta, m, bulk_modulus, shear_modulus)
#             index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
#             stride = blockDim().x * gridDim().x

#             for ind in index:stride:length(_N)
#                 _i = _N[ind]
#                 __i = __N[ind]
#                 for k in 1:size(family, 1)
#                     if intact[k, __i]
#                         j = family[k, __i]
#                         _j = j + _i - __i

#                         xij = sqrt((x[1,j]-x[1,__i])^2 + (x[2,j]-x[2,__i])^2 + (x[3,j]-x[3,__i])^2)

#                         Y1 = y[1,_j] - y[1,_i]
#                         Y2 = y[2,_j] - y[2,_i]
#                         Y3 = y[3,_j] - y[3,_i]

#                         _Y = sqrt((Y1)^2 + (Y2)^2 + (Y3)^2)
#                         extension = _Y - xij
#                         s = extension/xij

#                         wij = 1/xij
#                         wji = 1/xij

#                         type1 = types[__i] - tstart + 1
#                         type2 = types[j] - tstart + 1
#                         if (type1!=-1) && (type2!=-1) && (s < critical_stretch[type1, type2])
#                             K = bulk_modulus[type1, type2]
#                             G = shear_modulus[type1, type2]
#                             a_, b_ = wij/m[__i], wji/m[j]
#                             t =  (3*K-5*G) * (theta[__i]*xij*a_ + theta[j]*xij*b_)  + 15*G*extension*(a_ + b_)
#                             t = t * volume[j]
#                             force[1, _i] += t * Y1/_Y
#                             force[2, _i] += t * Y2/_Y
#                             force[3, _i] += t * Y3/_Y
#                         else
#                             intact[k,__i] = false
#                         end
#                     end
#                 end
#             end
#             return nothing
#         end

#         kernel = CUDA.@cuda launch=false cal_force(force, y, x, family, intact, volume, types, _N, __N,
#                                             critical_stretch, theta, m, bulk_modulus, shear_modulus)
#         config = launch_configuration(kernel.fun)
#         nthreads = Base.min(length(_N), config.threads)
#         nblocks =  cld(length(_N), nthreads)

#         function fn(args...)
#             CUDA.@sync kernel(args...; threads=nthreads, blocks=nblocks)
#         end
#         mat.specific.kernel = fn
#     end


#     mat.specific.kernel(force, y, x, family, intact, volume, types, _N, __N,
#                         critical_stretch, theta, m, bulk_modulus, shear_modulus)
# end



