export force_density_t_ij, BondBasedMaterial, BondBasedSpecific

"""
    BondBasedSpecific

Bond based specific material type.

# Fields
- `bond_stiffness::AbstractArray{T,2}`: Bond stiffness matrix.
- `bulk_modulus::AbstractArray{T,2}`: Bulk modulus matrix.
- `critical_stretch::AbstractArray{T,2}`: Critical stretch matrix.
- `density::AbstractArray{T,1}`: Density vector.
- `func::Function`: Bond force function.

# Constructor
```
BondBasedSpecific(K::AbstractArray, critical_stretch::AbstractArray, density::AbstractArray; horizon=nothing, func=nothing)
BondBasedSpecific(K::Real, critical_stretch::Real, density::Real; kwargs...)
```
"""
struct BondBasedSpecific <: SpecificMaterial
    bond_stiffness::AbstractArray{T1,2} where T1 <: QF
    bulk_modulus::AbstractArray{T2,2}   where T2 <: QF
    critical_stretch::AbstractArray{TN,2} where TN
    density::AbstractArray{T4,1} where T4 <: QF
    func::Function
end

"""
    BondBasedSpecific(K::Matrix, critical_stretch::Matrix, density::Vector; horizon=nothing, func=nothing)

Constructs `BondBasedSpecific` material type.

# Arguments
- `K::Vector`: Bulk modulus matrix.
- `critical_stretch::Vector`: Critical stretch matrix.
- `density::Vector`: Density vector.

# Keyword Arguments
- `horizon::Real`: Horizon.
- `func::Function`: Bond force function.

# Returns
- `BondBasedSpecific`: Bond based specific material type.
"""
function BondBasedSpecific(
            K::AbstractArray{T1,2} where T1 <: QF, 
            critical_stretch::AbstractArray{TN,2} where TN, 
            density::AbstractArray{T4,1} where T4 <: QF; 
            horizon=nothing, func=nothing)
    if isa(func, Nothing)
        func = bond_force
    end
    if isa(horizon, Nothing)
        bond_stiffness = zero(K)
        log_impinfo("Bond stiffness of K is used instead of 18K/πδ⁴. Please specify horizon.")
    else
        bond_stiffness = 18 * K ./ pi ./ horizon^4
        log_impinfo("Bond stiffness of 18K/πδ⁴ is used.")
    end
    BondBasedSpecific(bond_stiffness, K, critical_stretch, density, func)
end

"""
    BondBasedSpecific(K::Vector, critical_stretch::Vector, density::Vector; kwargs...)

Constructs `BondBasedSpecific` material type.

# Arguments
- `K::Vector`: Bulk modulus matrix.
- `critical_stretch::Vector`: Critical stretch matrix.
- `density::Vector`: Density vector.

# Keyword Arguments
- `horizon::Real`: Horizon.
- `func::Function`: Bond force function.

# Returns
- `BondBasedSpecific`: Bond based specific material type.
"""
function BondBasedSpecific(K::Vector, critical_stretch::Vector, density::Vector; kwargs...)
    BondBasedSpecific(make_matrix(K), make_matrix(critical_stretch), density; kwargs...)
end

"""
    BondBasedSpecific(K::Real, critical_stretch::Real, density::Real; kwargs...)

Constructs `BondBasedSpecific` material type.

# Arguments
- `K::Real`: Bulk modulus matrix.
- `critical_stretch::Real`: Critical stretch matrix.
- `density::Real`: Density vector.

# Keyword Arguments
- `horizon::Real`: Horizon.
- `func::Function`: Bond force function.

# Returns
- `BondBasedSpecific`: Bond based specific material type.
"""
function BondBasedSpecific(K::Real, critical_stretch::Real, density::Real; kwargs...)
    BondBasedSpecific([K], [critical_stretch], [density]; kwargs...)
end

"""
    init(spc::BondBasedSpecific, gen::GeneralMaterial)

Initializes `BondBasedSpecific` material type.

# Arguments
- `spc::BondBasedSpecific`: Bond based specific material type.
- `gen::GeneralMaterial`: General material type.

# Returns
- `BondBasedSpecific`: Bond based specific material type.
"""
function init(spc::BondBasedSpecific, gen::GeneralMaterial)
    horizon = gen.horizon
    BondBasedSpecific(spc.bulk_modulus, spc.critical_stretch, spc.density;
                        horizon=horizon, func=spc.func)
end

"""
    BondBasedMaterial <: PeridynamicsMaterial

Bond based material type.

# Fields
- `name::String`: Name of material.
- `type::Array`: Type of material particles.
- `bid::Int`: Material block id.
- `general::GeneralMaterial`: General material type.
- `specific::BondBasedSpecific`: Bond based specific material type.
"""
struct BondBasedMaterial <: PeridynamicsMaterial
    @PeridynamicsMaterial_gf
    # macro will fill the following fields
    # name::String
    # type::Array
    # bid::Int
    # general::GeneralMaterial
    specific::BondBasedSpecific
    function BondBasedMaterial(name, type, bid, gen, spc)
        new(name, type, bid, gen, init(spc, gen))
    end
end

"""
    PeridynamicsMaterial(name, type, bid, gen, spc::BondBasedSpecific)

Constructs `BondBasedMaterial` material type.
"""
function PeridynamicsMaterial(name, type, bid, gen, spc::BondBasedSpecific)
    BondBasedMaterial(name, type, bid, gen, spc)
end

"""
    out_of_loop!(y_, mat::BondBasedMaterial, device)

It is called once before the main loop for force calculation.
It does nothing for `BondBasedMaterial`.
"""
@inline function out_of_loop!(y_, mat::BondBasedMaterial, device)
    nothing
end

"""
    get_items(mat::BondBasedMaterial)

Returns items for force calculation.

# Arguments
- `mat::BondBasedMaterial`: Bond based material type.

# Returns
- `bond_stiffness::AbstractArray`: Bond stiffness.
"""
function get_items(mat::BondBasedMaterial)
    return mat.specific.bond_stiffness
end

"""
    bond_force(s, bond_stiffness)

Bond force function.

# Arguments
- `s::Real`: Bond stretch.
- `bond_stiffness::Real`: Bond stiffness.

# Returns
- `Real`: Bond force.
"""
@inline function bond_force(s, bond_stiffness)
    return s * bond_stiffness
end

"""
    force_density_t_ij(mat::BondBasedMaterial,
                    i, j, k, type1, type2,
                    xij, yij, extension, s,
                    wij, wji, items)

Calculates force density (actually acceleration) for bond based material type.

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
function force_density_t_ij(mat::BondBasedMaterial,
                    i, j, k, type1, type2,
                    xij, yij, extension, s,
                    wij, wji, items)
    bond_stiffness = items[1]
    t = bond_force(s, bond_stiffness)
    return t
end
 

# """
#     force_density_T!(force::AbstractArray{Float64,2}, y::AbstractArray{Float64,2}, limits, mat::BondBasedMaterial, ::Type{Val{:cuda}}; particles=nothing)

# Calculates force density (actually acceleration) for bond based material type.
# """
# function force_density_T!(force::AbstractArray{Float64,2}, y::AbstractArray{Float64,2},
#                             limits, mat::BondBasedMaterial, ::Type{Val{:cuda}};
#                             particles=nothing)

#     force = CUDA.CuArray(zeros(eltype(y), size(y)))
#     y = CUDA.CuArray(y)
#     types = CUDA.CuArray(mat.general.type)
#     tstart = mat.type.start
#     x = CUDA.CuArray(mat.general.x)
#     volume = CUDA.CuArray(mat.general.volume)
#     stiffness = CUDA.CuArray(mat.specific.bond_stiffness)
#     critical_stretch = CUDA.CuArray(mat.specific.critical_stretch)
#     skip_bb = CUDA.CuArray([mat.general.skip_bb])
#     func = mat.specific.func

#     if isnothing(particles)
#         _N = 1:size(family, 2)
#     else
#         _N = particles
#     end

#     function cal_force(y, x, force, _N, family, intact, stiffness, volume, types, critical_stretch, skip_bb, func)
#         index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
#         stride = blockDim().x * gridDim().x

#         for ind in index:stride:length(_N)
#             i = _N[ind]
#             for k in 1:size(family, 1)
#                 if intact[k, ind]
#                     j = family[k, ind]
#                     _X = sqrt((x[1,j]-x[1,i])^2 + (x[2,j]-x[2,i])^2 + (x[3,j]-x[3,i])^2)
#                     Y1 = y[1,j]-y[1,i]
#                     Y2 = y[2,j]-y[2,i]
#                     Y3 = y[3,j]-y[3,i]
#                     _Y = sqrt((Y1)^2 + (Y2)^2 + (Y3)^2)
#                     ext = _Y - _X
#                     s = ext/_X
#                     type1 = types[i] - tstart+ 1
#                     type2 = types[j] - tstart + 1
#                     if skip_bb[1] || ((type1!=-1) && (type2!=-1) && (s < critical_stretch[type1, type2]))
#                         t = func(s, stiffness[type1, type2])
#                         force[1, i] += t*volume[j] * Y1/_Y
#                         force[2, i] += t*volume[j] * Y2/_Y
#                         force[3, i] += t*volume[j] * Y3/_Y
#                     else
#                         intact[k,ind] = false
#                     end
#                 end
#             end
#         end
#         return nothing
#     end


#     if length(_N) > max_num

#         __N = _N[1:max_num]
#         intact = CUDA.CuArray(mat.general.intact[:, __N])
#         family = CUDA.CuArray(mat.general.family[:, __N])
#         _cuN = CUDA.CuArray(__N)

#         kernel = CUDA.@cuda launch=false cal_force(y, x, force, _cuN, family, intact, stiffness, volume, types, critical_stretch, skip_bb, func)
#         config = launch_configuration(kernel.fun)

#         nthreads = Base.min(length(__N), config.threads)
#         nblocks =  cld(length(__N), nthreads)
        
#         for i in 1:(1+(length(_N) // max_num))
#             n1 = 1 + (i-1)*max_num
#             n2 = min(max_num*i, length(_N))
#             __N = _N[n1:n2]
#             intact = CUDA.CuArray(mat.general.intact[:, __N])
#             family = CUDA.CuArray(mat.general.family[:, __N])
#             _cuN = CUDA.CuArray(__N)
#             CUDA.@sync kernel(y, x, force, _cuN, family, intact, stiffness, volume, types, critical_stretch, skip_bb, func; threads=nthreads, blocks=nblocks)
#             mat.general.intact .= Array(intact)
#         end
#     else
#         intact = CUDA.CuArray(mat.general.intact[:, _N])
#         family = CUDA.CuArray(mat.general.family[:, _N])
#         _cuN = CUDA.CuArray(_N)
#         CUDA.@sync kernel(y, x, force, _cuN, family, intact, stiffness, volume, types, critical_stretch, skip_bb, func; threads=nthreads, blocks=nblocks)
#         mat.general.intact .= Array(intact)
#     end
#     return  CUDA.Array(force)[:, _N]
# end



