# export PeridynamicsMaterial, ViscoElasticSolidMaterial,
#         ViscoElasticSolidSpecific, ViscoElasticSolidSpecificFull

# mutable struct ViscoElasticSolidSpecific <: SpecificMaterial
#     bulk_modulus::Array{Float64,2}
#     shear_modulus::Array{Float64,2}
#     critical_stretch::Array{Float64,2}
#     density::Array{Float64,1}
#     σγ::Array{Float64,2}
#     ψ::Array{Float64,2}
#     edp::Array{Float64,2}
#     td_norm2::Array{Float64,1}
#     theta::Array{Float64,1}
#     criteria::PlasticFailureCriteria
# end

# function init(spc::ViscoElasticSolidSpecific, gen::GeneralMaterial)
#     δ = gen.horizon
#     σe2 = effective_stress2(spc.σγ, 0.0, 0.0, spc.criteria)
#     spc.ψ = 75/8pi * σe2 / δ^5
#     spc.edp = zero(gen.family)
#     spc.td_norm2, spc.theta = zero(gen.volume), zero(gen.volume)
#     spc
# end

# function ViscoElasticSolidSpecific(
#     bulk_modulus::Array{Float64,1},
#     shear_modulus::Array{Float64,1},
#     critical_stretch::Array{Float64,1},
#     density::Array{Float64,1},
#     σγ;
#     criteria = VonMises())
#     return ViscoElasticSolidSpecific(
#         make_matrix(bulk_modulus),
#         make_matrix(shear_modulus),
#         make_matrix(critical_stretch),
#         density,
#         make_matrix(σγ),
#         make_matrix(σγ),
#         make_matrix(zero(σγ)),
#         [],
#         [],
#         criteria,
#     )
# end

# struct ViscoElasticSolidMaterial <: PeridynamicsMaterial
#     @PeridynamicsMaterial_gf
#     specific::ViscoElasticSolidSpecific
# end


# function PeridynamicsMaterial(name, type, bid, gen, spc::ViscoElasticSolidSpecific)
#     ViscoElasticSolidMaterial(name, type, bid, gen, init(spc, gen))
# end


# """
#     force_density_T!(mat::BondBasedMaterial)

# Calculates force density (actually acceleration) for bond based material type.
# """
# function force_density_T!(f, y::AbstractArray{Float64,2}, limits, mat::ViscoElasticSolidMaterial; kwargs...)
#     force_density_T!(f, y, limits, mat, :cpu; kwargs)
# end

# function out_of_loop!(y_, mat::ViscoElasticSolidMaterial, device)
#     @timeit timings "dilatation" dilatation!(y_, mat, device) # calculate dilatation
#     fill!(mat.specific.td_norm2, 0)
#     @timeit timings "td_norm2" loop_over_neighs!(y_, mat, td_norm2_ij!, device) # calculate td_norm2
# end

# function get_items(mat::ViscoElasticSolidMaterial)
#     return mat.general.weighted_volume,
#                     mat.specific.theta,
#                     mat.specific.td_norm2,
#                     mat.specific.bulk_modulus,
#                     mat.specific.shear_modulus,
#                     mat.specific.σγ,
#                     mat.specific.ψ,
#                     mat.specific.edp,
#                     mat.general.horizon,
#                     mat.general.volume,
#                     mat.specific.criteria
# end

# function force_density_t_ij(mat::ViscoElasticSolidMaterial,
#         i, j, k, type1, type2,
#         xij, yij, extension, s,
#         wij, wji, items)
#     m, theta, td_norm2, K, G, σγ, ψ, edp, δ, vol, criteria = items
#     Kij = K[type1, type2]
#     Gij = G[type1, type2]
#     ψij = ψ[type1, type2]
#     a, b = wij/m[i], wji/m[j]
#     θω_m = (theta[i] * a + theta[j] * b)
#     theta_ij = θω_m / (a + b)
#     ed = extension - theta_ij * xij / 3
#     ede = ed - edp[k, i]
#     td_trial = 15 * Gij * ede * (a + b)
#     pyf = td_norm2[i] - 2ψij
#     if pyf < 0
#         td = td_trial
#     else
#         td = sqrt(2 * ψij) * td_trial / sqrt(td_norm2[i])
#         td_extra = td_trial - td
#         plastic_ed = td_extra / (15 * Gij * (a + b))
#         edp[k, i] = edp[k, i] + plastic_ed
#     end
#     tn = (3 * Kij) * xij * θω_m
#     t = tn + td
#     return t
# end

# function td_norm2_ij!(mat, i, j, k, type1, type2,
#                         xij, yij, extension, s,
#                         wij, wji, items)
#     m, theta, td_norm2, K, G, σγ, ψ, edp, δ, vol, criteria = items
#     Gij = G[type1, type2]
#     edpij = edp[k, i]
#     a, b = wij/m[i], wji/m[j]
#     ω_m = (a + b)
#     θω_m = (theta[i] * a + theta[j] * b)
#     theta_ij = θω_m / ω_m
#     ed = extension - theta_ij * xij / 3
#     td = 15Gij * (ed - edpij) * ω_m
#     td_norm2[i] += td^2 * vol[j]
# end






# # """
# #     force_density_T!(y::Array{Float64,2},mat::ViscoElasticSolidMaterial)

# # Calculates force density (actually acceleration) for ordinary state based material type.
# # """
# # function force_density_T!(y::Array{Float64,2}, mat::ViscoElasticSolidMaterial)
# #     types = mat.general.type
# #     force = zeros(size(y))
# #     X = Float64[1.0, 0.0, 0.0]
# #     Y = Float64[1.0, 0.0, 0.0]
# #     S = mat.general
# #     m = mat.general.weighted_volume
# #     theta = dilatation(y, S, m)
# #     td_norm2 = td_norm2_fn(y, theta, mat)
# #     δ = S.horizon
# #     edp = mat.specific.edp
# #     N = size(S.x, 2)
# #     _max = 0.0
# #     for i = 1:N
# #         for k = 1:size(S.family, 1)
# #             j = S.family[k, i]::Int64
# #             if !S.intact[k, i] || i == j
# #                 # force[1,i] += 0.0
# #                 # force[2,i] += 0.0
# #                 # force[3,i] += 0.0
# #                 S.intact[k, i] = false
# #             else
# #                 type1 = types[i] - mat.type.start + 1
# #                 type2 = types[j] - mat.type.start + 1
# #                 K = mat.specific.bulk_modulus[type1, type2]
# #                 G = mat.specific.shear_modulus[type1, type2]
# #                 σ = mat.specific.σγ[type1, type2]
# #                 X[1], X[2], X[3] =
# #                     S.x[1, j] - S.x[1, i], S.x[2, j] - S.x[2, i], S.x[3, j] - S.x[3, i]
# #                 Y[1], Y[2], Y[3] = y[1, j] - y[1, i], y[2, j] - y[2, i], y[3, j] - y[3, i]
# #                 _Y, _X = @_magnitude(Y), @_magnitude(X)
# #                 xij = @_magnitude(X)::Float64
# #                 wij = influence_function(X)::Float64
# #                 wji = influence_function(-X)::Float64

# #                 e = _Y - _X
# #                 ed = e - theta[i] * _X / 3
# #                 ee = ed - edp[k, i]

# #                 td_trial = 15 * G * (ee * wij / m[i] + ee * wji / m[j])

# #                 pyf, ψ = plastic_yield_function(td_norm2[i], σ, δ)

# #                 if pyf < 0
# #                     td = td_trial
# #                 else
# #                     td = sqrt(2 * ψ) * td_trial / td_norm2[i]
# #                     td_extra = td_trial - td
# #                     plastic_ed = td_extra / (15 * G * (wij / m[i] + wji / m[j]))
# #                     edp[k, i] = edp[k, i] + plastic_ed
# #                 end
# #                 if (e / xij) < mat.specific.critical_stretch[type1, type2]
# #                     t =
# #                         (3 * K) *
# #                         (theta[i] * xij * wij / m[i] + theta[j] * xij * wji / m[j])
# #                     t = t + td
# #                     t = t * S.volume[i] / @_magnitude(Y)
# #                     force[1, i] += t * Y[1]
# #                     force[2, i] += t * Y[2]
# #                     force[3, i] += t * Y[3]
# #                 else
# #                     # force[1,i] += 0.0
# #                     # force[2,i] += 0.0
# #                     # force[3,i] += 0.0
# #                     S.intact[k, i] = false
# #                 end
# #             end
# #         end
# #     end
# #     return force
# # end


# # """
# #     _tij(x,mi,thetai,eij,K,G)::Float64

# # calculates tij force density magnitude (actually acceleration).
# # """
# # function _tij(x, mi, thetai, eij, K, G)::Float64
# #     xij = @_magnitude(x)::Float64
# #     wij = influence_function(x)::Float64
# #     return ((3 * K - 5 * G) * (thetai * xij) + 15 * G * eij) * wij / mi
# # end


# # """
# #     td_norm2_fn
# # """
# # function td_norm2_fn(y, theta, mat)
# #     S = mat.general
# #     types = mat.general.type
# #     intact = S.intact
# #     family = S.family
# #     x = S.x
# #     N = size(family, 2)
# #     M = size(family, 1)
# #     ARGS = map((i) -> (i, 1:M), 1:N)
# #     m = mat.general.weighted_volume
# #     function with_if_cal_td_2(i, k)
# #         if intact[k, i]
# #             j = family[k, i]
# #             type1 = types[i] - mat.type.start + 1
# #             type2 = types[j] - mat.type.start + 1
# #             G = mat.specific.shear_modulus[type1, type2]
# #             edp = mat.specific.edp[k, i]

# #             X = [x[1, j] - x[1, i], x[2, j] - x[2, i], x[3, j] - x[3, i]]
# #             Y = [y[1, j] - y[1, i], y[2, j] - y[2, i], y[3, j] - y[3, i]]

# #             _X = @_magnitude(X)
# #             _Y = @_magnitude(Y)

# #             e = _Y - _X

# #             ed = e - theta[i] * _X / 3

# #             ee = ed - edp

# #             s = e / _X

# #             wij = influence_function(X)::Float64
# #             wji = influence_function(-X)::Float64
# #             type1 = types[i] - mat.type.start + 1
# #             type2 = types[j] - mat.type.start + 1
# #             if s < mat.specific.critical_stretch[type1, type2]
# #                 K = mat.specific.bulk_modulus[type1, type2]
# #                 G = mat.specific.shear_modulus[type1, type2]
# #                 td = 15 * G * (ee * wij / m[i] + ee * wji / m[j])
# #                 return td^2
# #             else
# #                 intact[k, i] = false
# #                 return 0.0
# #             end
# #         else
# #             return 0.0
# #         end
# #     end

# #     inner_map(i, inds) = map_reduce((j)-> with_if_cal_force_ij(i,j), +, inds)

# #     # inner_map(i, inds) = mapreduce((j)-> with_if_cal_force_ij(i,j), +, inds)
# #     outer_map(ARGS) = map((x) -> inner_map(x[1], x[2]), ARGS)
# #     return outer_map(ARGS)
# # end

