export ElastoPlasticSolidMaterial,
    ElastoPlasticSolidSpecific,
    ElastoPlasticSolidSpecificFull,
    force_density_T,
    PeridynamicsMaterial



abstract type PlasticFailureCriteria end

struct VonMises <: PlasticFailureCriteria
end
    

"""
Specific Elasto Plastic Solid Material type.
"""
struct ElastoPlasticSolidSpecific <: SpecificMaterial
    bulk_modulus::Array{Float64,2}
    shear_modulus::Array{Float64,2}
    critical_stretch::Array{Float64,2}
    density::Array{Float64,1}
    σy::Array{Float64,2}
    criteria::PlasticFailureCriteria
end

"""
Specific Elasto Plastic Solid Material type.
"""
struct ElastoPlasticSolidSpecificFull <: SpecificMaterial
    bulk_modulus::Array{Float64,2}
    shear_modulus::Array{Float64,2}
    critical_stretch::Array{Float64,2}
    density::Array{Float64,1}
    σy::Array{Float64,2}
    edp::Array{Float64,2}
    criteria::PlasticFailureCriteria
end

"""
    ElastoPlasticSolidSpecific(bulk_modulus::Array{Float64, 1}, shear_modulus::Array{Float64,1}, critical_stretch::Array{Float64,1}, density::Array{Float64,1})
Specific Elasto Plastic Solid Material type.
"""
function ElastoPlasticSolidSpecific(
    bulk_modulus::Array{Float64,1},
    shear_modulus::Array{Float64,1},
    critical_stretch::Array{Float64,1},
    density::Array{Float64,1},
    σy;
    criteria = VonMises()
)
    return ElastoPlasticSolidSpecific(
        make_matrix(bulk_modulus),
        make_matrix(shear_modulus),
        make_matrix(critical_stretch),
        density,
        make_matrix(σy),
        criteria,
    )
end

"""
    ElastoPlasticSolidMaterial(type::UnitRange{Int64}, general::GeneralMaterial, specific::ElastoPlasticSolidSpecific)
"""
struct ElastoPlasticSolidMaterial <: PeridynamicsMaterial
    name::String
    type::UnitRange{Int64}
    general::GeneralMaterial
    specific::ElastoPlasticSolidSpecificFull
end

"""
    PeridynamicsMaterial(gen, spc::ElastoPlasticSolidSpecific)
"""
function PeridynamicsMaterial(name, type, gen, spc::ElastoPlasticSolidSpecific)
    edp = 0 * gen.family
    _spc = ElastoPlasticSolidSpecificFull(
        spc.bulk_modulus,
        spc.shear_modulus,
        spc.critical_stretch,
        spc.density,
        spc.σy,
        edp,
        spc.criteria,
    )
    ElastoPlasticSolidMaterial(name, type, gen, _spc)
end


"""
    force_density_T(y::Array{Float64,2},mat::ElastoPlasticSolidMaterial)

Calculates force density (actually acceleration) for ordinary state based material type.
"""
function force_density_T(y::Array{Float64,2}, mat::ElastoPlasticSolidMaterial)
    types = mat.general.type
    force = zeros(size(y))
    X = Float64[1.0, 0.0, 0.0]
    Y = Float64[1.0, 0.0, 0.0]
    S = mat.general
    m = mat.general.weighted_volume
    theta = dilatation(y, S, m)
    td_norm = td_norm_fn(y, theta, mat)
    δ = S.horizon
    edp = mat.specific.edp
    N = size(S.x, 2)
    _max = 0.0
    for i = 1:N
        for k = 1:size(S.family, 1)
            j = S.family[k, i]::Int64
            if !S.intact[k, i] || i == j
                # force[1,i] += 0.0
                # force[2,i] += 0.0
                # force[3,i] += 0.0
                S.intact[k, i] = false
            else
                type1 = types[i] - mat.type.start + 1
                type2 = types[j] - mat.type.start + 1
                K = mat.specific.bulk_modulus[type1, type2]
                G = mat.specific.shear_modulus[type1, type2]
                σ = mat.specific.σy[type1, type2]
                X[1], X[2], X[3] =
                    S.x[1, j] - S.x[1, i], S.x[2, j] - S.x[2, i], S.x[3, j] - S.x[3, i]
                Y[1], Y[2], Y[3] = y[1, j] - y[1, i], y[2, j] - y[2, i], y[3, j] - y[3, i]
                _Y, _X = @_magnitude(Y), @_magnitude(X)
                xij = @_magnitude(X)::Float64
                wij = influence_function(X)::Float64
                wji = influence_function(-X)::Float64

                e = _Y - _X
                ed = e - theta[i] * _X / 3
                ee = ed - edp[k, i]

                td_trial = 15 * G * (ee * wij / m[i] + ee * wji / m[j])

                pyf, ψ = plastic_yield_function(td_norm[i], σ, δ)

                if pyf < 0
                    td = td_trial
                else
                    td = sqrt(2 * ψ) * td_trial / td_norm[i]
                    td_extra = td_trial - td
                    plastic_ed = td_extra / (15 * G * (wij / m[i] + wji / m[j]))
                    edp[k, i] = edp[k, i] + plastic_ed
                end
                if (e / xij) < mat.specific.critical_stretch[type1, type2]
                    t =
                        (3 * K) *
                        (theta[i] * xij * wij / m[i] + theta[j] * xij * wji / m[j])
                    t = t + td
                    t = t * S.volume[i] / @_magnitude(Y)
                    force[1, i] += t * Y[1]
                    force[2, i] += t * Y[2]
                    force[3, i] += t * Y[3]
                else
                    # force[1,i] += 0.0
                    # force[2,i] += 0.0
                    # force[3,i] += 0.0
                    S.intact[k, i] = false
                end
            end
        end
    end
    return force
end


"""
    _tij(x,mi,thetai,eij,K,G)::Float64

calculates tij force density magnitude (actually acceleration).
"""
function _tij(x, mi, thetai, eij, K, G)::Float64
    xij = @_magnitude(x)::Float64
    wij = influence_function(x)::Float64
    return ((3 * K - 5 * G) * (thetai * xij) + 15 * G * eij) * wij / mi
end


"""
    td_norm_fn
"""
function td_norm_fn(y, theta, mat)
    S = mat.general
    types = mat.general.type
    intact = S.intact
    family = S.family
    x = S.x
    N = size(family, 2)
    M = size(family, 1)
    ARGS = map((i) -> (i, 1:M), 1:N)
    m = mat.general.weighted_volume
    function with_if_cal_td_2(i, k)
        if intact[k, i]
            j = family[k, i]
            type1 = types[i] - mat.type.start + 1
            type2 = types[j] - mat.type.start + 1
            G = mat.specific.shear_modulus[type1, type2]
            edp = mat.specific.edp[k, i]

            X = [x[1, j] - x[1, i], x[2, j] - x[2, i], x[3, j] - x[3, i]]
            Y = [y[1, j] - y[1, i], y[2, j] - y[2, i], y[3, j] - y[3, i]]

            _X = @_magnitude(X)
            _Y = @_magnitude(Y)

            e = _Y - _X

            ed = e - theta[i] * _X / 3

            ee = ed - edp

            s = e / _X

            wij = influence_function(X)::Float64
            wji = influence_function(-X)::Float64
            type1 = types[i] - mat.type.start + 1
            type2 = types[j] - mat.type.start + 1
            if s < mat.specific.critical_stretch[type1, type2]
                K = mat.specific.bulk_modulus[type1, type2]
                G = mat.specific.shear_modulus[type1, type2]
                td = 15 * G * (ee * wij / m[i] + ee * wji / m[j])
                return td^2
            else
                intact[k, i] = false
                return 0.0
            end
        else
            return 0.0
        end
    end
    
    inner_map(i, inds) = map_reduce((j)-> with_if_cal_force_ij(i,j), +, inds)
    
    # inner_map(i, inds) = mapreduce((j)-> with_if_cal_force_ij(i,j), +, inds)
    outer_map(ARGS) = map((x) -> inner_map(x[1], x[2]), ARGS)
    return outer_map(ARGS)
end



"""
    yield function
"""
function plastic_yield_function(td_norm, σ, δ)
    Ey = VonMises_effective_stress(σ, 0, 0)
    ψ = 75 * Ey^2 / (8 * pi * δ^5)
    pyf = td_norm^2 / 2 - ψ
    return pyf, ψ
end

function VonMises_effective_stress(σ₁, σ₂, σ₃)
    return sqrt(((σ₁ - σ₂)^2 + (σ₂ - σ₃)^2 + (σ₁ - σ₃)^2) / 6)
end

#
