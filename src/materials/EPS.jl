export ElastoPlasticSolidMaterial, ElastoPlasticSolidSpecific, ElastoPlasticSolidSpecificFull, force_density_T, PeridynamicsMaterial

"""
Specific Elasto Plastic Solid Material type.
"""
struct ElastoPlasticSolidSpecific <: SpecificMaterial
    bulk_modulus::Array{Float64, 2}
    shear_modulus::Array{Float64, 2}
    critical_stretch::Array{Float64, 2}
    density::Array{Float64, 1}
    σy::Array{Float64, 2}
end

"""
Specific Elasto Plastic Solid Material type.
"""
struct ElastoPlasticSolidSpecificFull <: SpecificMaterial
    bulk_modulus::Array{Float64, 2}
    shear_modulus::Array{Float64, 2}
    critical_stretch::Array{Float64, 2}
    density::Array{Float64, 1}
    σy::Array{Float64, 2}
    edp::Array{Float64, 2}
end

"""
    ElastoPlasticSolidSpecific(bulk_modulus::Array{Float64, 1}, shear_modulus::Array{Float64,1}, critical_stretch::Array{Float64,1}, density::Array{Float64,1})
Specific Elasto Plastic Solid Material type.
"""
function ElastoPlasticSolidSpecific(bulk_modulus::Array{Float64, 1}, shear_modulus::Array{Float64,1}, critical_stretch::Array{Float64,1}, density::Array{Float64,1}, σy)
    return ElastoPlasticSolidSpecific(make_matrix(bulk_modulus), make_matrix(shear_modulus), make_matrix(critical_stretch), density, make_matrix(σy))
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
function PeridynamicsMaterial(gen, spc::ElastoPlasticSolidSpecific; name="Default")
    type = minimum(gen.type):maximum(gen.type)
    edp = 0*gen.family
    _spc = ElastoPlasticSolidSpecificFull(spc.bulk_modulus, spc.shear_modulus, spc.critical_stretch, spc.density, spc.σy, edp)
    ElastoPlasticSolidMaterial(name, type, gen, _spc)
end


"""
    force_density_T(y::Array{Float64,2},mat::ElastoPlasticSolidMaterial)

Calculates force density (actually acceleration) for ordinary state based material type.
"""
function force_density_T(y::Array{Float64,2}, mat::ElastoPlasticSolidMaterial)
    types = mat.general.type
    force = zeros(size(y))
    X = Float64[1.0,0.0,0.0]
    Y = Float64[1.0,0.0,0.0]
    S = mat.general
    m = mat.general.weighted_volume
    theta = dilatation(y, S, m)
    td_trial = trial_force_density_deviatoric(y, theta, mat)
    δ = S.horizon
    N = size(S.x, 2)
    for i in 1:N
        for k in 1:size(S.family,1)
            j = S.family[k,i]::Int64
            if !S.intact[k,i] || i==j
                # force[1,i] += 0.0
                # force[2,i] += 0.0
                # force[3,i] += 0.0
                S.intact[k,i] = false
            else
                type1 = types[i] - mat.type.start + 1
                type2 = types[j] - mat.type.start + 1
                K = mat.specific.bulk_modulus[type1, type2]
                G = mat.specific.shear_modulus[type1, type2]
                σ = mat.specific.σy[type1, type2]
                pyf, ψ = plastic_yield_function(td_trial[i], σ, δ)
                if pyf < 0
                    td = td_trial[i]
                else
                    α = 15*G/m[i]
                    Δλ = 1/α*(mod(td_trial[i])/sqrt(2*ψ) -1)
                    td = sqrt(2*ψ) * td_trial[i] / mod(td_trial[i])
                    edp[i] = edp[i] + Δλ*(td[i])
                end
                X[1],X[2],X[3] = S.x[1,j]-S.x[1,i],S.x[2,j]-S.x[2,i],S.x[3,j]-S.x[3,i]
                Y[1],Y[2],Y[3] = y[1,j]-y[1,i],y[2,j]-y[2,i],y[3,j]-y[3,i]
                e = (_magnitude(Y) - _magnitude(X))
                xij = _magnitude(X)::Float64
                wij = influence_function(X)::Float64
                wji = influence_function(-X)::Float64
                if (e/xij) < mat.specific.critical_stretch[type1, type2]

                    t = (3*K)*(theta[i]*xij*wij/m[i]+theta[j]*xij*wji/m[j])  

                    t = t + td

                    t = t*S.volume[i]/_magnitude(Y)

                    force[1,i] += t*Y[1]
                    force[2,i] += t*Y[2]
                    force[3,i] += t*Y[3]
                else
                    # force[1,i] += 0.0
                    # force[2,i] += 0.0
                    # force[3,i] += 0.0
                    S.intact[k,i] = false
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
function _tij(x,mi,thetai,eij,K,G)::Float64
    xij = _magnitude(x)::Float64
    wij = influence_function(x)::Float64
    return ((3*K-5*G)*(thetai*xij) + 15*G*eij)*wij/mi
end


"""
    trial_force_density_deviatoric
"""
function trial_force_density_deviatoric(y, theta, mat)
    S = mat.general
    types = mat.general.type
    intact = S.intact
    family = S.family
    x = S.x
    N = size(family, 2)
    M = size(family, 1)
    ARGS = map((i) -> (i, 1:M), 1:N)    
    m = mat.general.weighted_volume
    function with_if_cal_td_norm(i, k)
        if intact[k,i]
            j = family[k, i]
            type1 = types[i] - mat.type.start + 1
            type2 = types[j] - mat.type.start + 1
            G = mat.specific.shear_modulus[type1, type2]
            edp = mat.specific.edp[k, i]

            X = [x[1,j]-x[1,i],x[2,j]-x[2,i],x[3,j]-x[3,i]]
            Y = [y[1,j]-y[1,i],y[2,j]-y[2,i],y[3,j]-y[3,i]]

            _X = _magnitude(X)
            _Y = _magnitude(Y)

            e = _Y - _X

            ed = e - theta[i] * _X / 3
            
            ee = ed - edp

            s = e/_X

            wij = influence_function(X)::Float64
            wji = influence_function(-X)::Float64
            type1 = types[i] - mat.type.start + 1
            type2 = types[j] - mat.type.start + 1
            if s < mat.specific.critical_stretch[type1, type2]
                K = mat.specific.bulk_modulus[type1, type2]
                G = mat.specific.shear_modulus[type1, type2]
                t = 15*G*(ee*wji/m[i]+ee*wji/m[j])
                t = t*S.volume[i]
                return t
            else
                intact[k,i] = false
                return 0.0
            end            
        else
            return 0.0
        end
    end

    inner_map(i, inds) = sum(map((j)-> with_if_cal_td_norm(i,j), inds))
    outer_map(ARGS) = map((x)->inner_map(x[1], x[2]), ARGS)
    return outer_map(ARGS)
end


"""
    yield function
"""
function plastic_yield_function(td, σ, δ)
    ψ = 25σ^2 / (8 * pi * δ^5)
    return mod(td)^2/2 - ψ, ψ 
end
#
