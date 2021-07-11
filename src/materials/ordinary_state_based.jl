export OrdinaryStateBasedMaterial, OrdinaryStateBasedSpecific, force_density_T, PeridynamicsMaterial

"""
Specific ordinary state based materail type.
"""
struct OrdinaryStateBasedSpecific <: SpecificMaterial
    bulk_modulus::Array{Float64, 2}
    shear_modulus::Array{Float64, 2}
    critical_stretch::Array{Float64, 2}
    density::Array{Float64, 1}
end

"""
    OrdinaryStateBasedSpecific(bulk_modulus::Array{Float64, 1}, shear_modulus::Array{Float64,1}, critical_stretch::Array{Float64,1}, density::Array{Float64,1})
Specific ordinary state based materail type.
"""
function OrdinaryStateBasedSpecific(bulk_modulus::Array{Float64, 1}, shear_modulus::Array{Float64,1}, critical_stretch::Array{Float64,1}, density::Array{Float64,1})
    return OrdinaryStateBasedSpecific(make_matrix(bulk_modulus), make_matrix(shear_modulus), make_matrix(critical_stretch), density)
end

"""
    OrdinaryStateBasedMaterial(type::UnitRange{Int64}, general::GeneralMaterial, specific::OrdinaryStateBasedSpecific)
"""
struct OrdinaryStateBasedMaterial <: PeridynamicsMaterial
    type::UnitRange{Int64}
    general::GeneralMaterial
    specific::OrdinaryStateBasedSpecific
end

"""
    PeridynamicsMaterial(gen, spc::OrdinaryStateBasedSpecific)
"""
function PeridynamicsMaterial(gen, spc::OrdinaryStateBasedSpecific)
    type = minimum(gen.type):maximum(gen.type)
    OrdinaryStateBasedMaterial(type, gen, spc)
end


"""
    force_density_T(y::Array{Float64,2},mat::OrdinaryStateBasedMaterial)

Calculates force density (actually acceleration) for ordinary state based material type.
"""
function force_density_T(y::Array{Float64,2}, mat::OrdinaryStateBasedMaterial)
    types = mat.general.type
    force = zeros(size(y))
    X = Float64[1.0,0.0,0.0]
    Y = Float64[1.0,0.0,0.0]
    S = mat.general
    m = mat.general.weighted_volume
    theta = dilatation(y,S,m)
    for i in 1:size(S.x,2)
        for k in 1:size(S.family,1)
            if !S.intact[k,i]
                force[1,i] += 0.0
                force[2,i] += 0.0
                force[3,i] += 0.0
            else
                j = S.family[k,i]::Int64
                X[1],X[2],X[3] = S.x[1,j]-S.x[1,i],S.x[2,j]-S.x[2,i],S.x[3,j]-S.x[3,i]
                Y[1],Y[2],Y[3] = y[1,j]-y[1,i],y[2,j]-y[2,i],y[3,j]-y[3,i]
                e = (_magnitude(Y) - _magnitude(X))
                xij = _magnitude(X)::Float64
                wij = influence_function(X)::Float64
                wji = influence_function(-X)::Float64
                type1 = types[i] - mat.type.start + 1
                type2 = types[j] - mat.type.start + 1
                    if (e/xij)<mat.specific.critical_stretch[type1, type2]
                        K = mat.specific.bulk_modulus[type1, type2]
                        G = mat.specific.shear_modulus[type1, type2]
                        t = ( (((3*K-5*G)*(theta[i]*xij*wij/m[i]+theta[j]*xij*wji/m[j]) + 15*G*(e*wji/m[i]+e*wji/m[j]))) )*S.volume[j]/_magnitude(Y)
                        force[1,i] += t*Y[1]
                        force[2,i] += t*Y[2]
                        force[3,i] += t*Y[3]
                    else
                        force[1,i] += 0.0
                        force[2,i] += 0.0
                        force[3,i] += 0.0
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



#
