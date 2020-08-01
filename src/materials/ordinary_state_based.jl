export OrdinaryStateBasedMaterial, OrdinaryStateBasedSpecific, force_density_T

"""
Specific ordinary state based materail type.
"""
struct OrdinaryStateBasedSpecific
    bulk_modulus::Float64
    shear_modulus::Float64
    weighted_volume::Array{Float64,1}
end

"""
    OrdinaryStateBasedSpecific(bulk_modulus::Float64, shear_modulus::Float64, mat_gen::GeneralMaterial)
Specific ordinary state based materail type.
"""
function OrdinaryStateBasedSpecific(bulk_modulus::Float64, shear_modulus::Float64, mat_gen::GeneralMaterial)
    m = weighted_volume(mat_gen)
    return OrdinaryStateBasedSpecific(bulk_modulus, shear_modulus, m)
end


struct OrdinaryStateBasedMaterial<:PeridynamicsMaterial
    type::Int64
    general::GeneralMaterial
    specific::OrdinaryStateBasedSpecific
end

"""
    force_density_T(y::Array{Float64,2},mat::OrdinaryStateBasedMaterial)

Calculates force density (actually acceleration) for ordinary state based material type.
"""
function force_density_T(y::Array{Float64,2},mat::OrdinaryStateBasedMaterial)
    force = zeros(size(y))
    X = Float64[1.0,0.0,0.0]
    Y = Float64[1.0,0.0,0.0]
    S = mat.general
    K = mat.specific.bulk_modulus
    G = mat.specific.shear_modulus
    m = mat.specific.weighted_volume
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
                e = (magnitude(Y) - magnitude(X))
                xij = magnitude(X)::Float64
                wij = influence_function(X)::Float64
                wji = influence_function(-X)::Float64
                    if (e/xij)<S.critical_stretch
                        t = ( (((3*K-5*G)*(theta[i]*xij*wij/m[i]+theta[j]*xij*wji/m[j]) + 15*G*(e*wji/m[i]+e*wji/m[j]))) )*S.volume[j]/magnitude(Y)

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
    xij = magnitude(x)::Float64
    wij = influence_function(x)::Float64
    return ((3*K-5*G)*(thetai*xij) + 15*G*eij)*wij/mi
end



#
