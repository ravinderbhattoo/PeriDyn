struct OrdinaryStateBasedSpecific
    bulk_modulus::Float64
    shear_modulus::Float64
    weighted_volume::Array{Float64,1}
end

function OrdinaryStateBasedSpecific(bulk_modulus, shear_modulus, mat_gen)
    m = weighted_volume(mat_gen)
    return OrdinaryStateBasedSpecific(bulk_modulus, shear_modulus,m)
end


struct OrdinaryStateBasedMaterial<:PeridynamicsMaterial
    type::Int64
    general::AbstractMaterial
    specific::OrdinaryStateBasedSpecific
end


function s_force_density_T(y::Array{Float64,2},mat::OrdinaryStateBasedMaterial)
    S = mat.general
    K = mat.specific.bulk_modulus
    G = mat.specific.shear_modulus
    m = mat.specific.weighted_volume
    theta = dilatation(y,S,m)
    force = y*0
    for i in 1:size(S.x,1)
        for k in 1:size(S.family,2)
            if !S.damage[i,k]
                force[i,:] .+= [0,0,0]
            else
                j = S.family[i,k]
                E = [i,j]
                E_ = [j,i]
                Y = s_Y(E,y)
                X = s_X(E,S.x)
                e = (s_magnitude(Y) - s_magnitude(X))
                    if (e/s_magnitude(X))<S.critical_stretch
                        M = Y./s_magnitude(Y)
                        t = _tij(E,S.x,m,theta,e,K,G)
                        t_ = _tij(E_,S.x,m,theta,e,K,G)
                        force[i,:] += ((t + t_).*M)*S.volume[j]
                    else
                        force[i,:] .+= [0,0,0]
                        S.damage[i,k] = true
                    end
            end
        end
    end
    return force
end

function _tij(E,x,m,theta,eij,K,G)
    xij = s_magnitude(s_X(E,x))
    wij = influence_function(s_X(E,x))
    mi = m[E[1]]
    thetai = theta[E[1]]
    return (3*K-5*G)*(thetai*wij*xij)/mi + 15*G*wij*eij/mi
end



#
