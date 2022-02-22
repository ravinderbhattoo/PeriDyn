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
    name::String
    type::UnitRange{Int64}
    general::GeneralMaterial
    specific::OrdinaryStateBasedSpecific
end

function PeridynamicsMaterial(name, type, gen, spc::OrdinaryStateBasedSpecific)
    OrdinaryStateBasedMaterial(name, type, gen, spc) 
end


"""
    force_density_T(y::Array{Float64,2},mat::OrdinaryStateBasedMaterial)

Calculates force density (actually acceleration) for ordinary state based material type.
"""
function PeriDyn.force_density_T(y::Array{Float64,2}, mat::OrdinaryStateBasedMaterial; particles=nothing)
    force = zeros(size(y))
    types = mat.general.type
    S = mat.general
    m = mat.general.weighted_volume
    intact = S.intact
    family = S.family

    if isnothing(particles) 
        _N = 1:size(family, 2)
    else
        _N = particles
    end

    theta = dilatation(y, S, m)
    N = size(S.x, 2)
    Threads.@threads for i in _N
        for k in 1:size(S.family,1)
            j = S.family[k,i]
            if !S.intact[k,i]
                S.intact[k,i] = false
            else
                X = get_ij(j,i,S.x)
                Y = get_ij(j,i,y)
                yij = get_magnitude(Y)
                xij = get_magnitude(X)
                
                extention = yij - xij
                wij = influence_function(X)
                wji = influence_function(-X)
                type1 = types[i] - mat.type.start + 1
                type2 = types[j] - mat.type.start + 1
                if (extention/xij) < mat.specific.critical_stretch[type1, type2]
                    K = mat.specific.bulk_modulus[type1, type2]
                    G = mat.specific.shear_modulus[type1, type2]
                    a_, b_ = wij/m[i], wji/m[j]
                    t =  (3*K-5*G) * (theta[i]*xij*a_ + theta[j]*xij*b_)  + 15*G*extention*(a_ + b_) 
                    t = t*S.volume[j] / yij
                    force[1,i] += t*Y[1]
                    force[2,i] += t*Y[2]
                    force[3,i] += t*Y[3]
                else
                    S.intact[k,i] = false
                end
            end
        end
    end
    return force
end



#
