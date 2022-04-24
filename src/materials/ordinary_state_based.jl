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


# """
#     force_density_T(y::Array{Float64,2},mat::OrdinaryStateBasedMaterial)

# Calculates force density (actually acceleration) for ordinary state based material type.
# """
# function PeriDyn.force_density_T(y::Array{Float64,2}, mat::OrdinaryStateBasedMaterial; particles=nothing)
    
#     force = Zygote.Buffer(y)
        
#     types = mat.general.type
#     gen = mat.general
#     m = gen.weighted_volume
#     intact = gen.intact
    
#     intact = Zygote.Buffer(intact)
    
#     family = gen.family
#     x = gen.x
#     volume = gen.volume
    
#     if isnothing(particles) 
#         _N = 1:size(family, 2)
#     else
#         _N = particles
#     end
    
#     theta = dilatation(y, gen, m)
    
#     N = size(x, 2)
#     for i in _N
#         force[1,i] = 0.0
#         force[2,i] = 0.0
#         force[3,i] = 0.0
#         for k in 1:size(family,1)
#             j = family[k,i]
#             if !intact[k,i]
#                 intact[k,i] = false
#             else
#                 X = get_ij(j,i,x)
#                 Y = get_ij(j,i,y)
#                 yij = get_magnitude(Y)
#                 xij = get_magnitude(X)
                
#                 extention = yij - xij
#                 wij = influence_function(X)
#                 wji = influence_function(-X)
#                 type1 = types[i] - mat.type.start + 1
#                 type2 = types[j] - mat.type.start + 1
#                 if (extention/xij) < mat.specific.critical_stretch[type1, type2]
#                     K = mat.specific.bulk_modulus[type1, type2]
#                     G = mat.specific.shear_modulus[type1, type2]
#                     a_, b_ = wij/m[i], wji/m[j]
#                     t =  (3*K-5*G) * (theta[i]*xij*a_ + theta[j]*xij*b_)  + 15*G*extention*(a_ + b_) 
#                     t = t*volume[j] / yij
#                     force[1,i] += t*Y[1]
#                     force[2,i] += t*Y[2]
#                     force[3,i] += t*Y[3]
#                 else
#                     intact[k,i] = false
#                 end
#             end
#         end
#     end
#     return copy(force)
# end



"""
    force_density_T(y::Array{Float64,2},mat::OrdinaryStateBasedMaterial)

Calculates force density (actually acceleration) for ordinary state based material type.
"""
function PeriDyn.force_density_T(y::Array{Float64,2}, mat::OrdinaryStateBasedMaterial; particles=nothing)
            
    types = mat.general.type
    gen = mat.general
    m = gen.weighted_volume
    intact = gen.intact    
    family = gen.family
    x = gen.x
    volume = gen.volume
    
    if isnothing(particles) 
        _N = 1:size(family, 2)
    else
        _N = particles
    end
    
    theta = dilatation(y, gen, m)
    
    M = size(family, 1)
    ARGS = map((i) -> (i, 1:M), _N)    
    
    function with_if_cal_force_ij(i, k)
        if intact[k,i]
            j = family[k,i]            
            X = get_ij(j,i,x)
            Y = get_ij(j,i,y)
            yij = get_magnitude(Y)
            xij = get_magnitude(X)
            
            extention = yij - xij
            s = extention/xij
            
            wij = influence_function(X)
            wji = influence_function(-X)
            type1 = types[i] - mat.type.start + 1
            type2 = types[j] - mat.type.start + 1
            if s < mat.specific.critical_stretch[type1, type2]
                K = mat.specific.bulk_modulus[type1, type2]
                G = mat.specific.shear_modulus[type1, type2]
                a_, b_ = wij/m[i], wji/m[j]
                t =  (3*K-5*G) * (theta[i]*xij*a_ + theta[j]*xij*b_)  + 15*G*extention*(a_ + b_) 
                t = t*volume[j] / yij
                return t*Y
            else
                intact[k,i] = false
                return [0.0, 0.0, 0.0]
            end
        else
            return [0.0, 0.0, 0.0]
        end
    end
    
    inner_map(i, inds) = map_reduce((j)-> with_if_cal_force_ij(i,j), +, inds)
    
    # inner_map(i, inds) = mapreduce((j)-> with_if_cal_force_ij(i,j), +, inds)
    outer_map(ARGS) = map((x)->inner_map(x[1], x[2]), ARGS)
    
    return hcat(outer_map(ARGS)...)
    
 
end



#
