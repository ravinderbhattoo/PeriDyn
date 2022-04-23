"""
This module contains LJ replusive model definitions.
"""

export LJRepulsionModel, LJRepulsionModel11, LJRepulsionModel12, repulsion_force


struct LJRepulsionModel12<:RepulsionModel12
    pair::Vector{UnitRange{Int64}}
    alpha::Float64
    epsilon::Float64
    equi_dist::Float64
    neighs::Array{Int64,2}
    distance::Float64
    max_neighs::Int64
end

struct LJRepulsionModel11<:RepulsionModel11
    type::UnitRange{Int64}
    material::GeneralMaterial
    alpha::Float64
    epsilon::Float64
    neighs::Array{Int64,2}
    distance::Float64
    max_neighs::Int64
end

"""
    LJRepulsionModel(alpha::Float64,epsilon::Float64, mat1::PeridynamicsMaterial,mat2::PeridynamicsMaterial;distanceX=5,max_neighs=50)

LJ repulsive model for 1-2 material blocks.
"""
function LJRepulsionModel(alpha, epsilon, mat1::PeridynamicsMaterial,mat2::PeridynamicsMaterial;distanceX=5,max_neighs=50)
    p_size = (mat1.general.particle_size+mat2.general.particle_size)/2
    neighs = zeros(min(max_neighs,size(mat2.general.x,2)), size(mat1.general.x,2))
    LJRepulsionModel12([mat1.type,mat2.type],alpha,epsilon,p_size,neighs,p_size*distanceX,min(max_neighs,size(mat2.general.x,2)))
end

"""
    LJRepulsionModel(alpha::Float64,epsilon::Float64, mat1::PeridynamicsMaterial;distanceX=5,max_neighs=50)

LJ repulsive model for 1-1 material blocks.
"""
function LJRepulsionModel(alpha, epsilon, mat1::PeridynamicsMaterial; distanceX=5, max_neighs=50)
    max_neighs = min(max_neighs,size(mat1.general.x,2))
    neighs = zeros(max_neighs, size(mat1.general.x,2))
    LJRepulsionModel11(mat1.type,mat1.general,alpha,epsilon,neighs, mat1.general.particle_size*distanceX,max_neighs)
end

"""
    repulsion_force(dr, RepMod::LJRepulsionModel12)

Calculates repulsive acceleration for 1-2 materials block interaction.
"""
function repulsion_force(dr, RepMod::LJRepulsionModel12)
    mag_dr = @_magnitude(dr) + 1.0e-10
    del_x = RepMod.equi_dist - mag_dr
    if del_x<0
        return zeros(size(dr)...)
    else
        return -(RepMod.epsilon*del_x^(RepMod.alpha-1)).*dr/mag_dr
    end
end


"""
    repulsion_force(dr, RepMod::LJRepulsionModel11)

Calculates repulsive acceleration for 1-1 materials block interaction.
"""
function repulsion_force(dr, RepMod::LJRepulsionModel11)
    del_x = -1+RepMod.material.particle_size/@_magnitude(dr)
    if del_x<0
        return zeros(size(dr)...)
    else
        return -(RepMod.epsilon*del_x^(RepMod.alpha-1)).*dr
    end
end
