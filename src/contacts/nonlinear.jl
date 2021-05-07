"""
This module contains NonLinear replusive model definitions.
"""

export NonLinearRepulsionModel, NonLinearRepulsionModel11, NonLinearRepulsionModel12, repulsion_acc


struct NonLinearRepulsionModel12<:RepulsionModel12
    pair::Array{Int64,1}
    exponent::Float64
    stifness::Float64
    equi_dist::Float64
    neighs::Array{Int64,2}
    distance::Float64
    max_neighs::Int64
end

struct NonLinearRepulsionModel11<:RepulsionModel11
    type::UnitRange{Int64}
    material::GeneralMaterial
    exponent::Float64
    stifness::Float64
    neighs::Array{Int64,2}
    distance::Float64
    max_neighs::Int64
end

"""
    NonLinearRepulsionModel(exponent::Float64,stifness::Float64, mat1::PeridynamicsMaterial,mat2::PeridynamicsMaterial;distanceX=5,max_neighs=50)

NonLinear repulsive model for 1-2 material blocks.
"""
function NonLinearRepulsionModel(exponent::Float64,stifness::Float64, mat1::PeridynamicsMaterial,mat2::PeridynamicsMaterial;distanceX=5,max_neighs=50)
    p_size = (mat1.general.particle_size+mat2.general.particle_size)/2
    neighs = zeros(min(max_neighs,size(mat2.general.x,2)), size(mat1.general.x,2))
    NonLinearRepulsionModel12([mat1.type,mat2.type],exponent,stifness,p_size,neighs,p_size*distanceX,min(max_neighs,size(mat2.general.x,2)))
end

"""
    NonLinearRepulsionModel(exponent::Float64,stifness::Float64, mat1::PeridynamicsMaterial;distanceX=5,max_neighs=50)

NonLinear repulsive model for 1-1 material blocks.
"""
function NonLinearRepulsionModel(exponent::Float64,stifness::Float64, mat1::PeridynamicsMaterial; distanceX=5, max_neighs=50)
    max_neighs = min(max_neighs,size(mat1.general.x,2))
    neighs = zeros(max_neighs, size(mat1.general.x,2))
    NonLinearRepulsionModel11(mat1.type,mat1.general,exponent,stifness,neighs, mat1.general.particle_size*distanceX,max_neighs)
end

"""
    repulsion_acc(dr,den_i,RepMod::NonLinearRepulsionModel12)

Calculates repulsive acceleration for 1-2 materials block interaction.
"""
function repulsion_acc(dr,den_i,RepMod::NonLinearRepulsionModel12)
    mag_dr = magnitude(dr) + 1.0e-10
    del_x = RepMod.equi_dist - mag_dr
    if del_x<0
        return zeros(size(dr)...)
    else
        return -(RepMod.stifness*del_x^RepMod.exponent).*dr/mag_dr
    end
end


"""
    repulsion_acc(dr,den_i,RepMod::NonLinearRepulsionModel11)

Calculates repulsive acceleration for 1-1 materials block interaction.
"""
function repulsion_acc(dr,RepMod::NonLinearRepulsionModel11)
    del_x = -1+RepMod.material.particle_size/magnitude(dr)
    if del_x<0
        return zeros(size(dr)...)
    else
        return -(RepMod.stifness*del_x^RepMod.exponent).*dr
    end
end
