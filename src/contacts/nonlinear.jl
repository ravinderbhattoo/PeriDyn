"""
This module contains NonLinear replusive model definitions.
"""

export NonLinearRepulsionModel, NonLinearRepulsionModel11, NonLinearRepulsionModel12, repulsion_force


"""
"""
struct NonLinearRepulsionModel12<:RepulsionModel12
    pair::Vector{UnitRange{Int64}}
    exponent::Float64
    stifness::Float64
    equi_dist::Float64
    neighs::AbstractArray{Int64,2}
    distance::Float64
    max_neighs::Int64
end

struct NonLinearRepulsionModel11<:RepulsionModel11
    type::UnitRange{Int64}
    bid::Int64
    material::GeneralMaterial
    exponent::Float64
    stifness::Float64
    neighs::AbstractArray{Int64,2}
    distance::Float64
    max_neighs::Int64
end

"""
    NonLinearRepulsionModel(exponent::Float64,stifness::Float64, mat1::PeridynamicsMaterial,mat2::PeridynamicsMaterial;distanceX=5,max_neighs=50)

NonLinear repulsive model for 1-2 material blocks.
"""
function NonLinearRepulsionModel(exponent,stifness, mat1::PeridynamicsMaterial,mat2::PeridynamicsMaterial;distanceD=1.0,distanceX=3.0,max_neighs=50)
    p_size = (mat1.general.particle_size+mat2.general.particle_size)/2
    neighs = zeros(min(max_neighs,size(mat2.general.x,2)), size(mat1.general.x,2))
    NonLinearRepulsionModel12([mat1.type,mat2.type],exponent,stifness,p_size*distanceD,neighs,p_size*distanceX,min(max_neighs,size(mat2.general.x,2)))
end

"""
    NonLinearRepulsionModel(exponent::Float64,stifness::Float64, mat1::PeridynamicsMaterial;distanceX=5,max_neighs=50)

NonLinear repulsive model for 1-1 material blocks.
"""
function NonLinearRepulsionModel(exponent,stifness, mat1::PeridynamicsMaterial; distanceX=5, max_neighs=50)
    max_neighs = min(max_neighs,size(mat1.general.x,2))
    neighs = zeros(max_neighs, size(mat1.general.x,2))
    NonLinearRepulsionModel11(mat1.type,mat1.blockid,mat1.general,exponent,stifness,neighs, mat1.general.particle_size*distanceX,max_neighs)
end

"""
    repulsion_force(dr,RepMod::NonLinearRepulsionModel12)

Calculates repulsive acceleration for 1-2 materials block interaction.
"""
function repulsion_force(dr,RepMod::NonLinearRepulsionModel12)
    mag_dr = get_magnitude(dr) + 1.0e-10
    del_x = RepMod.equi_dist - mag_dr
    strain = del_x / RepMod.equi_dist
    if del_x<0
        return zeros(size(dr)...)
    else
        return ( RepMod.stifness * strain^RepMod.exponent )  .*dr/mag_dr
    end
end


"""
    repulsion_force(dr,RepMod::NonLinearRepulsionModel11)

Calculates repulsive acceleration for 1-1 materials block interaction.
"""
function repulsion_force(dr, RepMod::NonLinearRepulsionModel11)
    mag_dr = get_magnitude(dr) + 1.0e-10
    del_x = RepMod.material.particle_size - mag_dr
    strain = del_x / RepMod.material.particle_size
    if del_x<0
        return zeros(size(dr)...)
    else
        return ( RepMod.stifness * strain^RepMod.exponent )  .*dr/mag_dr
    end
end


function get_repulsion_force_fn(RepMod::NonLinearRepulsionModel11)
    equi_size = RepMod.material.particle_size
    expo = RepMod.exponent
    K = RepMod.stifness

    function fn(dr)
        mag_dr = get_magnitude(dr) + 1.0e-10
        del_x = equi_size - mag_dr
        strain = del_x / equi_size
        if del_x < 0
            return (0.0, 0.0, 0.0)
        else
            s = ( K * strain^expo )  / mag_dr
            return (s*dr[1], s*dr[2], s*dr[3])
        end
    end
    return fn
end
