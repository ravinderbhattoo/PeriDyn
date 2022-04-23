"""
This module contains ShortRange replusive model definitions.
"""

export ShortRangeRepulsionModel, ShortRangeRepulsionModel11, ShortRangeRepulsionModel12, repulsion_force

function short_range_model_definition(spring_const, del_x, hor, dir)
    return 9.0 * spring_const / (hor^4 * pi) * (del_x / hor) .*  dir
end


struct ShortRangeRepulsionModel12 <: RepulsionModel12
    @RepulsionModel12_gf
    spring_const::Float64
end

struct ShortRangeRepulsionModel11 <: RepulsionModel11
    @RepulsionModel11_gf
    spring_const::Float64
end

"""
    ShortRangeRepulsionModel(spring_const::Float64, mat1::PeridynamicsMaterial,mat2::PeridynamicsMaterial;distanceX=5,max_neighs=50)

ShortRange repulsive model for 1-2 material blocks.
"""
function ShortRangeRepulsionModel(spring_const, mat1::PeridynamicsMaterial, mat2::PeridynamicsMaterial; distanceX=5, max_neighs=50)
    out = RepulsionModel12_gcal(mat1, mat2, distanceX, max_neighs)
    ShortRangeRepulsionModel12(out..., spring_const)
end

"""
    ShortRangeRepulsionModel(spring_const::Float64, mat1::PeridynamicsMaterial;distanceX=5,max_neighs=50)

ShortRange repulsive model for 1-1 material blocks.
"""
function ShortRangeRepulsionModel(spring_const, mat1::PeridynamicsMaterial; distanceX=5, max_neighs=50)
    out = RepulsionModel11_gcal(mat1, distanceX, max_neighs)
    ShortRangeRepulsionModel11(out..., spring_const)
end

"""
    repulsion_force(dr, RepMod::ShortRangeRepulsionModel12)

Calculates repulsive acceleration for 1-2 materials block interaction.
"""
function repulsion_force(dr, RepMod::ShortRangeRepulsionModel12)
    mag_dr = @_magnitude(dr) + 1.0e-10
    del_x = RepMod.equi_dist - mag_dr
    if del_x<0
        return zeros(size(dr)...)
    else
        return short_range_model_definition(RepMod.spring_const, del_x, RepMod.hor, dr/mag_dr)
    end
end


"""
    repulsion_force(dr, RepMod::ShortRangeRepulsionModel11)

Calculates repulsive acceleration for 1-1 materials block interaction.
"""
function repulsion_force(dr, RepMod::ShortRangeRepulsionModel11)
    mag_dr = @_magnitude(dr) + 1.0e-10
    del_x = RepMod.material.particle_size - mag_dr
    if del_x<0
        return zeros(size(dr)...)
    else
        return short_range_model_definition(RepMod.spring_const, del_x, RepMod.hor, dr/mag_dr)
    end
end
