"""
This module contains ShortRange replusive model definitions.
"""

export ShortRangeRepulsionModel, ShortRangeRepulsionModel11, ShortRangeRepulsionModel12, repulsion_force, short_range_rm_force

"""
    short_range_rm_force(spring_const, del_x, hor, dir)

Calculates the repulsive force using the ShortRange repulsion model formula.
The formula is: 18 * spring_const / (hor^4 * pi) * (del_x / hor) .* dir

Arguments:
- `spring_const`: Float64, the spring constant.
- `del_x`: Float64, the difference in position.
- `hor`: Float64, the characteristic length scale.
- `dir`: Vector{Float64}, the direction vector.

Returns:
- Vector{Float64}: The calculated repulsive force vector.
"""
function short_range_rm_force(spring_const, del_x, hor, dir)
    return 18 * spring_const / (hor^4 * pi) * (del_x / hor) .*  dir
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
    ShortRangeRepulsionModel(spring_const::Float64, mat1::PeridynamicsMaterial, mat2::PeridynamicsMaterial; distanceX=5, max_neighs=50)

Create a ShortRangeRepulsionModel for 1-2 material blocks.

Arguments:
- `spring_const`: Float64, the spring constant for repulsion.
- `mat1`: PeridynamicsMaterial, the first material block.
- `mat2`: PeridynamicsMaterial, the second material block.
- `distanceX`: Float64, the distance factor for neighbor searching. Default is 5.
- `max_neighs`: Int64, the maximum number of neighbors. Default is 50.

Returns:
- ShortRangeRepulsionModel12: The ShortRangeRepulsionModel for 1-2 material blocks.
"""
function ShortRangeRepulsionModel(spring_const, mat1::PeridynamicsMaterial, mat2::PeridynamicsMaterial; distanceX=5, max_neighs=50)
    out = RepulsionModel12_gcal(mat1, mat2, distanceX, max_neighs)
    ShortRangeRepulsionModel12(out..., spring_const)
end

"""
    ShortRangeRepulsionModel(spring_const::Float64, mat1::PeridynamicsMaterial; distanceX=5, max_neighs=50)

Create a ShortRangeRepulsionModel for 1-1 material blocks.

Arguments:
- `spring_const`: Float64, the spring constant for repulsion.
- `mat1`: PeridynamicsMaterial, the material block.
- `distanceX`: Float64, the distance factor for neighbor searching. Default is 5.
- `max_neighs`: Int64, the maximum number of neighbors. Default is 50.

Returns:
- ShortRangeRepulsionModel11: The ShortRangeRepulsionModel for 1-1 material blocks.
"""
function ShortRangeRepulsionModel(spring_const, mat1::PeridynamicsMaterial; distanceX=5, max_neighs=50)
    out = RepulsionModel11_gcal(mat1, distanceX, max_neighs)
    ShortRangeRepulsionModel11(out..., spring_const)
end

"""
    repulsion_force(dr, RepMod::ShortRangeRepulsionModel12)

Calculates the repulsive acceleration for 1-2 materials block interaction.

Arguments:
- `dr`: Vector{Float64}, the displacement vector between particles.
- `RepMod`: ShortRangeRepulsionModel12, the repulsion model for 1-2 material blocks.

Returns:
- Vector{Float64}: The repulsive acceleration vector.
"""
function repulsion_force(dr, RepMod::ShortRangeRepulsionModel12)
    mag_dr = @_magnitude(dr) + 1.0e-10
    del_x = RepMod.equi_dist - mag_dr
    if del_x<0
        return zeros(size(dr)...)
    else
        return short_range_rm_force(RepMod.spring_const, del_x, RepMod.hor, dr/mag_dr)
    end
end


"""
    repulsion_force(dr, RepMod::ShortRangeRepulsionModel11)

Calculates the repulsive acceleration for 1-1 materials block interaction.

Arguments:
- `dr`: Vector{Float64}, the displacement vector between particles.
- `RepMod`: ShortRangeRepulsionModel11, the repulsion model for 1-1 material blocks.

Returns:
- Vector{Float64}: The repulsive acceleration vector.
"""
function repulsion_force(dr, RepMod::ShortRangeRepulsionModel11)
    mag_dr = @_magnitude(dr) + 1.0e-10
    del_x = RepMod.material.particle_size - mag_dr
    if del_x<0
        return zeros(size(dr)...)
    else
        return short_range_rm_force(RepMod.spring_const, del_x, RepMod.hor, dr/mag_dr)
    end
end
