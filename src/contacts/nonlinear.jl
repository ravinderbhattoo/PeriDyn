"""
This module contains definitions for nonlinear repulsive models.

The following types and functions are exported: `NonLinearRepulsionModel`, `NonLinearRepulsionModel12`, `NonLinearRepulsionModel11`, `repulsion_force`, and `get_repulsion_force_fn`.
"""

export NonLinearRepulsionModel, NonLinearRepulsionModel12, NonLinearRepulsionModel11, repulsion_force, get_repulsion_force_fn


"""
    NonLinearRepulsionModel12(exponent::Float64, stifness::Float64, pair::Vector{AbstractVector{Int64}}, name::String, equi_dist::Float64, distance::Float64, neighs::AbstractArray{Int64, 2}, max_neighs::Int64)

Nonlinear repulsive model for 1-2 material blocks.
"""
struct NonLinearRepulsionModel12<:RepulsionModel12
    exponent::Float64
    stifness::Float64
    # pair::Vector{AbstractVector{Int64}}
    # name::String
    # equi_dist::Float64
    # distance::Float64
    # neighs::AbstractArray{Int64,2}
    # max_neighs::Int64
    @RepulsionModel12_gf
end

"""
    NonLinearRepulsionModel11(exponent::Float64, stifness::Float64, type::AbstractVector{Int64}, bid::Int64, material::GeneralMaterial, name::String, equi_dist::Float64, distance::Float64, neighs::AbstractArray{Int64, 2}, max_neighs::Int64)

Nonlinear repulsive model for 1-1 material blocks.
"""
struct NonLinearRepulsionModel11<:RepulsionModel11
    exponent::Float64
    stifness::Float64
    # type::AbstractVector{Int64}
    # bid::Int64
    # material::GeneralMaterial
    # name::String
    # equi_dist::Float64
    # distance::Float64
    # neighs::AbstractArray{Int64,2}
    # max_neighs::Int64
    @RepulsionModel11_gf
end

"""
    NonLinearRepulsionModel(exponent::Float64, stifness::Float64, mat1::PeridynamicsMaterial, mat2::PeridynamicsMaterial; distanceD=1.0, distanceX=3.0, max_neighs=50)

Constructs a nonlinear repulsive model for 1-2 material blocks.

Arguments:
- `exponent`: The exponent for the repulsive force calculation.
- `stifness`: The stiffness coefficient for the repulsive force calculation.
- `mat1`: The first PeridynamicsMaterial.
- `mat2`: The second PeridynamicsMaterial.
- `distanceD`: The distance factor for determining the search distance of neighbors within the same material.
- `distanceX`: The distance factor for determining the search distance of neighbors between different materials.
- `max_neighs`: The maximum number of neighbors to consider.

Returns:
A NonLinearRepulsionModel12 object.
"""
function NonLinearRepulsionModel(exponent, stifness,
                    mat1::PeridynamicsMaterial,
                    mat2::PeridynamicsMaterial;
                    distanceD=1.0,
                    distanceX=3.0,
                    max_neighs=50
                    )
    # Order of arguments:
    # exponent::Float64
    # stifness::Float64
    # generated using RepulsionModel12_gcal
    args = RepulsionModel12_gcal(mat1, mat2, distanceD, distanceX, max_neighs)
    NonLinearRepulsionModel12(exponent, stifness, args...)
end

"""
    NonLinearRepulsionModel(exponent::Float64, stifness::Float64, mat1::PeridynamicsMaterial; distanceX=5, max_neighs=50)

Constructs a nonlinear repulsive model for 1-1 material blocks.

Arguments:
- `exponent`: The exponent for the repulsive force calculation.
- `stifness`: The stiffness coefficient for the repulsive force calculation.
- `mat1`: The PeridynamicsMaterial.
- `distanceX`: The distance factor for determining the search distance of neighbors.
- `max_neighs`: The maximum number of neighbors to consider.

Returns:
A NonLinearRepulsionModel11 object.
"""
function NonLinearRepulsionModel(exponent, stifness,
                                mat1::PeridynamicsMaterial;
                                distanceD=1.0,
                                distanceX=5,
                                max_neighs=50)
    # Order of arguments:
    # exponent::Float64
    # stifness::Float64
    # type::AbstractVector{Int64}
    # bid::Int64
    # material::GeneralMaterial
    # name::String
    # equi_dist::Float64
    # distance::Float64
    # neighs::AbstractArray{Int64,2}
    # max_neighs::Int64
    args = RepulsionModel11_gcal(mat1, distanceD, distanceX, max_neighs)
    NonLinearRepulsionModel11(exponent, stifness, args...)
end

"""
    repulsion_force(dr, RepMod::NonLinearRepulsionModel12)

Calculates the repulsive acceleration for 1-2 material block interaction.

Arguments:
- `dr`: The vector representing the distance between particles.
- `RepMod`: The NonLinearRepulsionModel12 object.

Returns:
The repulsive acceleration as a vector.
"""
function repulsion_force(dr, RepMod::NonLinearRepulsionModel12)
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
    repulsion_force(dr, RepMod::NonLinearRepulsionModel11)

Calculates the repulsive acceleration for 1-1 material block interaction.

Arguments:
- `dr`: The vector representing the distance between particles.
- `RepMod`: The NonLinearRepulsionModel11 object.

Returns:
The repulsive acceleration as a vector.
"""
function repulsion_force(dr, RepMod::NonLinearRepulsionModel11)
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
    get_repulsion_force_fn(RepMod::NonLinearRepulsionModel11)

Returns a repulsion force function for a NonLinearRepulsionModel11 object.

Arguments:
- `RepMod`: The NonLinearRepulsionModel11 object.

Returns:
A repulsion force function that takes a distance vector and returns the repulsive acceleration.
"""
function get_repulsion_force_fn(RepMod::NonLinearRepulsionModel11)
    equi_size = RepMod.equi_dist
    expo = RepMod.exponent
    K = RepMod.stifness

    @inline function fn(dr)
        mag_dr = get_magnitude(dr) + 1.0e-10
        del_x = equi_size - mag_dr
        strain = del_x / equi_size
        if del_x < 0
            return [0.0, 0.0, 0.0]
        else
            s = ( K * strain^expo )  / mag_dr
            return s * dr
        end
    end
    return fn
end


