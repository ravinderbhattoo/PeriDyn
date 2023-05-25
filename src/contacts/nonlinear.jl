"""
This module contains definitions for nonlinear repulsive models.

The following types and functions are exported: `NonLinearRepulsionModel`, `NonLinearRepulsionModel12`, `NonLinearRepulsionModel11`, `repulsion_force`, and `get_repulsion_force_fn`.
"""

export NonLinearRepulsionModel, NonLinearRepulsionModel12, NonLinearRepulsionModel11, repulsion_force, get_repulsion_force_fn


"""
    NonLinearRepulsionModel12(pair::Vector{UnitRange{Int64}}, exponent::Float64, stifness::Float64, equi_dist::Float64, neighs::AbstractArray{Int64, 2}, distance::Float64, max_neighs::Int64)

Nonlinear repulsive model for 1-2 material blocks.
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

"""
    NonLinearRepulsionModel11(type::UnitRange{Int64}, bid::Int64, material::GeneralMaterial, exponent::Float64, stifness::Float64, neighs::AbstractArray{Int64, 2}, distance::Float64, max_neighs::Int64)

Nonlinear repulsive model for 1-1 material blocks.
"""
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
function NonLinearRepulsionModel(exponent,stifness, mat1::PeridynamicsMaterial,mat2::PeridynamicsMaterial;distanceD=1.0,distanceX=3.0,max_neighs=50)
    p_size = (mat1.general.particle_size+mat2.general.particle_size)/2
    neighs = zeros(min(max_neighs,size(mat2.general.x,2)), size(mat1.general.x,2))
    NonLinearRepulsionModel12([mat1.type,mat2.type],exponent,stifness,p_size*distanceD,neighs,p_size*distanceX,min(max_neighs,size(mat2.general.x,2)))
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
function NonLinearRepulsionModel(exponent,stifness, mat1::PeridynamicsMaterial; distanceX=5, max_neighs=50)
    max_neighs = min(max_neighs,size(mat1.general.x,2))
    neighs = zeros(max_neighs, size(mat1.general.x,2))
    NonLinearRepulsionModel11(mat1.type,mat1.blockid,mat1.general,exponent,stifness,neighs, mat1.general.particle_size*distanceX,max_neighs)
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
    del_x = RepMod.material.particle_size - mag_dr
    strain = del_x / RepMod.material.particle_size
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


