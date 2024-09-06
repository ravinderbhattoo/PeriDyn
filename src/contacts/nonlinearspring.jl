"""
This module contains definitions for nonlinear contact models.

The following types and functions are exported: `NonLinearSpringContactModel`, `NonLinearSpringContactModel12`, `NonLinearSpringContactModel11`, `contact_force`, and `get_contact_force_fn`.
"""

export NonLinearSpringContactModel, NonLinearSpringContactModel12, NonLinearSpringContactModel11, contact_force, get_contact_force_fn


"""
    NonLinearSpringContactModel12(exponent::Float64, stifness::Float64, pair::Vector{AbstractVector{Int64}}, name::String, equi_dist::Float64, distance::Float64, neighs::AbstractArray{Int64, 2}, max_neighs::Int64)

Nonlinear contact model for 1-2 material blocks.
"""
struct NonLinearSpringContactModel12<:ContactModel12
    exponent::Float64
    stifness::QF
    # pair::Vector{AbstractVector{Int64}}
    # name::String
    # equi_dist::QF
    # distance::QF
    # neighs::AbstractArray{Int64,2}
    # max_neighs::Int64
    @ContactModel12_gf
end

"""
    NonLinearSpringContactModel11(exponent::Float64, stifness::Float64, type::AbstractVector{Int64}, bid::Int64, material::GeneralMaterial, name::String, equi_dist::Float64, distance::Float64, neighs::AbstractArray{Int64, 2}, max_neighs::Int64)

Nonlinear contact model for 1-1 material blocks.
"""
struct NonLinearSpringContactModel11<:ContactModel11
    exponent::Float64
    stifness::QF
    # type::AbstractVector{Int64}
    # bid::Int64
    # material::GeneralMaterial
    # name::String
    # equi_dist::QF
    # distance::QF
    # neighs::AbstractArray{Int64,2}
    # max_neighs::Int64
    @ContactModel11_gf
end

"""
    NonLinearSpringContactModel(exponent::Float64, stifness::Float64, mat1::PeridynamicsMaterial, mat2::PeridynamicsMaterial; distanceD=1.0, distanceX=3.0, max_neighs=50)

Constructs a nonlinear contact model for 1-2 material blocks.

Arguments:
- `exponent`: The exponent for the contact force calculation.
- `stifness`: The stiffness coefficient for the contact force calculation.
- `mat1`: The first PeridynamicsMaterial.
- `mat2`: The second PeridynamicsMaterial.
- `distanceD`: The distance factor for determining the search distance of neighbors within the same material.
- `distanceX`: The distance factor for determining the search distance of neighbors between different materials.
- `max_neighs`: The maximum number of neighbors to consider.

Returns:
A NonLinearSpringContactModel12 object.
"""
function NonLinearSpringContactModel(exponent, stifness,
                    mat1::PeridynamicsMaterial,
                    mat2::PeridynamicsMaterial;
                    distanceD=1.0,
                    distanceX=3.0,
                    max_neighs=50
                    )
    # Order of arguments:
    # exponent::Float64
    # stifness::Float64
    # generated using ContactModel12_gcal
    args = ContactModel12_gcal(mat1, mat2, distanceD, distanceX, max_neighs)
    NonLinearSpringContactModel12(exponent, stifness, args...)
end

"""
    NonLinearSpringContactModel(exponent::Float64, stifness::Float64, mat1::PeridynamicsMaterial; distanceX=5, max_neighs=50)

Constructs a nonlinear contact model for 1-1 material blocks.

Arguments:
- `exponent`: The exponent for the contact force calculation.
- `stifness`: The stiffness coefficient for the contact force calculation.
- `mat1`: The PeridynamicsMaterial.
- `distanceX`: The distance factor for determining the search distance of neighbors.
- `max_neighs`: The maximum number of neighbors to consider.

Returns:
A NonLinearSpringContactModel11 object.
"""
function NonLinearSpringContactModel(exponent, stifness,
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
    args = ContactModel11_gcal(mat1, distanceD, distanceX, max_neighs)
    NonLinearSpringContactModel11(exponent, stifness, args...)
end

"""
    get_contact_force_fn(RepMod::NonLinearSpringContactModel11)

Returns a contact force function for a NonLinearSpringContactModel11 object.

Arguments:
- `RepMod`: The NonLinearSpringContactModel11 object.

Returns:
A contact force function that takes a distance vector and returns the contact acceleration.
"""
function get_contact_force_fn(RepMod::T) where
        T<:Union{NonLinearSpringContactModel11,NonLinearSpringContactModel12}
    equi_size = RepMod.equi_dist
    expo = RepMod.exponent
    K = RepMod.stifness
    # (force/vol/vol)
    # 18K/πδ⁴ = force/area/length^4 = force/vol/vol
    zero_ = zero(K)

    @inline function fn(dr)
        mag_dr = get_magnitude(dr)
        del_x = equi_size - mag_dr
        strain = del_x / equi_size
        if del_x < zero(del_x)
            return zero_
        else
            s = ( K * strain^expo )  / (mag_dr + 1e-10*unit(mag_dr))
            return s .* dr
        end
    end
    return fn
end


