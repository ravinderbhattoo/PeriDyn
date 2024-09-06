"""
The `LinearSpringContactModel` module contains definitions for linear contact models.
It provides constructors for creating linear contact models for 1-2 and 1-1 material blocks.
The constructors accept the stiffness coefficient and PeridynamicsMaterial objects representing the material blocks.
Additional optional keyword arguments can be provided to customize the distance parameter and the maximum number of neighbors.
"""

export LinearSpringContactModel


"""
    LinearSpringContactModel(stiffness::Float64, mat1::PeridynamicsMaterial, mat2::PeridynamicsMaterial; distanceX=5, max_neighs=50)

Constructs a linear contact model for 1-2 material blocks. The constructors internally call the `NonLinearSpringContactModel` constructor with a nonlinearity parameter of 1 and the provided arguments.


Arguments:
- `stiffness`: The stiffness coefficient for the contact model.
- `mat1`: PeridynamicsMaterial representing the first material block.
- `mat2`: PeridynamicsMaterial representing the second material block.
- `distanceX`: Optional keyword argument specifying the distance parameter (default: 5).
- `max_neighs`: Optional keyword argument specifying the maximum number of neighbors (default: 50).

Returns:
- An instance of the LinearSpringContactModel.

"""
function LinearSpringContactModel(stifness, mat1::PeridynamicsMaterial, mat2::PeridynamicsMaterial; kwargs...)
    NonLinearSpringContactModel(1, stifness, mat1::PeridynamicsMaterial, mat2::PeridynamicsMaterial; kwargs...)
end

"""
    LinearSpringContactModel(stiffness::Float64, mat1::PeridynamicsMaterial; distanceX=5, max_neighs=50)

Constructs a linear contact model for 1-1 material blocks. The constructors internally call the `NonLinearSpringContactModel` constructor with a nonlinearity parameter of 1 and the provided arguments.

Arguments:
- `stiffness`: The stiffness coefficient for the contact model.
- `mat1`: PeridynamicsMaterial representing the material block.
- `distanceX`: Optional keyword argument specifying the distance parameter (default: 5).
- `max_neighs`: Optional keyword argument specifying the maximum number of neighbors (default: 50).

Returns:
- An instance of the LinearSpringContactModel.

"""
function LinearSpringContactModel(stifness, mat1::PeridynamicsMaterial; kwargs...)
    NonLinearSpringContactModel(1, stifness, mat1::PeridynamicsMaterial; kwargs...)
end
