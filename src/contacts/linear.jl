"""
The `LinearRepulsionModel` module contains definitions for linear repulsive models.
It provides constructors for creating linear repulsive models for 1-2 and 1-1 material blocks.
The constructors accept the stiffness coefficient and PeridynamicsMaterial objects representing the material blocks.
Additional optional keyword arguments can be provided to customize the distance parameter and the maximum number of neighbors.
"""

export LinearRepulsionModel


"""
    LinearRepulsionModel(stiffness::Float64, mat1::PeridynamicsMaterial, mat2::PeridynamicsMaterial; distanceX=5, max_neighs=50)

Constructs a linear repulsive model for 1-2 material blocks. The constructors internally call the `NonLinearRepulsionModel` constructor with a nonlinearity parameter of 1 and the provided arguments.


Arguments:
- `stiffness`: The stiffness coefficient for the repulsive model.
- `mat1`: PeridynamicsMaterial representing the first material block.
- `mat2`: PeridynamicsMaterial representing the second material block.
- `distanceX`: Optional keyword argument specifying the distance parameter (default: 5).
- `max_neighs`: Optional keyword argument specifying the maximum number of neighbors (default: 50).

Returns:
- An instance of the LinearRepulsionModel.

"""
function LinearRepulsionModel(stifness, mat1::PeridynamicsMaterial, mat2::PeridynamicsMaterial; kwargs...)
    NonLinearRepulsionModel(1, stifness, mat1::PeridynamicsMaterial, mat2::PeridynamicsMaterial; kwargs...)
end

"""
    LinearRepulsionModel(stiffness::Float64, mat1::PeridynamicsMaterial; distanceX=5, max_neighs=50)

Constructs a linear repulsive model for 1-1 material blocks. The constructors internally call the `NonLinearRepulsionModel` constructor with a nonlinearity parameter of 1 and the provided arguments.

Arguments:
- `stiffness`: The stiffness coefficient for the repulsive model.
- `mat1`: PeridynamicsMaterial representing the material block.
- `distanceX`: Optional keyword argument specifying the distance parameter (default: 5).
- `max_neighs`: Optional keyword argument specifying the maximum number of neighbors (default: 50).

Returns:
- An instance of the LinearRepulsionModel.

"""
function LinearRepulsionModel(stifness, mat1::PeridynamicsMaterial; kwargs...)
    NonLinearRepulsionModel(1, stifness, mat1::PeridynamicsMaterial; kwargs...)
end
