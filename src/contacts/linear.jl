"""
This module contains Linear replusive model definitions.
"""

export LinearRepulsionModel

"""
    LinearRepulsionModel(stifness::Float64, mat1::PeridynamicsMaterial,mat2::PeridynamicsMaterial;distanceX=5,max_neighs=50)

Linear repulsive model for 1-2 material blocks.
"""
function LinearRepulsionModel(stifness, mat1::PeridynamicsMaterial, mat2::PeridynamicsMaterial; kwargs...)
    NonLinearRepulsionModel(1, stifness, mat1::PeridynamicsMaterial, mat2::PeridynamicsMaterial; kwargs...)
end

"""
    LinearRepulsionModel(stifness::Float64, mat1::PeridynamicsMaterial;distanceX=5,max_neighs=50)

Linear repulsive model for 1-1 material blocks.
"""
function LinearRepulsionModel(stifness, mat1::PeridynamicsMaterial; kwargs...)
    NonLinearRepulsionModel(1, stifness, mat1::PeridynamicsMaterial; kwargs...)
end
