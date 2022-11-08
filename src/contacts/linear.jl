"""
This module contains Linear replusive model definitions.
"""

export LinearRepulsionModel

"""
    LinearRepulsionModel(stifness::Float64, mat1::PeridynamicsMaterial,mat2::PeridynamicsMaterial;distanceX=5,max_neighs=50)

Linear repulsive model for 1-2 material blocks.
"""
function LinearRepulsionModel(stifness, mat1::PeridynamicsMaterial,mat2::PeridynamicsMaterial;distanceX=5,max_neighs=50)
    p_size = (mat1.general.particle_size+mat2.general.particle_size)/2
    neighs = zeros(min(max_neighs,size(mat2.general.x,2)), size(mat1.general.x,2))
    NonLinearRepulsionModel12([mat1.type,mat2.type],1,stifness,p_size,neighs,p_size*distanceX,min(max_neighs,size(mat2.general.x,2)))
end

"""
    LinearRepulsionModel(stifness::Float64, mat1::PeridynamicsMaterial;distanceX=5,max_neighs=50)

Linear repulsive model for 1-1 material blocks.
"""
function LinearRepulsionModel(stifness, mat1::PeridynamicsMaterial; distanceX=5, max_neighs=50)
    max_neighs = min(max_neighs,size(mat1.general.x,2))
    neighs = zeros(max_neighs, size(mat1.general.x,2))
    NonLinearRepulsionModel11(mat1.type,mat1.blockid,mat1.general,1,stifness,neighs, mat1.general.particle_size*distanceX,max_neighs)
end
