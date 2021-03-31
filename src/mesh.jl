"""
This module contains standard shape point mesh functions.
"""
# imports
using ParticlesMesh

# exports
export create_block

"""
    create_block(lattice::Array{Float64,1}, N::Array{Int64,1})::Array{Float64,2}

Create a rectangular point mesh.
"""
function create_block(lattice::Array{Float64,1}, N::Array{Int64,1})
    mesh = rectangular(lattice,N)
    return mesh, zeros(size(mesh)), copy(mesh), ones(prod(N))
end
