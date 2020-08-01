"""
This module contains standard shape point mesh functions.
"""

export create_block

"""
    create_block(lattice::Array{Float64,1}, N::Array{Int64,1})::Array{Float64,2}

Create a rectangular point mesh.
"""
function create_block(lattice::Array{Float64,1}, N::Array{Int64,1})
    mesh = zeros(3,prod(N))
    a = 1
    for i in 1:N[1]
        for j in 1:N[2]
            for k in 1:N[3]
                mesh[1,a] = i*lattice[1]
                mesh[2,a] = j*lattice[2]
                mesh[3,a] = k*lattice[3]
                a += 1
            end
        end
    end
    return mesh, zeros(size(mesh)), copy(mesh), ones(prod(N))
end
