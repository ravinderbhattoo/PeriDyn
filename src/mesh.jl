# exports
export create_block

"""
    create_block(lattice::Array{Float64,1}, N::Array{Int64,1}; rand_=0.0)

Create a rectangular point mesh.

=============
Return:
    - x
    - v
    - y
    - vol
"""
function create_block(lattice::Array{Float64,1}, N::Array{Int64,1}; rand_=0.0)
    mesh = zeros(3,prod(N))
    a = 1
    for i in 1:N[1]
        for j in 1:N[2]
            for k in 1:N[3]
                mesh[1,a] = (i + rand_*rand())*lattice[1]
                mesh[2,a] = (j + rand_*rand())*lattice[2]
                mesh[3,a] = (k + rand_*rand())*lattice[3]
                a += 1
            end
        end
    end
    # x, v, y, vol
    return mesh, zeros(size(mesh)), copy(mesh), ones(prod(N))
end
