function create_block(lattice,N)
    mesh = zeros(prod(N),3)
    a = 1
    for i in 1:N[1]
        for j in 1:N[2]
            for k in 1:N[3]
                mesh[a,:] = [i,j,k].*lattice
                a += 1
            end
        end
    end
    return mesh
end
