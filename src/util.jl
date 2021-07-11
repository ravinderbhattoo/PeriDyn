
export horizon_correction, dilatation, dilatation_theta, influence_function, weighted_volume, cal_family!

mod(X) = sqrt(sum(X.^2))
ϵ(Y, X) = mod(Y) - mod(X)
ϵⁱ(Y, X) = θ(Y, S) * mod(X) / 3
ϵᵈ(Y, X) = ϵ(Y, X) - ϵⁱ(Y, X)

function _θ(y::Array{Float64,2}, S::GeneralMaterial)
    theta = zeros(Float64, size(S.volume)...)
    j = 1::Int64
    E_n = 0.0::Float64
    E = 0.0::Float64
    x = [1.0,0,0]::Array{Float64}
    e = 0.0::Float64
    for i in 1:size(S.x,2)
        for k in 1:size(S.family,1)
            if S.intact[k,i]
                j = S.family[k,i]
                if (j>0) & (i!=j)
                    x = _ij(j,i,S.x)
                    E = _magnitude(x)
                    E_n = _magnitude(_ij(j,i,y))
                    e = E_n - E
                    theta[i] += influence_function(x)*E*e*S.volume[j] *horizon_correction(x,S.particle_size,S.horizon)
                end
            end
        end
    end
    return theta
end

function θ(Y, S)
    3 * _θ(Y, S) ./ m
end



"""
    horizon_correction(dr,r,hr)

It gives horizon correction factor (It will give 1 as of now).
"""
function horizon_correction(dr,r,hr)
    """
    Horizon correction for perdynamic model.
    Args
        dr vector_ij.
    """
    return 1
end

"""
    influence_function(dr)

It gives influence function factor (It will give 1/r as of now).
"""
function influence_function(dr)
    return 1/_magnitude(dr)
end

"""
    dilatation(y::Array{Float64,2},S::GeneralMaterial,m::Array{Float64,1})

It gives dilatation as given ordinary state material model.
"""
function dilatation(y::Array{Float64,2},S::GeneralMaterial,m::Array{Float64,1})
    return 3*dilatation_theta(y,S)./m
end

"""
    dilatation_theta(y::Array{Float64,2},S::GeneralMaterial)

It gives dilatation as given ordinary state material model.
"""
function dilatation_theta(y::Array{Float64,2},S::GeneralMaterial)
    theta = zeros(Float64,size(S.volume)...)
    j = 1::Int64
    E_n = 0.0::Float64
    E = 0.0::Float64
    x = [1.0,0,0]::Array{Float64}
    e = 0.0::Float64
    for i in 1:size(S.x,2)
        for k in 1:size(S.family,1)
            if S.intact[k,i]
                j = S.family[k,i]
                if (j>0) & (i!=j)
                    x = _ij(j,i,S.x)
                    E = _magnitude(x)
                    E_n = _magnitude(_ij(j,i,y))
                    e = E_n - E
                    theta[i] += influence_function(x)*E*e*S.volume[j] *horizon_correction(x,S.particle_size,S.horizon)
                end
            end
        end
    end
    return theta
end


"""
    weighted_volume(S::GeneralMaterial)

It gives weighted volume.
"""
function weighted_volume(x, volume, particle_size, family, horizon)
    m = 0*volume::Array{Float64,1}
    dr = zeros(3)::Array{Float64,1}
    j = 0::Int64
    for i in 1:size(x,2)
        for k in 1:size(family,1)
            j = family[k,i]
            if (j>0) & (i!=j)
                dr[1] = x[1,j]-x[1,i]
                dr[2] = x[2,j]-x[2,i]
                dr[3] = x[3,j]-x[3,i]
                m[i] += influence_function(dr)*_magnitude(dr)^2*horizon_correction(dr, particle_size,horizon)*volume[j]
            end
        end
    end
    return m
end


"""
    neigh_cells(i,j,k,N)

Calculate all neighboring cells for i,j,k cell.
"""
function neigh_cells(i,j,k,N)
    a = Vector{Int}()
    for kk in k-1:k+1
        for jj in j-1:j+1
            for ii in i-1:i+1
                if (0<ii*jj*kk) && (ii<=N[1]) && (jj<=N[2]) && (kk<=N[3])
                    push!(a,cell_number(ii,jj,kk,N))
                end
            end
        end
    end
    return a
end

"""
    cell_number(i,j,k,N)

Calculates cell number for given i,j,k cell indices (when in a list).
"""
function cell_number(i,j,k,N)
    return i+(j-1)*N[1]+(k-1)*N[1]*N[2]
end

"""
    get_cells(x::Array{Float64,2},horizon::Float64)

Fill cells with material points.
"""
function get_cells(x::Array{Float64,2}, horizon::Float64)
    _min = minimum(x,dims=2)
    _max = maximum(x,dims=2)
    N = Int.(1 .+ floor.((_max-_min)/horizon))
    cells = [Vector{Int}() for i in 1:prod(N)]
    cell_neighs = Vector{Vector{Int}}(undef,prod(N))
    for k in 1:N[3]
        for j in 1:N[2]
            for i in 1:N[1]
                cell_neighs[cell_number(i,j,k,N)] = neigh_cells(i,j,k,N)
            end
        end
    end
    for i in 1:size(x,2)
        ii,jj,kk = Int.(1 .+ floor.((x[:,i].-_min)/horizon))
        push!(cells[cell_number(ii,jj,kk,N)],i)
    end
    return cells, cell_neighs
end


"""
    cal_family!(family::Array{Int64,2},x::Array{Float64,2}, horizon::Float64)

Calculate family members for each material point (inplace).
"""
function cal_family!(family::Array{Int64,2},x::Array{Float64,2}, horizon::Float64)

    cells, cell_neighs = get_cells(x,horizon)
    for cell_i in 1:length(cells)
        for ca_id in 1:length(cells[cell_i])
            ind = 1
            ca = cells[cell_i][ca_id]
            a1,b1,c1 = x[1,ca],x[2,ca],x[3,ca]
            for neigh_id in 1:length(cell_neighs[cell_i])
                neighs = cell_neighs[cell_i][neigh_id]
                for fa_id in 1:length(cells[neighs])
                    fa = cells[neighs][fa_id]
                    a2,b2,c2 = x[1,fa], x[2,fa], x[3,fa]
                    if 1e-4<((a1-a2)*(a1-a2)+(b1-b2)*(b1-b2)+(c1-c2)*(c1-c2))<horizon^2
                        family[ind,ca] = fa
                        ind += 1
                    end
                end
            end
        end
    end
    sort!(family,dims=1)
end


"""
    cal_family!(family::Array{Int64,2},x::Array{Float64,2}, horizon::Float64)

Calculate family members for each material point.
"""
function cal_family(x::Array{Float64,2}, horizon::Float64, max_neigh::Int64)::Array{Int64,2}
    max_neigh = min(max_neigh,size(x,2))
    family = zeros(Int64,max_neigh,size(x,2))
    cal_family!(family,x,horizon)
    return family
end

"""
    write_data(filename,every,type,pos)

It writes the data file.
"""
function write_data(filename, every, type, pos)
    file = open(filename, "w+")
    for i in 1:every:size(pos,3)
        y = pos[:,:,i]
        N = size(y,1)
        write(file, "$N \n\n")
        for j in 1:size(y,1)
            t = type[j]
            a,b,c = y[j,:]
            write(file, "$t, $a, $b, $c ")
            write(file,"\n")
        end
    end
    close(file)
end

"""
    write_data(filename::String,type::Array{Int64,1},y::Array{Float64,2})

It writes the data file.
"""
function write_data(filename::String, type::Array{Int64,1}, y::Array{Float64,2})
    file = open(filename, "w+")
    N = size(y,2)
    write(file, "$N \n\n")
    for j in 1:size(y,2)
        t = type[j]
        a,b,c = y[1,j],y[2,j],y[3,j]
        write(file, "$t, $a, $b, $c ")
        write(file,"\n")
    end
    close(file)
end


"""
    write_data(filename::String,type::Array{Int64,1},y::Array{Float64,2}, v::Array{Float64,2},f::Array{Float64,2},)

It writes the data file.
"""
function write_data(filename::String, type::Array{Int64,1}, y::Array{Float64,2}, v::Array{Float64,2}, f::Array{Float64,2},)
    file = open(filename, "w+")
    N = size(y,2)
    write(file, "$N \n\n")
    for j in 1:size(y,2)
        t = type[j]
        a,b,c = y[1,j],y[2,j],y[3,j]
        a1,b1,c1 = v[1,j],v[2,j],v[3,j]
        a2,b2,c2 = f[1,j],f[2,j],f[3,j]
        write(file, "$t, $a, $b, $c, $a1, $b1, $c1, $a2, $b2, $c2 ")
        write(file,"\n")
    end
    close(file)
end


"""
    write_data_cell_ids(filename::String, y::Array{Float64,2}, cells::Any)

It writes the data file.
"""
function write_data_cell_ids(filename::String, y::Matrix, horizon)
    cells, _ = get_cells(y, horizon)
    file = open(filename, "w+")
    N = size(y,2)
    write(file, "$N \n\n")
    ind = 1
    for i1 in 1:length(cells)
        for i2 in 1:length(cells[i1])
            for i3 in 1:length(cells[i1][i2])
                ids = cells[i1][i2][i3]
                for j in 1:length(ids)
                    a,b,c = y[1,ids[j]],y[2,ids[j]],y[3,ids[j]]
                    write(file, "$ind, $a, $b, $c \n")
                end
                ind += 1
            end
        end
    end
    close(file)
end


"""
    make_matrix(S::Array{T,1})

Create an symmetrical NxN matrix from a vector of length N(N+1)/2.    

================
## Returns
    M :: Matrix
"""
function make_matrix(S::Array{T,1}) where T
    n = length(S)
    N = Int64(sqrt(2*n + 0.5^2) - 0.5)
    bs = zeros(T, N, N)
    a = 1
    for i in 1:N 
        for j in i:N 
            bs[i, j] = S[a]
            bs[j, i] = S[a]
            a += 1 
        end 
    end 
    return bs
end


"""
    make_matrix(S::Array{T,1})

Create an symmetrical NxN matrix from a vector of length N(N+1)/2.    

================
## Returns
    M :: Matrix
"""
function make_matrix_gm(S::Array{T,1}) where T
    N = length(S)
    bs = zeros(T, N, N)
    for i in 1:N 
        for j in i:N 
            bs[i, j] = (S[i]*S[j])^0.5
            bs[j, i] = b[i, j]
        end 
    end 
    return bs
end
