export horizon_correction, dilatation, dilatation_theta, influence_function, weighted_volume, cal_family!
export _O, _I, _unit

@inline function _O(x::AbstractArray)
    return 0 * x
end

function _I(x::AbstractArray)
    return x
end

function _unit(x::Array{T,1}) where {T}
    return x / get_magnitude(x)
end

mod(X) = sqrt(mapreduce((x) -> x^2, +, X))


"""
    horizon_correction(dr,r,hr)

It gives horizon correction factor (It will give 1 as of now).
"""
function horizon_correction(dr, r, hr)
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
    return 1 / get_magnitude(dr)
end

"""
    dilatation(y::Array{Float64,2},S::GeneralMaterial,m::Array{Float64,1})

It gives dilatation as given ordinary state material model.
"""
function dilatation(y::Array{Float64,2}, S::GeneralMaterial, m::Array{Float64,1})
    return 3 * dilatation_theta(y, S) ./ m
end

"""
    dilatation_theta(y::Array{Float64,2},S::GeneralMaterial)

It gives dilatation as given ordinary state material model.
"""
function dilatation_theta(y::Array{Float64,2}, S::GeneralMaterial)
    intact = S.intact
    family = S.family
    N = size(family, 2)

    function with_if_cal_theta_ij(i, j)
        X = get_ij(j, i, S.x)
        Y = get_ij(j, i, y)
        _X = get_magnitude(X)
        _Y = get_magnitude(Y)
        extention = _Y - _X
        return influence_function(X) * _X * extention * S.volume[j] * horizon_correction(X, S.particle_size, S.horizon)
    end

    inner_map(i, js) = mapreduce((j) -> with_if_cal_theta_ij(i, j), +, js, init = 0.0)
    return Folds.map((i) -> inner_map(i, family[intact[:, i], i]), 1:N)

end


"""
    weighted_volume(S::GeneralMaterial)

It gives weighted volume.
"""
function weighted_volume(x, volume, particle_size, family, horizon)
    m = 0 * volume
    dr = zeros(3)
    for i in 1:size(x, 2)
        for k in 1:size(family, 1)
            j = family[k, i]
            if (j > 0) & (i != j)
                dr[1] = x[1, j] - x[1, i]
                dr[2] = x[2, j] - x[2, i]
                dr[3] = x[3, j] - x[3, i]
                m[i] += influence_function(dr) * get_magnitude(dr)^2 * horizon_correction(dr, particle_size, horizon) * volume[j]
            end
        end
    end
    return m
end

"""
    neigh_cells(i,j,k,N)

Calculate all neighboring cells for i,j,k cell.
"""
function neigh_cells(i, j, k, N)
    a = Vector{Int}()
    for kk in k-1:k+1
        for jj in j-1:j+1
            for ii in i-1:i+1
                if (0 < ii * jj * kk) && (ii <= N[1]) && (jj <= N[2]) && (kk <= N[3])
                    push!(a, cell_number(ii, jj, kk, N))
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
function cell_number(i, j, k, N)
    return i + (j - 1) * N[1] + (k - 1) * N[1] * N[2]
end

"""
    get_cells(x::Array{Float64,2},horizon::Float64)

Fill cells with material points.
"""
function get_cells(x::Array{Float64,2}, horizon::Number; max_part=200)
    _min = minimum(x, dims = 2)
    _max = maximum(x, dims = 2)
    N = @. Int(max(1, min(max_part, fld((_max - _min), horizon))))
    cells = [Vector{Int}() for i in 1:prod(N)]
    cell_neighs = Vector{Vector{Int}}(undef, prod(N))
    for k in 1:N[3]
        for j in 1:N[2]
            for i in 1:N[1]
                cell_neighs[cell_number(i, j, k, N)] = neigh_cells(i, j, k, N)
            end
        end
    end
    for i in 1:size(x, 2)
        ii, jj, kk = @. Int(max(1, min(max_part, fld((x[:, i] - _min), horizon))))
        # ii, jj, kk = Int.(1 .+ floor.((x[:, i] .- _min) ./ horizon))
        push!(cells[cell_number(ii, jj, kk, N)], i)
    end
    return cells, cell_neighs
end


"""
    cal_family!(family::Array{Int64,2},x::Array{Float64,2}, horizon::Float64)

Calculate family members for each material point (inplace).
"""
function cal_family!(family::Array{Int64,2}, x::Array{Float64,2}, horizon::Float64)
    cells, cell_neighs = get_cells(x, horizon)
    for cell_i in 1:length(cells)
        Threads.@threads for ca_id in 1:length(cells[cell_i])
            ind = 1
            ca = cells[cell_i][ca_id]
            a1, b1, c1 = x[1, ca], x[2, ca], x[3, ca]
            for neigh_id in 1:length(cell_neighs[cell_i])
                neighs = cell_neighs[cell_i][neigh_id]
                for fa_id in 1:length(cells[neighs])
                    fa = cells[neighs][fa_id]
                    a2, b2, c2 = x[1, fa], x[2, fa], x[3, fa]
                    if (ca != fa) & (((a1 - a2) * (a1 - a2) + (b1 - b2) * (b1 - b2) + (c1 - c2) * (c1 - c2)) < horizon * horizon)
                        family[ind, ca] = fa
                        ind += 1
                    end
                end
            end
        end
    end
    Threads.@threads for i in 1:size(family, 2)
        col = @view family[:, i]
        family[:, i] .= sort(col)
    end
    family
end


"""
    cal_family!(family::Array{Int64,2},x::Array{Float64,2}, horizon::Float64)

Calculate family members for each material point.
"""
function cal_family(x::Array{Float64,2}, horizon::Float64, max_neigh::Int64)::Array{Int64,2}
    max_neigh = min(max_neigh, size(x, 2))
    family = zeros(Int64, max_neigh, size(x, 2))
    cal_family!(family, x, horizon)
    return family
end

"""
    check directory if present.
"""


"""
    write_data(filename,every,type,pos)

It writes the data file.
"""
function write_data(filename, every, type, pos)
    file = open(filename, "w+")
    for i in 1:every:size(pos, 3)
        y = pos[:, :, i]
        N = size(y, 1)
        write(file, "$N \n\n")
        for j in 1:size(y, 1)
            t = type[j]
            a, b, c = y[j, :]
            write(file, "$j $t, $a, $b, $c ")
            write(file, "\n")
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
    N = size(y, 2)
    write(file, "$N \n\n")
    for j in 1:size(y, 2)
        t = type[j]
        a, b, c = y[1, j], y[2, j], y[3, j]
        write(file, "$j $t, $a, $b, $c ")
        write(file, "\n")
    end
    close(file)
end


"""
    write_data(filename::String,type::Array{Int64,1},y::Array{Float64,2}, v::Array{Float64,2},f::Array{Float64,2},)

It writes the data file.
"""
function write_data(filename::String, type::Array{Int64,1}, y::Array{Float64,2}, v::Array{Float64,2}, f::Array{Float64,2},)
    file = open(filename, "w+")
    N = size(y, 2)
    write(file, "$N \n\n")
    for j in 1:size(y, 2)
        t = type[j]
        a, b, c = y[1, j], y[2, j], y[3, j]
        a1, b1, c1 = v[1, j], v[2, j], v[3, j]
        a2, b2, c2 = f[1, j], f[2, j], f[3, j]
        write(file, "$j $t, $a, $b, $c, $a1, $b1, $c1, $a2, $b2, $c2 ")
        write(file, "\n")
    end
    close(file)
end

"""
    write_data(filename::String, args...)

It writes the data file.
"""
function write_data(filename::String, args...)
    file = open(filename, "w+")
    if length(size(args[1])) == 1
        N = length(args[1])
    else
        N = size(args[1], 2)
    end
    function getitems(i, arg)
        if length(size(arg)) == 1
            return arg[i]
        else
            return arg[:, i]
        end
    end

    write(file, "$N \n\n")
    for j in 1:N
        row = foldr((a, b) -> "$a $b", foldr((a, b) -> (a..., b...), [getitems(j, arg) for arg in args]))
        write(file, "$row \n")
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
    N = size(y, 2)
    write(file, "$N \n\n")
    ind = 1
    for i1 in 1:length(cells)
        for i2 in 1:length(cells[i1])
            for i3 in 1:length(cells[i1][i2])
                ids = cells[i1][i2][i3]
                for j in 1:length(ids)
                    a, b, c = y[1, ids[j]], y[2, ids[j]], y[3, ids[j]]
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
function make_matrix(S::Array{T,1}) where {T}
    n = length(S)
    N = Int64(sqrt(2 * n + 0.5^2) - 0.5)
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
function make_matrix_gm(S::Array{T,1}) where {T}
    N = length(S)
    bs = zeros(T, N, N)
    for i in 1:N
        for j in i:N
            bs[i, j] = (S[i] * S[j])^0.5
            bs[j, i] = b[i, j]
        end
    end
    return bs
end

"""
    make_NN(layers::Tuple{T}, N) where T

Create an symmetrical NxN matrix from a vector of length N(N+1)/2.    

================
## Returns
    M :: Matrix
"""
function make_NN(layers::T, N) where {T}
    L = length(layers)
    bs = []
    for i in 1:N
        for j in i:N
            push!(bs, make_NN(layers))
        end
    end
    SymMat(bs)
end

struct SymMat
    v
    N
end

function SymMat(x)
    L = length(x)
    SymMat(x, Int((-1 + sqrt(1 + 4 * 2 * L)) / 2))
end

function Base.getindex(m::SymMat, i::Int)
    m.v[i]
end

function Base.getindex(m::SymMat, i::Int, j::Int)
    @assert((i <= m.N) & (j <= m.N))
    if j > i
        i, j = j, i
    end
    ind = Int(i * (i - 1) / 2 + j)
    m.v[ind]
end

function make_NN(layers::T; act = Flux.relu) where {T}
    L = length(layers)
    Chain([Dense(layers[i], layers[i+1], act) for i in 1:L-2]..., Dense(layers[end-1], layers[end]))
end




# """
#     insert_slit!(family::Array{Int64, 2}, x::Array{T, 2}, )

# Remove connections from family members.
# """
# function insert_slit!(family::Array{Int64, 2}, x::Array{T, 2}, cutout::CutOut)
#     N, M = size(family)
#     for i in 1:M
#         x0 = x[:, i]
#         for k in 1:N
#             if family[k, i] != 0
#                 j = family[k, i]
#                 if cut(cutout, x0, x[:, j])
#                     family[k, i] = 0
#                 end
#             end
#         end
#     end
# end


# abstract type CutOut end

# struct Slit<:CutOut
#     origin
#     v1
#     v2
#     n
#     d
# end


# """
# """
# function cut(cutout::Slit, x0, x1)
#     v1 = cutout.v1
#     v2 = cutout.v2
#     o = cutout.origin
#     n = ((v1[2]*v2[3] - v1[3]*v2[2]), - (v1[1]*v2[3] - v1[3]*v2[1]), (v1[1]*v2[2] - v1[2]*v2[1])) 
#     d = - (n[1]*o[1] + n[2]*o[2] + n[3]*o[3])
#     l_, m_, n_ = x1 - x0

#     r = - (n[1]*x0[1] + n[2]*x0[2] + n[3]*x0[3] + d) / (n[1]*l_ + n[2]*m_ + n[3]*n_)
#     p = [r*l_, r*m_, r*n_] .+ x0

#     k2 = (v1[2]*(p[1] - o[1])  - v1[1]*(p[2] - o[2]))/(v2[1]*v1[2]  - v2[2]*v1[1])
#     k1 = (p[1] - o[1] - k2*v2[1]) / v1[1]
#     false
# end
