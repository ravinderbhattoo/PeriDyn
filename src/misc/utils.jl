"""
    make_matrix(S::Array{T,1})

Create an symmetrical NxN matrix from a vector of length N(N+1)/2.

# Arguments
- `S::Array`: Vector of length N(N+1)/2.

# Returns
- `M::Matrix`: Symmetrical NxN matrix.
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
    make_matrix(S::Matrix)

Return the matrix itself.
"""


function make_matrix(S::Matrix)
    return S
end

"""
    make_vector(M::Array{T,2})

Create an array of upper triangle of an symmetrical NxN matrix.

# Arguments
- `M::Matrix`: Symmetrical NxN matrix.

# Returns
- `S::Array`: Vector of length N(N+1)/2.
"""
function make_vector(bs::Array{T, 2}) where {T}
    n = size(bs, 1)
    N = Int64(n*(n+1)/2)
    S = zeros(T, N)
    a = 1
    for i in 1:n
        for j in i:n
            S[a] = bs[i, j]
            a += 1
        end
    end
    return S
end


"""
    make_matrix(S::Array{T,1})

Create an symmetrical NxN matrix from a vector of length N(N+1)/2.

# Arguments
- `S::Array`: Vector of length N(N+1)/2.

# Returns
- `M::Matrix`: Symmetrical NxN matrix.
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
    make_NN(layers::T, N; act=Flux.relu) where {T}

Create an symmetrical NxN matrix from a vector of length N(N+1)/2.

# Arguments
- `layers::Tuple`: Tuple of layers.
- `N::Int`: Number of layers.

# Keyword Arguments
- `act::Function`: Activation function.

# Returns
- `M::Matrix`: Symmetrical NxN matrix.
"""
function make_NN(layers::T, N; act=Flux.relu) where {T}
    L = length(layers)
    bs = []
    for i in 1:N
        for j in i:N
            push!(bs, make_NN(layers; act=act))
        end
    end
    SymMat(bs)
end

"""
    SymMat

Create an symmetrical NxN matrix from a vector of length N(N+1)/2.

# Arguments
- `v::Array`: Vector of length N(N+1)/2.
- `N::Int`: Number of layers.

# Returns
- `M::Matrix`: Symmetrical NxN matrix.
"""
struct SymMat
    v
    N
end


"""
    SymMat

Create an symmetrical NxN matrix from a vector of length N(N+1)/2.

# Arguments
- `v::Array`: Vector of length N(N+1)/2.

# Returns
- `M::Matrix`: Symmetrical NxN matrix.
"""
function SymMat(x)
    L = length(x)
    SymMat(x, Int((-1 + sqrt(1 + 4 * 2 * L)) / 2))
end

function Base.setindex!(m::SymMat, a, i::Int64, j::Int64)
    @assert((i <= m.N) & (j <= m.N))
    if j > i
        i, j = j, i
    end
    ind = Int(i * (i - 1) / 2 + j)
    m.v[ind] = a
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

"""
    make_NN(layers::T; act = Flux.relu) where {T}

Create an neural network with given layers and activation function.

# Arguments
- `layers::Tuple`: Tuple of layers.

# Keyword Arguments
- `act::Function`: Activation function.

# Returns
- `M::Matrix`: Symmetrical NxN matrix.
"""
function make_NN(layers::T; act = Flux.relu) where {T}
    L = length(layers)
    Chain([Dense(layers[i], layers[i+1], act) for i in 1:L-2]..., Dense(layers[end-1], layers[end]))
end


"""
    perturb(x::AbstractArray; e=1.0e-6)

Perturb an array with Gaussian noise.

# Arguments
- `x::AbstractArray`: Array to be perturbed.

# Keyword Arguments
- `e::Float64`: Standard deviation of the Gaussian noise.

# Returns
- `x::AbstractArray`: Perturbed array.
"""
perturb(x::AbstractArray; e=1.0e-6) = e*randn(size(x)...) .+ x

"""
    perturb(i::Int64, j::Int64; e=1.0e-6)

Create a Matrix of size i x j with Gaussian noise.

# Arguments
- `i::Int64`: Number of rows.
- `j::Int64`: Number of columns.

# Keyword Arguments
- `e::Float64`: Standard deviation of the Gaussian noise.

# Returns
- `x::AbstractArray`: Matrix of size i x j with Gaussian noise.
"""
perturb(i::Int64, j::Int64; e=1.0e-6) = e*randn(i, j)


"""
    unpack(d::Dict)

Unpack a dictionary into its components.

# Arguments
- `d::Dict`: Dictionary to be unpacked.

# Returns
- `x::Array{Float64, 1}`: x coordinates.
- `v::Array{Float64, 1}`: v coordinates.
- `y::Array{Float64, 1}`: y coordinates.
- `volume::Array{Float64, 1}`: Volume of the mesh.
- `type::Array{Int64, 1}`: Type of the mesh.
"""
function unpack(d::Dict)
    return d[:x], d[:v], d[:y], d[:volume], d[:type]
end

"""
    repack(args...; keys_ = (:x, :v, :y, :volume, :type))

Repack a dictionary from its components.

# Arguments
- `args...`: Components to be packed.

# Keyword Arguments
- `keys_ = (:x, :v, :y, :volume, :type)`: Keys of the dictionary.

# Returns
- `d::Dict`: Dictionary containing the components.
"""
function repack(args...; keys_ = (:x, :v, :y, :volume, :type))
    d = Dict()
    for i in 1:5
        d[keys_[i]] = args[i]
    end
    d
end

"""
    repack!(d::Dict, keys_, vals)

Repack a dictionary from its components inplace.

# Arguments
- `d::Dict`: Dictionary to be repacked.
- `keys_`: Keys of the dictionary.
- `vals`: Values of the dictionary.

# Returns
- `d::Dict`: Dictionary containing the components.
"""
function repack!(d::Dict, keys_, vals)
    if length(keys_)==length(vals)
        for i in 1:length(vals)
            d[keys_[i]] = vals[i]
        end
    else
        error("Lengths are not same.")
    end
    d
end

function meminfo_julia()
  # @printf "GC total:  %9.3f MiB\n" Base.gc_total_bytes(Base.gc_num())/2^20
  # Total bytes (above) usually underreports, thus I suggest using live bytes (below)
  Panel("GC live:      $(Base.gc_live_bytes()/2^20) MiB\n" *
        "JIT:          $(Base.jit_total_bytes()/2^20) MiB\n" *
        "Max. RSS:     $(Sys.maxrss()/2^20) MiB\n",
        title = "MEM USAGE"
  )
end

