
export horizon_correction, dilatation, dilatation_theta, dilatation!,
       influence_function, weighted_volume, cal_family!, loop_over_neighs!
export __O, __I, __unit, __norm


"""
    __O(x::AbstractArray)

Returns zero of the same type as of x.

# Arguments
- `x::AbstractArray`: input array

# Returns
- `zero(x)`: zero of the same type as of x
"""
@inline function __O(x::AbstractArray)
    return zero(x)
end

"""
    __I(x::AbstractArray)

Returns x of the same type as of x (identity function). It is not a true identity
function as it does not return a copy of x but returns x itself.

# Arguments
- `x::AbstractArray`: input array

# Returns
- `x`: x of the same type as of x (identity function)
"""
@inline function __I(x::AbstractArray)
    return x
end

"""
    __norm(x::AbstractArray)

Returns magnitude of x (``\\lVert x \\rVert``).

# Arguments
- `x::AbstractArray`: input array

# Returns
- `y`: magnitude of x
"""
@inline function __norm(x::AbstractArray)
    @assert N == SPATIAL_DIMENSIONS_REF[]
    return get_magnitude(x)
end

"""
    __unit(x::AbstractArray)

Returns unit vector of x.

# Arguments
- `x::AbstractArray`: input array

# Returns
- `y`: unit vector of x
"""
@inline function __unit(x::Array{T, N}) where {T, N}
    return x / __norm(x)
end


"""
    horizon_correction(dr, ps, hr)

It gives horizon correction factor (It will give 1 as of now).

# Arguments
- `dr`: ``X_{ij}``
- `ps`: particle size
- `hr`: horizon

# Returns
- `s`: horizon correction factor
"""
@inline function horizon_correction(dr, ps, hr)
    # r = ps / 2
    # s = (hr - get_magnitude(dr)) / r
    # return max(-1.0, min(1.0, s)) * 0.5 + 0.5
    return 1
end

"""
    influence_function(dr)

It gives influence function as given ordinary state material model.
Return ``\\frac{1}{ \\lVert dr \\rVert}`` as of now.

# Arguments
- `dr`: ``X_{ij}``

# Returns
- ``\\omega_{ij}``: influence function
"""
function influence_function(dr)
    return 1 / get_magnitude(dr)
end

"""
    dilatation!(y, mat, device)

It gives dilatation as given by ordinary state material model.
Calls `dilatation!(theta, y, x, intact, family, volume, m, particle_size, horizon, device)`.

# Arguments
- `y`: deformed positions
- `mat`: material model
- `device`: device type

# Returns
- `nothing`: In-place operation
"""
function dilatation!(y, mat, device)
    gen = mat.general
    dilatation!(mat.specific.theta, y, gen.x, gen.intact, gen.family, gen.volume, gen.weighted_volume, gen.particle_size, gen.horizon, device)
end

"""
    dilatation!(theta, y, x, intact, family, volume, m, particle_size, horizon, device::Type{Val{:cpu}})

It gives dilatation as given ordinary state material model.

``\\theta_i = \\frac{3}{m_i} \\sum_{j \\in \\mathcal{N}_i} \\omega_{ij} \\times \\lVert X_{ij} \\rVert \\times e_{ij} \\times vol_j \\times s_{ij}``

where,
``\\omega_{ij}`` is influence function,
``X_{ij}`` is ``x_j - x_i``,
``e_{ij}`` is bond extention,
``vol_j`` is particle volume,
``s_{ij}``is horizon correction factor.

# Arguments
- `theta`: dilatation
- `y`: deformed positions
- `x`: initial positions
- `intact`: intact bonds
- `family`: family of particles
- `volume`: particle volume
- `m`: weighted volume
- `particle_size`: particle size
- `horizon`: horizon
- `device::Type{Val{:cpu}}`: device type

# Returns
- `nothing`: `theta` is updated in-place

"""
function dilatation!(theta, y, x, intact, family, volume, m, particle_size, horizon, device::Type{Val{:cpu}})
    M, N = size(family)
    fill!(theta, zero(eltype(theta)))
    Threads.@threads for i in 1:N
        c_i = 3 / m[i]
        for k in 1:M
            if intact[k, i] == 1
                j = family[k, i]
                X = get_ij(j, i, x)
                Y = get_ij(j, i, y)
                xij = get_magnitude(X)
                yij = get_magnitude(Y)
                extension = yij - xij
                theta[i] += c_i * influence_function(X) * xij * extension * volume[j] * horizon_correction(X, particle_size, horizon)
            end
        end
    end
end

"""
    dilatation(y, x, intact, family, volume, m, particle_size, horizon, device::Type{Val{:cpu}})

It gives dilatation as given ordinary state material model.

``\\theta_i = \\frac{3}{m_i} \\sum_{j \\in \\mathcal{N}_i} \\omega_{ij} \\times \\lVert X_{ij} \\rVert \\times e_{ij} \\times vol_j \\times s_{ij}``

where,
``\\omega_{ij}`` is influence function,
``X_{ij}`` is ``x_j - x_i``,
``e_{ij}`` is bond extention,
``vol_j`` is particle volume,
``s_{ij}``is horizon correction factor.

# Arguments
- `y`: deformed positions
- `x`: initial positions
- `intact`: intact bonds
- `family`: family of particles
- `volume`: particle volume
- `m`: mass
- `particle_size`: particle size
- `horizon`: horizon
- `device::Type{Val{:cpu}}`: device type

# Returns
- `theta`: dilatation

"""
function dilatation(y, x, intact, family, volume, m, particle_size, horizon, device::Type{Val{:cpu}})
    N = size(family, 2)
    function with_if_cal_theta_ij(i, j)
        X = get_ij(j, i, x)
        Y = get_ij(j, i, y)
        xij = get_magnitude(X)
        yij = get_magnitude(Y)
        extension = _Y - _X
        return influence_function(X) * xij * extension * volume[j] * horizon_correction(X, particle_size, horizon)
    end

    inner_map(i, js) = map_reduce((j) -> with_if_cal_theta_ij(i, j), +, js)
    out = map((i) -> 3/m[i]*inner_map(i, family[intact[:, i], i]), 1:N)
    return out
end


"""
    dilatation(y, x, intact, family, volume, m, particle_size, horizon, device::Type{Val{:cuda}})

It gives dilatation as given ordinary state material model.

``\\theta_i = \\frac{3}{m_i} \\sum_{j \\in \\mathcal{N}_i} \\omega_{ij} \\times \\lVert X_{ij} \\rVert \\times e_{ij} \\times vol_j \\times s_{ij}``

where,
``\\omega_{ij}`` is influence function,
``X_{ij}`` is ``x_j - x_i``,
``e_{ij}`` is bond extention,
``vol_j`` is particle volume,
``s_{ij}``is horizon correction factor.

# Arguments
- `y`: deformed positions
- `x`: initial positions
- `intact`: intact bonds
- `family`: family of particles
- `volume`: particle volume
- `m`: mass
- `particle_size`: particle size
- `horizon`: horizon
- `device::Type{Val{:cuda}}`: device type

# Returns
- `theta`: dilatation
"""
function dilatation(y, x, intact, family, volume, m, particle_size, horizon, device::Type{Val{:cuda}})
    theta = CUDA.CuArray(zeros(eltype(volume), size(volume)))
    dilatation!(theta, y, x, intact, family, volume, m, particle_size, horizon, device)
    return theta
end

"""
    dilatation!(theta, y, x, intact, family, volume, m, particle_size, horizon, device::Type{Val{:cuda}})

It gives dilatation as given ordinary state material model.

``\\theta_i = \\frac{3}{m_i} \\sum_{j \\in \\mathcal{N}_i} \\omega_{ij} \\times \\lVert X_{ij} \\rVert \\times e_{ij} \\times vol_j \\times s_{ij}``

where,
``\\omega_{ij}`` is influence function,
``X_{ij}`` is ``x_j - x_i``,
``e_{ij}`` is bond extention,
``vol_j`` is particle volume,
``s_{ij}``is horizon correction factor.

# Arguments
- `theta`: dilatation
- `y`: deformed positions
- `x`: initial positions
- `intact`: intact bonds
- `family`: family of particles
- `volume`: particle volume
- `m`: mass
- `particle_size`: particle size
- `horizon`: horizon
- `device::Type{Val{:cuda}}`: device type

# Returns
- `nothing`: `theta` is updated in-place
"""
function dilatation!(theta, l1, y, x, intact, family, volume, m, particle_size, horizon, device::Type{Val{:cuda}})
    function cal_theta(theta, y, x, intact, family, volume, m)
        index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        stride = blockDim().x * gridDim().x
        for i in index:stride:length(volume)
            theta[i] = 0.0
            _i = i
            for k in 1:size(family, 1)
                if intact[k,i]
                    j = family[k,i]
                    _j = j
                    _X = sqrt((x[1,j]-x[1,i])^2 + (x[2,j]-x[2,i])^2 + (x[3,j]-x[3,i])^2)
                    _Y = sqrt((y[1,_j]-y[1,_i])^2 + (y[2,_j]-y[2,_i])^2 + (y[3,_j]-y[3,_i])^2)
                    extension = _Y - _X
                    theta[i] += 1/_X * _X * volume[j] * extension  #* horizon_correction(X, particle_size, horizon)
                end
            end
            theta[i] *= 3 / m[i]
        end
        return nothing
    end

    kernel = CUDA.@cuda launch=false cal_theta(theta, y, x, intact, family, volume, m)
    config = launch_configuration(kernel.fun)
    nthreads = Base.min(length(volume), config.threads)
    nblocks =  cld(length(volume), nthreads)

    theta .= 0.0
    CUDA.@sync kernel(theta, y, x, intact, family, volume, m; threads=nthreads, blocks=nblocks)

    @check_nan theta "theta"

end


"""
    weighted_volume(x, volume, particle_size, family, horizon)

It gives weighted volume as given ordinary state material model.

``m_i = \\sum_{j \\in \\mathcal{N}_i} \\omega_{ij} \\times \\lVert X_{ij} \\rVert^2 \\times vol_j \\times s_{ij}``

where,

- ``\\omega_{ij}``: influence function
- ``X_{ij}``: ``x_j - x_i``
- ``vol_j``: particle volume
- ``s_{ij}``: horizon correction

# Arguments
- `x`: initial positions
- `volume`: particle volume
- `particle_size`: particle size
- `family`: family of particles
- `horizon`: horizon

# Returns
- `m`: weighted volume
"""
function weighted_volume(x, volume, particle_size, family, horizon)
    m = zero(volume) * unit(eltype(x))
    Threads.@threads for i in 1:size(x, 2)
        for k in 1:size(family, 1)
            j = family[k, i]
            if (j > 0) & (i != j)
                dr = get_ij(j, i, x)
                m[i] += influence_function(dr) * get_magnitude(dr)^2 * horizon_correction(dr, particle_size, horizon) * volume[j]
            end
        end
    end
    return m
end


"""
    loop_over_neighs!(y_, mat, fn_map, device; fn_reduce=(args...)->())

The function loops over neighbors and applies `fn_map` to each neighbor and `fn_reduce` to reduce the results.

The function `fn_map` should have the following signature:

```
fn_map(mat, i, j, k, type1, type2,
                                xij, yij, extension, s,
                                wij, wji, items)
```

where,

- `mat`: material model
- `i`: particle index
- `j`: neighbor index
- `k`: family index
- `type1`: type of particle `i`
- `type2`: type of particle `j`
- `xij`: `x_j - x_i`
- `yij`: `y_j - y_i`
- `extension`: `y_j - y_i`
- `s`: `s_ij` bond stretch
- `wij`: influence function `w_ij`
- `wji`: influence function `w_ji`
- `items`: output of `get_items(mat)` function

The function `fn_reduce` should have the following signature:

```
fn_reduce(mat, i, items)
```

where,

- `mat`: material model
- `i`: particle index
- `items`: output of `get_items(mat)` function


# Arguments
- `y_`: deformed positions
- `mat`: material model
- `fn_map`: function to be applied to each neighbor
- `device::Type{Val{:cuda}}`: device type

# Keyword Arguments
- `fn_reduce=(args...)->()`: function to reduce the results

# Returns
- `nothing`: `y_` is updated in-place
"""
function loop_over_neighs!(y_, mat, fn_map, device; fn_reduce=(args...)->())
    start = mat.type.start
    type = mat.general.type .- (start - 1)
    intact = mat.general.intact
    family = mat.general.family
    x = mat.general.x
    critical_stretch = mat.specific.critical_stretch
    items = get_items(mat)
    M, N = size(family)
    Threads.@threads for i in 1:N
        for k in 1:M
            if intact[k, i]
                j = family[k, i]
                X = get_ij(j, i, x)
                Y = get_ij(j, i, y_)
                yij = get_magnitude(Y)
                xij = get_magnitude(X)
                extension = yij - xij
                s = extension/xij
                type1 = type[i]
                type2 = type[j]
                if s < critical_stretch[type1, type2]
                    wij = influence_function(X)
                    wji = influence_function(.-X)
                    fn_map(mat, i, j, k, type1, type2,
                                xij, yij, extension, s,
                                wij, wji, items)
                else
                    intact[k, i] = false
                end
            end
        end
        fn_reduce(mat, i, items)
    end
end

"""
    neigh_cells(i, j, k, N)

Calculates neighbor cells for given i,j,k cell indices (when in a list).

# Arguments
- `i`: cell index in x direction
- `j`: cell index in y direction
- `k`: cell index in z direction
- `N`: number of cells in each direction

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
    cell_number(i, j, k, N)

Calculates cell number for given i,j,k cell indices (when in a list).

# Arguments
- `i`: cell index in x direction
- `j`: cell index in y direction
- `k`: cell index in z direction
- `N`: number of cells in each direction
"""
function cell_number(i, j, k, N)
    return i + (j - 1) * N[1] + (k - 1) * N[1] * N[2]
end

"""
    get_cells(x, horizon; max_part=30)

Calculates cells for given positions and horizon.

# Arguments
- `x`: positions
- `horizon`: horizon

# Keyword Arguments
- `max_part=30`: number of partitions in each direction
"""
function get_cells(x::AbstractArray{T,2}, horizon::T1; max_part=30) where {T,T1}
    x = Array(x)
    _min = Array(minimum(x, dims = 2))
    _max = Array(maximum(x, dims = 2))
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
    cal_family!(family, x, horizon)

Calculates family for given positions and horizon.

# Arguments
- `family`: family
- `x`: positions
- `horizon`: horizon
"""
function cal_family!(family::Array{Int64,2}, x::Array{T,2}, horizon::T1) where {T,T1}
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
    cal_family(x, horizon, max_neigh)

Calculates family for given positions and horizon.

# Arguments
- `x`: positions
- `horizon`: horizon
- `max_neigh`: maximum number of neighbors
"""
function cal_family(x::Array{T,2}, horizon::T1, max_neigh::Int64) where {T, T1}
    max_neigh = min(max_neigh, size(x, 2))
    family = zeros(Int64, max_neigh, size(x, 2))
    cal_family!(family, x, horizon)
    return family
end


"""
    cal_damage(env)

Calculates damage for each bond for a given environment.

# Arguments
- `env`: environment
"""
function cal_damage(env)
    intact = sum(env.material_blocks[1].general.intact, dims=1)
    for mat in env.material_blocks[2:end]
        intact = hcat(intact, sum(mat.general.intact, dims=1))
    end
    damage = 1 .- reshape(intact, :) ./ env.intact0
    return damage
end
