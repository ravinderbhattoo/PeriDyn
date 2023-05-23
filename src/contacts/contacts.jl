abstract type RepulsionModel11 end
abstract type RepulsionModel12 end

function _cudaconvert(x::Vector{T}) where T <: Union{RepulsionModel11,RepulsionModel12}
    _cudaconvert.(x)
end

function _cudaconvert(x::T) where T <: Union{RepulsionModel11,RepulsionModel12}
    function fn(x, k)
        if k!=:material
            return _cudaconvert(getfield(x, k))
        else
            return getfield(x, k)
        end
    end
    T((fn(x, k) for k in fieldnames(T))...)
end


@def RepulsionModel_gf begin
    equi_dist::Float64
    distance::Float64
    horizons::Tuple
    hor::Float64
    neighs::AbstractArray{Int64,2}
    max_neighs::Int64
end

@def RepulsionModel12_gf begin
    @RepulsionModel_gf
    pair::Vector{AbstractVector{Int64}}
end

@def RepulsionModel11_gf begin
    @RepulsionModel_gf
    type::AbstractVector{Int64}
    material::GeneralMaterial
end

function RepulsionModel12_gcal(mat1, mat2, distanceX, max_neighs)
    equi_dist = (mat1.general.particle_size + mat2.general.particle_size) / 2
    distance = distanceX * equi_dist
    horizons = (mat1.general.horizon, mat2.general.horizon)
    hor = sum(horizons) / 2
    neighs = zeros(min(max_neighs,size(mat2.general.x,2)), size(mat1.general.x,2))
    max_neighs = min(max_neighs,size(mat2.general.x,2))
    pair = [mat1.type,mat2.type]
    return (equi_dist, distance, horizons, hor, neighs, max_neighs, pair)
end

function RepulsionModel11_gcal(mat1, distanceX, max_neighs)
    equi_dist = mat1.general.particle_size
    distance = distanceX * equi_dist
    horizons = (mat1.general.horizon, )
    hor = mat1.general.horizon
    neighs = zeros(max_neighs, size(mat1.general.x,2))
    max_neighs = min(max_neighs,size(mat1.general.x,2))
    type = mat1.type
    material = mat1.general
    return (equi_dist, distance, horizons, hor, neighs, max_neighs, type, material)
end


function Base.show(io::IO, i::Union{RepulsionModel11, RepulsionModel12})
    println(io, typeof(i))
    for j in fieldnames(typeof(i))
        if j in [:neighs]
        else
        println(io, j, ": ", getproperty(i, j))
        end
    end
end

export short_range_repulsion!, collision_box, update_repulsive_neighs!


"""
    short_range_repulsion!(y,f,type,RepusionModel)

Updates (inplace) the repulsive acceleration of material points.

``α``

**`1-1 repulsion`**

# Input Args:
- `y :: Positions of material point`
- `f :: Acceleration of material points`
- `type :: Type of material points`
- `RepulsionModel :: Repulsion model (see contacts.jl for more details)`

# Output Args:
- `Noting (Inplace updation of f (acceleration))`

# Examples
```jldoctest
julia> short_range_repulsion!(y,f,type,RepusionModel)
```
"""
function short_range_repulsion!(y, f, type, bid, vol, RM::RepulsionModel12)
    mask1 = false
    for i in RM.pair[1]
        mask1 = mask1 .| (type .== i)
    end
    mask2 = false
    for i in RM.pair[2]
        mask2 = mask2 .| (type .== i)
    end
    f1 = Array(f[:, mask1])
    f2 = Array(f[:, mask2])
    x1 = Array(y[:, mask1])
    x2 = Array(y[:, mask2])
    vol1 = Array(vol[mask1])
    vol2 = Array(vol[mask2])
    neighs = Array(RM.neighs)
    for i in 1:size(neighs, 2)
        for k in 1:size(neighs, 1)
            j = neighs[k,i]
            if j>0
                force = repulsion_force(x1[:, i] .- x2[:, j], RM)
                f1[:,i] .+= force * vol2[j]
                f2[:,j] .-=  force * vol1[i]
            end
            if j==0 break end
        end
    end
    f[:, mask1] .= typeof(f)(f1)
    f[:, mask2] .= typeof(f)(f2)
    return nothing
end



function short_range_repulsion!(y, f, type, bid, vol, RM::RepulsionModel11)
    device = DEVICE[]
    short_range_repulsion!(y, f, type, bid, vol, RM::RepulsionModel11, Val{device})
end

"""
    short_range_repulsion!(y,f,type,RepusionModel)

Updates (inplace) the repulsive acceleration of material points.
"""
function short_range_repulsion!(y, f, type, bid, vol, RM::RepulsionModel11, device::Type{Val{:cuda}})
    force_fn = get_repulsion_force_fn(RM)
    mask1 = bid .== RM.bid

    type = type[mask1]
    f1 = f[:,mask1]
    x1 = y[:,mask1]
    vol1 = vol[mask1]
    neighs = RM.neighs

    function cal_force(x1, vol1, neighs, type, f1)
        index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        stride = blockDim().x * gridDim().x

        for i in index:stride:size(neighs, 2)
            for k in 1:size(neighs, 1)
                j = neighs[k,i]
                if (j>0) && (type[i]!=-1) && (type[j]!=-1)
                    force = force_fn((
                        x1[1, i]-x1[1, j],
                        x1[2, i]-x1[2, j],
                        x1[3, i]-x1[3, j]
                    ))

                    f1[1,i] += force[1] / vol1[j]
                    f1[2,i] += force[2] / vol1[j]
                    f1[3,i] += force[3] / vol1[j]

                    f1[1,j] -= force[1] / vol1[i]
                    f1[2,j] -= force[2] / vol1[i]
                    f1[3,j] -= force[3] / vol1[i]
                end
            end
        end
        return nothing
    end

    kernel = CUDA.@cuda launch=false cal_force(x1, vol1, neighs, type, f1)
    config = launch_configuration(kernel.fun)
    nthreads = Base.min(length(vol1), config.threads)
    nblocks =  cld(length(vol1), nthreads)

    CUDA.@sync kernel(x1, vol1, neighs, type, f1; threads=nthreads, blocks=nblocks)
    # f1 = apply_kernel(cal_force, x1, vol1, neighs, type, f1)[end]
    f[:, mask1] .= typeof(f)(f1)
    return nothing
end


"""
    short_range_repulsion!(y,f,type,RepusionModel)

Updates (inplace) the repulsive acceleration of material points.
"""
function short_range_repulsion!(y, f, type, bid, vol, RM::RepulsionModel11, device::Type{Val{:cpu}})
    mask1 = bid .== RM.bid
    f1 = f[:, mask1]
    x1 = y[:, mask1]
    vol1 = vol[mask1]
    for i in 1:size(RM.neighs, 2)
        for k in 1:size(RM.neighs, 1)
            j = RM.neighs[k,i]
            if j>0
                force = repulsion_force(x1[:, i] .- x2[:, j], RM)
                f1[:,i] .+= force * vol2[j]
                f2[:,j] .-=  force * vol1[i]
            end
            if j==0 break end
        end
    end
    f[:, mask1] .= f1
    return nothing
end


"""
    collision_box(x1::Array{Float64,2}, x2::Array{Float64,2}, skin::Float64)

Calculates collsion box between two material blocks.

# Input Args:
- `x1 :: Positions of material point (block 1)`
- `x2 :: Positions of material point (block 2)`
- `skin :: Extra distance need to consider (usually >= particle size)`

# Output Args:
- `box_min :: Minimum position limits for overlap`
- `box_max :: Miximum position limits for overlap`
- `ifoverlap :: Boolean (true if overlap)`

# Examples
```jldoctest
julia> collision_box(x1, x2, skin)
```
"""
function collision_box(x1::AbstractArray, x2::AbstractArray, skin::Float64)
    min1 = minimum(x1,dims=2).-skin
    min2 = minimum(x2,dims=2).-skin
    max1 = maximum(x1,dims=2).+skin
    max2 = maximum(x2,dims=2).+skin
    box_min = max.(min1,min2)
    box_max = min.(max1,max2)
    vol = prod(box_max-box_min)
    if vol<0.0
        return box_min, box_max, false
    else
        return box_min, box_max, true
    end
end


"""
    update_repulsive_neighs!(y,type,RepulsionModel12)

Update neighbour list for repulsive force calculation (1-2 interaction).

"""
function update_repulsive_neighs!(y, type, RM::RepulsionModel12; max_part=nothing)
    log_info("Updating repulsive neighs ...")
    _mask1 = false
    for i in RM.pair[1]
        _mask1 = _mask1 .| (type .== i)
    end

    _mask2 = false
    for i in RM.pair[2]
        _mask2 = _mask2 .| (type .== i)
    end

    x11 = Array(y[:, _mask1])
    x22 = Array(y[:, _mask2])

    RM.neighs .= 0

    box_min, box_max, ifcheck = collision_box(x11, x22, RM.distance)
    if ifcheck
        mask1 = reshape(prod(x11.>box_min,dims=1) .* prod(x11.<box_max,dims=1),:)
        mask2 = reshape(prod(x22.>box_min,dims=1) .* prod(x22.<box_max,dims=1),:)
        x1 = x11[1:end, mask1]
        x2 = x22[1:end, mask2]
        a_id = collect(1:size(x22,2))[mask2]
        family = zeros(Float64, max(RM.max_neighs, size(x2,2)), size(x1,2))
        for i in 1:size(x1,2)
            a1,b1,c1 = x1[1,i],x1[2,i],x1[3,i]
            for j in 1:size(x2,2)
                a2,b2,c2 = x2[1,j],x2[2,j],x2[3,j]
                if (a1-a2)*(a1-a2)+(b1-b2)*(b1-b2)+(c1-c2)*(c1-c2)<RM.distance^2
                    family[j,i] = a_id[j]
                end
            end
        end
        family = sort(family,dims=1)
        RM.neighs[1:end, mask1] = family[end:-1:end+1-RM.max_neighs,1:end]
        log_info("Average repulsive neighs: $(sum(RM.neighs .> 0.5) / (1+size(x1, 2)))")
    else
        RM.neighs[1:end, 1:end] .= 0
        log_info("No repulsive neighs.")
    end
    log_info("Done")
end

function update_repulsive_neighs!(y, type, RM::RepulsionModel11; kwargs...)
    device = DEVICE[]
    update_repulsive_neighs!(y, type, RM::RepulsionModel11, device; kwargs...)
end

function update_repulsive_neighs!(y, type, RM::RepulsionModel11, device::Symbol; kwargs...)
    update_repulsive_neighs!(y, type, RM::RepulsionModel11, Val{device}; kwargs...)
end



function update_repulsive_neighs!(neighbors, x, search_distance, equi_dist, family, intact, max_part)
    cells, cell_neighs = get_cells(x, search_distance; max_part=max_part)
    Threads.@threads for cell_i in 1:length(cells)
        for ca_id in 1:length(cells[cell_i])
            ind = 1
            ca = cells[cell_i][ca_id]
            a1,b1,c1 = x[1,ca],x[2,ca],x[3,ca]
            for neigh_id in 1:length(cell_neighs[cell_i])
                neighs = cells[cell_neighs[cell_i][neigh_id]]
                for fa_id in 1:length(neighs)
                    fa = neighs[fa_id]
                    if fa != ca
                        ifbroken = true
                        # Check if family is intact
                        for k in size(family, 1):-1:1
                            if family[k, ca]==0
                                break
                            else ((fa == family[k, ca]) && intact[k, ca]==1)
                                ifbroken = false
                                break
                            end
                        end
                        if ifbroken
                            a2,b2,c2 = x[1,fa], x[2,fa], x[3,fa]
                            dr2 = ((a1-a2)*(a1-a2)+(b1-b2)*(b1-b2)+(c1-c2)*(c1-c2))
                            if dr2 < equi_dist^2
                                neighbors[ind, ca] = fa
                                ind += 1
                            end
                        end
                    end
                end
            end
        end
    end
    log_info("Average repulsive neighs: $(sum(neighbors .> 0.5)/size(x, 2))")
end

"""
    update_repulsive_neighs!(y,type,RepulsionModel11)

Update neighbour list for repulsive force calculation (1-1 interaction).

"""
function update_repulsive_neighs!(y, type, RM::RepulsionModel11, device::Type{Val{:cpu}}; max_part=30)
    mask = false
    for j in RM.type
        mask = mask .| (type .== j)
    end
    x = y[:, mask]
    neighbors = RM.neighs
    search_distance = RM.distance
    equi_dist = RM.material.particle_size
    family = RM.material.family
    intact = RM.material.intact
    fill!(neighbors,0)
    update_repulsive_neighs!(neighbors, x, search_distance, equi_dist, family, intact, max_part)
end

"""
    update_repulsive_neighs!(y,type,RepulsionModel11)

Update neighbour list for repulsive force calculation (1-1 interaction).

"""
function update_repulsive_neighs!(y, type, RM::RepulsionModel11, device::Type{Val{:cuda}}; max_part=30)
    mask = false
    for j in RM.type
        mask = mask .| (type .== j)
    end
    x = Array(y[:, mask])
    neighbors = Array(RM.neighs)
    search_distance = RM.distance
    equi_dist = RM.material.particle_size
    family = Array(RM.material.family)
    intact = Array(RM.material.intact)
    fill!(neighbors,0)
    update_repulsive_neighs!(neighbors, x, search_distance, equi_dist, family, intact, max_part)
    RM.neighs .= CuArray(neighbors)
    nothing
end

include("./LJ.jl")
include("./shortrange.jl")
include("./nonlinear.jl")
include("./linear.jl")


#
