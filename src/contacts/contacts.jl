"""
This module contains functions for calculating the repulsive forces between material points.
"""

export RepulsionModel11, RepulsionModel12
export RepulsionModel11_gcal, RepulsionModel12_gcal
export short_range_repulsion!, collision_box, update_repulsive_neighs!

"""
    RepulsionModel11

Abstract type for repulsion model for a single material type.
"""
abstract type RepulsionModel11 end


"""
    RepulsionModel12

Abstract type for repulsion model for two material types.
"""
abstract type RepulsionModel12 end

@def RepulsionModel_gf begin
    name::String
    equi_dist::Float64
    distance::Float64
    neighs::AbstractArray{Int64,2}
    max_neighs::Int64
end

@def RepulsionModel12_gf begin
    pair::Vector{AbstractVector{Int64}}
    # name::String
    # equi_dist::Float64
    # distance::Float64
    # neighs::AbstractArray{Int64,2}
    # max_neighs::Int64
    @RepulsionModel_gf
end

@def RepulsionModel11_gf begin
    type::AbstractVector{Int64}
    bid::Int64
    material::GeneralMaterial
    # name::String
    # equi_dist::Float64
    # distance::Float64
    # neighs::AbstractArray{Int64,2}
    # max_neighs::Int64
    @RepulsionModel_gf
end

"""
    RepulsionModel12_gcal(mat1, mat2, distanceX, max_neighs)

Calculate the parameters for `RepulsionModel12` based on the input materials, distance, and maximum neighbors.

Arguments:
- `mat1`: The first material (type `RepulsionModel11`).
- `mat2`: The second material (type `RepulsionModel11`).
- `distanceX`: The distance factor.
- `max_neighs`: The maximum number of neighbors.

Returns a tuple containing the calculated parameters for `RepulsionModel12`.
"""
function RepulsionModel12_gcal(mat1, mat2, distanceD, distanceX, max_neighs)
    pair = [mat1.type, mat2.type]
    name = "$(mat1.name)-$(mat2.name)"
    equi_dist = (mat1.general.particle_size + mat2.general.particle_size) / 2 * distanceD
    distance = distanceX * equi_dist
    neighs = zeros(min(max_neighs,size(mat2.general.x,2)), size(mat1.general.x,2))
    max_neighs = min(max_neighs,size(mat2.general.x,2))
    # Order of arguments is important
    # pair::Vector{AbstractVector{Int64}}
    # name::String
    # equi_dist::Float64
    # distance::Float64
    # neighs::AbstractArray{Int64,2}
    # max_neighs::Int64
    return pair, name, equi_dist, distance, neighs, max_neighs
end

"""
    RepulsionModel11_gcal(mat1, distanceX, max_neighs)

Calculate the parameters for `RepulsionModel11` based on the input material, distance, and maximum neighbors.

Arguments:
- `mat1`: The material (type `RepulsionModel11`).
- `distanceX`: The distance factor.
- `max_neighs`: The maximum number of neighbors.

Returns a tuple containing the calculated parameters for `RepulsionModel11`.
"""
function RepulsionModel11_gcal(mat1, distanceD, distanceX, max_neighs)
    type = mat1.type
    bid = mat1.blockid
    material = mat1.general
    name = "$(mat1.name)-$(mat1.name)"
    equi_dist = mat1.general.particle_size * distanceD
    distance = distanceX * equi_dist
    neighs = zeros(max_neighs, size(mat1.general.x,2))
    max_neighs = min(max_neighs,size(mat1.general.x,2))
    # Order of arguments is important
    # type::AbstractVector{Int64}
    # bid::Int64
    # material::GeneralMaterial
    # name::String
    # equi_dist::Float64
    # distance::Float64
    # neighs::AbstractArray{Int64,2}
    # max_neighs::Int64
    return type, bid, material, name, equi_dist, distance, neighs, max_neighs
end


function Base.show(io::IO, i::Union{RepulsionModel11, RepulsionModel12})
    """
    Print the details of a RepulsionModel object.

    # Arguments
    - `io::IO`: The output IO stream.
    - `i::Union{RepulsionModel11, RepulsionModel12}`: The RepulsionModel object to print.
    """
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
    short_range_repulsion!(y, f, type, RepulsionModel)

Updates (inplace) the repulsive acceleration of material points.

## Arguments
- `y`: Positions of material points.
- `f`: Acceleration of material points.
- `type`: Type of material points.
- `RepulsionModel`: Repulsion model (see contacts.jl for more details).

## Output
- None (Inplace update of f (acceleration)).
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


"""
    short_range_repulsion!(y, f, type, bid, vol, RM::RepulsionModel11)

Updates (inplace) the repulsive acceleration of material points.

# Arguments
- `y`: Positions of material points
- `f`: Acceleration of material points
- `type`: Type of material points
- `bid`: BID values
- `vol`: Volume values
- `RM::RepulsionModel11`: Repulsion model

# Output
- No return value. The function updates `f` in place.

"""
function short_range_repulsion!(y, f, type, bid, vol, RM::RepulsionModel11)
    device = DEVICE[]
    short_range_repulsion!(y, f, type, bid, vol, RM::RepulsionModel11, Val{device})
end

"""
    short_range_repulsion!(y, f, type, bid, vol, RM::RepulsionModel11, device::Type{Val{:cuda}})

Updates (inplace) the repulsive acceleration of material points.

# Arguments
- `y`: Positions of material points
- `f`: Acceleration of material points
- `type`: Type of material points
- `bid`: BID values
- `vol`: Volume values
- `RM::RepulsionModel11`: Repulsion model
- `device::Type{Val{:cuda}}`: Device type for CUDA acceleration

# Output
- No return value. The function updates `f` in place.

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
    short_range_repulsion!(y, f, type, bid, vol, RM::RepulsionModel11, device::Type{Val{:cpu}})

Updates (inplace) the repulsive acceleration of material points.

# Arguments
- `y`: Positions of material points
- `f`: Acceleration of material points
- `type`: Type of material points
- `bid`: BID values
- `vol`: Volume values
- `RM::RepulsionModel11`: Repulsion model
- `device::Type{Val{:cpu}}`: Device type for CPU acceleration

# Output
- No return value. The function updates `f` in place.

"""
function short_range_repulsion!(y, f, type, bid, vol, RM::RepulsionModel11, device::Type{Val{:cpu}})
    force_fn = get_repulsion_force_fn(RM)
    mask1 = bid .== RM.bid
    f1 = @view f[:, mask1]
    x1 = @view y[:, mask1]
    vol1 = @view vol[mask1]
    for i in 1:size(RM.neighs, 2)
        Threads.@threads for k in 1:size(RM.neighs, 1)
            j = RM.neighs[k,i]
            if j != 0
                force = force_fn(x1[:, i] .- x1[:, j])
                f1[:,i] .+= force * vol1[j]
                f1[:,j] .-=  force * vol1[i]
            end
        end
    end
    # f[:, mask1] .= f1
    return nothing
end


"""
    collision_box(x1::Array{Float64,2}, x2::Array{Float64,2}, skin::Float64)

Calculates collision box between two material blocks.

# Arguments
- `x1`: Positions of material points (block 1)
- `x2`: Positions of material points (block 2)
- `skin`: Extra distance to consider (usually >= particle size)

# Output
- `box_min`: Minimum position limits for overlap
- `box_max`: Maximum position limits for overlap
- `ifoverlap`: Boolean indicating if there is an overlap

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
    update_repulsive_neighs!(y, type, RM::RepulsionModel12; max_part=nothing)

Update neighbor list for repulsive force calculation (1-2 interaction).

# Arguments
- `y`: Positions of material points
- `type`: Type of material points
- `RM::RepulsionModel12`: Repulsion model
- `max_part=nothing`: Maximum number of particles (optional)

# Output
- No return value. The function updates `RM.neighs` in place.

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


"""
    update_repulsive_neighs!(y, type, RM::RepulsionModel11; kwargs...)

Update neighbor list for repulsive force calculation (1-1 interaction).

# Arguments
- `y`: Positions of material points
- `type`: Type of material points
- `RM::RepulsionModel11`: Repulsion model
- `kwargs...`: Additional keyword arguments

# Output
- No return value. The function updates `RM.neighs` in place.

"""
function update_repulsive_neighs!(y, type, RM::RepulsionModel11; kwargs...)
    device = DEVICE[]
    update_repulsive_neighs!(y, type, RM::RepulsionModel11, device; kwargs...)
end


"""
    update_repulsive_neighs!(y, type, RM::RepulsionModel11, device::Symbol; kwargs...)

Update neighbor list for repulsive force calculation (1-1 interaction) on a specific device.

# Arguments
- `y`: Positions of material points
- `type`: Type of material points
- `RM::RepulsionModel11`: Repulsion model
- `device::Symbol`: Device type for acceleration
- `kwargs...`: Additional keyword arguments

# Output
- No return value. The function updates `RM.neighs` in place.

"""
function update_repulsive_neighs!(y, type, RM::RepulsionModel11, device::Symbol; kwargs...)
    update_repulsive_neighs!(y, type, RM::RepulsionModel11, Val{device}; kwargs...)
end


"""
    update_repulsive_neighs!(neighbors, x, search_distance, equi_dist, family, intact, max_part)

Update the neighbor list for repulsive force calculation.

# Arguments
- `neighbors`: Array storing the neighbor indices
- `x`: Positions of material points
- `search_distance`: Maximum search distance for neighbors
- `equi_dist`: Equilibrium distance for repulsion
- `family`: Array indicating the family relationship between material points
- `intact`: Array indicating if the family relationship is intact
- `max_part`: Maximum number of particles in a cell

# Output
- No return value. The function updates the `neighbors` array in place.

"""
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
                            if (fa == family[k, ca])
                                if intact[k, ca] == 1
                                    ifbroken = false
                                end
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
    update_repulsive_neighs!(y, type, RM::RepulsionModel11, device::Type{Val{:cpu}}; max_part=30)

Update the neighbor list for repulsive force calculation (1-1 interaction).

# Arguments
- `y`: Positions of material points
- `type`: Type of material points
- `RM::RepulsionModel11`: Repulsion model for 1-1 interaction
- `device::Type{Val{:cpu}}`: Device type (CPU)
- `max_part=30`: Maximum number of particles in a cell

# Output
- No return value. The function updates the neighbor list in the `RM` object.

"""
function update_repulsive_neighs!(y, type, RM::RepulsionModel11, device::Type{Val{:cpu}}; max_part=30)
    mask = false
    for j in RM.type
        mask = mask .| (type .== j)
    end
    x = y[:, mask]
    neighbors = RM.neighs
    search_distance = RM.distance
    equi_dist = RM.equi_dist
    family = RM.material.family
    intact = RM.material.intact
    fill!(neighbors,0)
    update_repulsive_neighs!(neighbors, x, search_distance, equi_dist, family, intact, max_part)
end

"""
    update_repulsive_neighs!(y, type, RM::RepulsionModel11, device::Type{Val{:cuda}}; max_part=30)

Update neighbor list for repulsive force calculation (1-1 interaction).

# Arguments
- `y`: Positions of material points
- `type`: Type of material points
- `RM::RepulsionModel11`: Repulsion model for 1-1 interaction
- `device::Type{Val{:cuda}}`: Device type (CUDA)
- `max_part=30`: Maximum number of particles in a cell

# Output
- No return value. The function updates the neighbor list in the `RM` object.

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


"""
    _cudaconvert(x::Vector{T}) where T <: Union{RepulsionModel11,RepulsionModel12}

Converts a vector of `RepulsionModel11` or `RepulsionModel12` objects to CUDA-compatible types.
"""
function _cudaconvert(x::Vector{T}) where T <: Union{RepulsionModel11,RepulsionModel12}
    _cudaconvert.(x)
end


"""
    _cudaconvert(x::T) where T <: Union{RepulsionModel11,RepulsionModel12}

Converts a single `RepulsionModel11` or `RepulsionModel12` object to a CUDA-compatible type.
"""
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


include("./nonlinearspring.jl")
include("./linearspring.jl")
include("./shortrange.jl")
include("./LJ.jl")


#
