"""
This module contains functions for calculating the contact forces between material points.
"""

export ContactModel, ContactModel11, ContactModel12
export ContactModel11_gcal, ContactModel12_gcal
export short_range_contact!
export collision_box, update_contact_neighs!


"""
    ContactModel

Abstract type for contact model.
"""
abstract type ContactModel end

"""
    ContactModel11

Abstract type for contact model for a single material block.
"""
abstract type ContactModel11 <: ContactModel end


"""
    ContactModel12

Abstract type for contact model between two material blocks.
"""
abstract type ContactModel12 <: ContactModel end

@def ContactModel_gf begin
    name::String
    equi_dist::QF
    distance::QF
    neighs::AbstractArray{Int64,2}
    max_neighs::Int64
end

@def ContactModel12_gf begin
    pair::Pair{AbstractVector{Int64}, AbstractVector{Int64}}
    # name::String
    # equi_dist::QF
    # distance::QF
    # neighs::AbstractArray{Int64,2}
    # max_neighs::Int64
    @ContactModel_gf
end

@def ContactModel11_gf begin
    type::AbstractVector{Int64}
    bid::Int64
    material::GeneralMaterial
    # name::String
    # equi_dist::QF
    # distance::QF
    # neighs::AbstractArray{Int64,2}
    # max_neighs::Int64
    @ContactModel_gf
end

"""
    ContactModel12_gcal(mat1, mat2, distanceX, max_neighs)

Calculate the parameters for `ContactModel12` based on the input materials, distance, and maximum neighbors.

Arguments:
- `mat1`: The first material (type `ContactModel11`).
- `mat2`: The second material (type `ContactModel11`).
- `distanceX`: The distance factor.
- `max_neighs`: The maximum number of neighbors.

Returns a tuple containing the calculated parameters for `ContactModel12`.
"""
function ContactModel12_gcal(mat1, mat2, distanceD, distanceX, max_neighs)
    pair = Pair(mat1.type, mat2.type)
    name = "$(mat1.name)-$(mat2.name)"
    equi_dist = (mat1.general.particle_size + mat2.general.particle_size) / 2 * distanceD
    distance = distanceX * equi_dist
    neighs = zeros(min(max_neighs,size(mat2.general.x,2)), size(mat1.general.x,2))
    max_neighs = min(max_neighs,size(mat2.general.x,2))
    # Order of arguments is important
    # pair::Vector{AbstractVector{Int64}}
    # name::String
    # equi_dist::QF
    # distance::QF
    # neighs::AbstractArray{Int64,2}
    # max_neighs::Int64
    return pair, name, equi_dist, distance, neighs, max_neighs
end

"""
    ContactModel11_gcal(mat1, distanceX, max_neighs)

Calculate the parameters for `ContactModel11` based on the input material, distance, and maximum neighbors.

Arguments:
- `mat1`: The material (type `ContactModel11`).
- `distanceX`: The distance factor.
- `max_neighs`: The maximum number of neighbors.

Returns a tuple containing the calculated parameters for `ContactModel11`.
"""
function ContactModel11_gcal(mat1, distanceD, distanceX, max_neighs)
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
    # equi_dist::QF
    # distance::QF
    # neighs::AbstractArray{Int64,2}
    # max_neighs::Int64
    return type, bid, material, name, equi_dist, distance, neighs, max_neighs
end


"""
Print the details of a ContactModel object.

# Arguments
- `io::IO`: The output IO stream.
- `i::Union{ContactModel11, ContactModel12}`: The ContactModel object to print.
"""
function Base.show(io::IO, i::Union{ContactModel11, ContactModel12})
    print(io, getPanel(i))
end

function getPanel(i::ContactModel; ptype=SPanel, width=Term.default_width())
    txt = @bold(@blue("$(typeof(i))")) * "\n"
    for j in fieldnames(typeof(i))
        item = getproperty(i, j)
        txt = txt * "$(variable_color(j)): "
        if j in [:neighs]
            txt = txt * "$(size(getproperty(i, j)))" * "\n"
        elseif j==:material
            txt = txt * "\n$(getPanel(getproperty(i, j); ptype=DPanel, width=width-6))"
        else
            # check if iterable
            if isa(item, AbstractArray)
                txt = txt * array_repr(item)
            else
                txt = txt * "$(item)\n"
            end
        end
    end
    txt = txt[1:end-1]
    return ptype(txt,
        title=i.name,
        justify=:left,
        # subtitle=@blue("$(typeof(i))"),
        # subtitle_justify=:center,
        width=width)
end


###############################################
# short range contact
###############################################

export short_range_contact!, collision_box, update_contact_neighs!

"""
    short_range_contact!(y, f, type, bid, vol, RM::ContactModel12)

Updates (inplace) the contact acceleration of material points.

# Arguments
- `y`: Positions of material points
- `f`: Acceleration of material points
- `type`: Type of material points
- `bid`: BID values
- `vol`: Volume values
- `RM::ContactModel12`: Repulsion model

# Output
- No return value. The function updates `f` in place.

"""
function short_range_contact!(y, f, type, bid, vol, den, RM::ContactModel12)
    device = DEVICE[]
    device = :cpu
    log_info("short_range_contact! 12 used ($device).")
    short_range_contact!(y, f, type, bid, vol, den, RM::ContactModel12, Val{device})
end

"""
    short_range_contact!(y, f, type, ContactModel)

Updates (inplace) the contact acceleration of material points.

## Arguments
- `y`: Positions of material points.
- `f`: Acceleration of material points.
- `type`: Type of material points.
- `ContactModel`: Repulsion model (see contacts.jl for more details).

## Output
- None (Inplace update of f (acceleration)).
"""
function short_range_contact!(y, f, type, bid, vol, den, RM::ContactModel12, device::Type{Val{:cpu}})
    force_fn = get_contact_force_fn(RM)
    mask1 = false
    for i in RM.pair[1]
        mask1 = mask1 .| (type .== i)
    end
    mask2 = false
    for i in RM.pair[2]
        mask2 = mask2 .| (type .== i)
    end
    f1 = @view f[:, mask1]
    f2 = @view f[:, mask2]
    x1 = @view y[:, mask1]
    x2 = @view y[:, mask2]
    vol1 = @view vol[mask1]
    vol2 = @view vol[mask2]
    den1 = @view den[mask1]
    den2 = @view den[mask2]
    neighs = RM.neighs
    M, N = size(neighs)
    Threads.@threads for i in 1:N
        _den1_i = 1/den1[i]
        for k in 1:M
            j = neighs[k, i]
            if j>0
                _den2_j = 1/den2[j]
                dr = (x1[1, i] - x2[1, j],
                        x1[2, i] - x2[2, j],
                        x1[3, i] - x2[3, j])
                force = force_fn(dr)
                f1[:,i] .+= force .* (vol2[j] * _den1_i)
                f2[:,j] .-=  force .* (vol1[i] * _den2_j)
            end
            if j==0 break end
        end
    end
    return nothing
end


"""
    short_range_contact!(y, f, type, bid, vol, RM::ContactModel11)

Updates (inplace) the contact acceleration of material points.

# Arguments
- `y`: Positions of material points
- `f`: Acceleration of material points
- `type`: Type of material points
- `bid`: BID values
- `vol`: Volume values
- `RM::ContactModel11`: Repulsion model

# Output
- No return value. The function updates `f` in place.

"""
function short_range_contact!(y, f, type, bid, vol, den, RM::ContactModel11)
    device = DEVICE[]
    short_range_contact!(y, f, type, bid, vol, den, RM::ContactModel11, Val{device})
end

"""
    short_range_contact!(y, f, type, bid, vol, RM::ContactModel11, device::Type{Val{:cuda}})

Updates (inplace) the contact acceleration of material points.

# Arguments
- `y`: Positions of material points
- `f`: Acceleration of material points
- `type`: Type of material points
- `bid`: BID values
- `vol`: Volume values
- `RM::ContactModel11`: Repulsion model
- `device::Type{Val{:cuda}}`: Device type for CUDA acceleration

# Output
- No return value. The function updates `f` in place.

"""
function short_range_contact!(y, f, type, bid, vol, den, RM::ContactModel11, device::Type{Val{:cuda}})
    force_fn = get_contact_force_fn(RM)
    mask1 = bid .== RM.bid

    type = @view type[mask1]
    f1 = @view f[:, mask1]
    x1 = @view y[:, mask1]
    vol1 = @view vol[mask1]
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
    # f[:, mask1] .= typeof(f)(f1)
    return nothing
end


"""
    short_range_contact!(y, f, type, bid, vol, RM::ContactModel11, device::Type{Val{:cpu}})

Updates (inplace) the contact acceleration of material points.

# Arguments
- `y`: Positions of material points
- `f`: Acceleration of material points
- `type`: Type of material points
- `bid`: BID values
- `vol`: Volume values
- `RM::ContactModel11`: Repulsion model
- `device::Type{Val{:cpu}}`: Device type for CPU acceleration

# Output
- No return value. The function updates `f` in place.

"""
function short_range_contact!(y, f, type, bid, vol, den, RM::ContactModel11, device::Type{Val{:cpu}})
    force_fn = get_contact_force_fn(RM)
    mask1 = bid .== RM.bid
    f1 = @view f[:, mask1]
    x1 = @view y[:, mask1]
    vol1 = @view vol[mask1]
    den1 = @view den[mask1]
    M, N = size(RM.neighs)
    neighs = RM.neighs
    Threads.@threads for i in 1:N
        _den_i = den1[i]
        for k in 1:M
            j = neighs[k,i]
            if j != 0
                _den_j = den1[j]
                dr = get_ij(i, j, x1)
                force = force_fn(dr)
                # (force_density[force/vol/vol] * volume)
                f1[:, i] .+= force .* (vol1[j] / _den_i)
                f1[:, j] .-=  force .* (vol1[j] / _den_j)
            end
        end
    end
    return nothing
end

###############################################
# update contact neibours
###############################################

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
function collision_box(x1::AbstractArray, x2::AbstractArray, skin::T) where T
    min1 = minimum(x1,dims=2).-skin
    min2 = minimum(x2,dims=2).-skin
    max1 = maximum(x1,dims=2).+skin
    max2 = maximum(x2,dims=2).+skin
    box_min = max.(min1,min2)
    box_max = min.(max1,max2)
    vol = prod(box_max-box_min)
    if vol < zero(vol)
        return box_min, box_max, false
    else
        return box_min, box_max, true
    end
end


"""
    update_contact_neighs!(y, type, RM::ContactModel12; max_part=nothing)

Update neighbor list for contact force calculation (1-2 interaction).

# Arguments
- `y`: Positions of material points
- `type`: Type of material points
- `RM::ContactModel12`: Repulsion model
- `max_part=nothing`: Maximum number of particles (optional)

# Output
- No return value. The function updates `RM.neighs` in place.

"""
function update_contact_neighs!(y, type, RM::ContactModel12; max_part=nothing)
    log_info("Updating contact neighs ...")
    _mask1 = false
    for i in RM.pair[1]
        _mask1 = _mask1 .| (type .== i)
    end

    _mask2 = false
    for i in RM.pair[2]
        _mask2 = _mask2 .| (type .== i)
    end

    x11 = @view y[:, _mask1]
    x22 = @view y[:, _mask2]

    RM.neighs .= 0

    box_min, box_max, ifcheck = collision_box(x11, x22, RM.distance)
    if ifcheck
        mask1 = reshape(prod(x11.>box_min,dims=1) .* prod(x11.<box_max,dims=1),:)
        mask2 = reshape(prod(x22.>box_min,dims=1) .* prod(x22.<box_max,dims=1),:)
        x1 = @view x11[1:end, mask1]
        x2 = @view x22[1:end, mask2]
        a_id = collect(1:size(x22,2))[mask2]
        family = zeros(Int64, max(RM.max_neighs, size(x2,2)), size(x1,2))
        Threads.@threads for i in 1:size(x1,2)
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
        log_info("Average contact neighs: $(sum(RM.neighs .> 0.5) / (1+size(x1, 2)))")
    else
        RM.neighs[1:end, 1:end] .= 0
        log_info("No contact neighs.")
    end
    log_info("Done")
end


"""
    update_contact_neighs!(y, type, RM::ContactModel11; kwargs...)

Update neighbor list for contact force calculation (1-1 interaction).

# Arguments
- `y`: Positions of material points
- `type`: Type of material points
- `RM::ContactModel11`: Repulsion model
- `kwargs...`: Additional keyword arguments

# Output
- No return value. The function updates `RM.neighs` in place.

"""
function update_contact_neighs!(y, type, RM::ContactModel11; kwargs...)
    device = DEVICE[]
    update_contact_neighs!(y, type, RM::ContactModel11, device; kwargs...)
end


"""
    update_contact_neighs!(y, type, RM::ContactModel11, device::Symbol; kwargs...)

Update neighbor list for contact force calculation (1-1 interaction) on a specific device.

# Arguments
- `y`: Positions of material points
- `type`: Type of material points
- `RM::ContactModel11`: Repulsion model
- `device::Symbol`: Device type for acceleration
- `kwargs...`: Additional keyword arguments

# Output
- No return value. The function updates `RM.neighs` in place.

"""
function update_contact_neighs!(y, type, RM::ContactModel11, device::Symbol; kwargs...)
    update_contact_neighs!(y, type, RM::ContactModel11, Val{device}; kwargs...)
end


"""
    update_contact_neighs!(neighbors, x, search_distance, equi_dist, family, intact, max_part)

Update the neighbor list for contact force calculation.

# Arguments
- `neighbors`: Array storing the neighbor indices
- `x`: Positions of material points
- `search_distance`: Maximum search distance for neighbors
- `equi_dist`: Equilibrium distance for contact
- `family`: Array indicating the family relationship between material points
- `intact`: Array indicating if the family relationship is intact
- `max_part`: Maximum number of particles in a cell

# Output
- No return value. The function updates the `neighbors` array in place.

"""
function update_contact_neighs!(mask, neighbors, x, search_distance, equi_dist, family, intact, max_part)
    cells, cell_neighs = get_cells(x, search_distance; max_part=max_part)
    FS = size(family, 1)
    equi_dist2 = equi_dist^2
    for cell_i in 1:length(cells)
        Threads.@threads for ca_id in 1:length(cells[cell_i])
            ind = 1
            ca = cells[cell_i][ca_id]
            if mask[ca]
                a1,b1,c1 = x[1,ca],x[2,ca],x[3,ca]
                for neigh_id in 1:length(cell_neighs[cell_i])
                    neighs = sort!(cells[cell_neighs[cell_i][neigh_id]])
                    fms = @view family[:, ca]
                    fm_ind = 1
                    for neigh in neighs
                        while (fms[fm_ind] < neigh) && (fm_ind < FS)
                            fm_ind += 1
                        end
                        if fms[fm_ind] == neigh
                            if ~intact[fm_ind, ca]
                                a2,b2,c2 = x[1,neigh], x[2,neigh], x[3,neigh]
                                dr2 = ((a1-a2)*(a1-a2)+(b1-b2)*(b1-b2)+(c1-c2)*(c1-c2))
                                if dr2 < equi_dist2
                                    neighbors[ind, ca] = neigh
                                    ind += 1
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    log_info("Average contact neighs: $(sum(neighbors .> 0.5)/size(x, 2))")
end

"""
    update_contact_neighs!(y, type, RM::ContactModel11, device::Type{Val{:cpu}}; max_part=30)

Update the neighbor list for contact force calculation (1-1 interaction).

# Arguments
- `y`: Positions of material points
- `type`: Type of material points
- `RM::ContactModel11`: Repulsion model for 1-1 interaction
- `device::Type{Val{:cpu}}`: Device type (CPU)
- `max_part=30`: Maximum number of particles in a cell

# Output
- No return value. The function updates the neighbor list in the `RM` object.

"""
function update_contact_neighs!(y, type, RM::ContactModel11, device::Type{Val{:cpu}}; max_part=30)
    mask = false
    for j in RM.type
        mask = mask .| (type .== j)
    end
    x = @view y[:, mask]
    neighbors = RM.neighs
    search_distance = RM.distance
    equi_dist = RM.equi_dist
    family = RM.material.family
    intact = RM.material.intact
    mask = any((family .> 0) .& map(~, intact) ; dims=1)
    fill!(neighbors, 0)
    update_contact_neighs!(mask, neighbors, x, search_distance, equi_dist, family, intact, max_part)
end

"""
    update_contact_neighs!(y, type, RM::ContactModel11, device::Type{Val{:cuda}}; max_part=30)

Update neighbor list for contact force calculation (1-1 interaction).

# Arguments
- `y`: Positions of material points
- `type`: Type of material points
- `RM::ContactModel11`: Repulsion model for 1-1 interaction
- `device::Type{Val{:cuda}}`: Device type (CUDA)
- `max_part=30`: Maximum number of particles in a cell

# Output
- No return value. The function updates the neighbor list in the `RM` object.

"""
function update_contact_neighs!(y, type, RM::ContactModel11, device::Type{Val{:cuda}}; max_part=30)
    mask = false
    for j in RM.type
        mask = mask .| (type .== j)
    end
    mask = Array(mask)
    x = Array(y[:, mask])
    neighbors = Array(RM.neighs)
    search_distance = RM.distance
    equi_dist = RM.material.particle_size
    family = Array(RM.material.family)
    intact = Array(RM.material.intact)
    fill!(neighbors,0)
    update_contact_neighs!(mask, neighbors, x, search_distance, equi_dist, family, intact, max_part)
    RM.neighs .= CuArray(neighbors)
    nothing
end

###############################################
# _cudaconvert
###############################################

"""
    _cudaconvert(x::Vector{T}) where T <: Union{ContactModel11,ContactModel12}

Converts a vector of `ContactModel11` or `ContactModel12` objects to CUDA-compatible types.
"""
function _cudaconvert(x::Vector{T}) where T <: Union{ContactModel11,ContactModel12}
    _cudaconvert.(x)
end


"""
    _cudaconvert(x::T) where T <: Union{ContactModel11,ContactModel12}

Converts a single `ContactModel11` or `ContactModel12` object to a CUDA-compatible type.
"""
function _cudaconvert(x::T) where T <: Union{ContactModel11} # ContactModel12 not implemented yet
    function fn(x, k)
        if k!=:material
            return _cudaconvert(getfield(x, k))
        else
            return getfield(x, k)
        end
    end
    T((fn(x, k) for k in fieldnames(T))...)
end

###############################################
# unit conversion
###############################################

function uconvert_to(DIMS, RM::T) where T <: ContactModel
    args = (uconvert_to(DIMS, getfield(RM, k)) for k in fieldnames(T))
    return T(args...)
end

function ustrip(RM::T) where T <: ContactModel
    args = (ustrip(getfield(RM, k)) for k in fieldnames(T))
    return T(args...)
end

###############################################
# contact models
###############################################

include("./nonlinearspring.jl")
include("./linearspring.jl")
include("./shortrange.jl")


#
