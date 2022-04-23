abstract type RepulsionModel11 end
abstract type RepulsionModel12 end

@def RepulsionModel_gf begin
    equi_dist::Float64
    distance::Float64
    horizons::Tuple
    hor::Float64
    neighs::Array{Int64,2}
    max_neighs::Int64
end

@def RepulsionModel12_gf begin
    @RepulsionModel_gf
    pair::Vector{UnitRange{Int64}}
end

@def RepulsionModel11_gf begin
    @RepulsionModel_gf
    type::UnitRange{Int64}
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
function short_range_repulsion!(y, f, type, vol, RM::RepulsionModel12)
    mask1 = false
    for i in RM.pair[1]
        mask1 = mask1 .| (type .== i)
    end
    mask2 = false
    for i in RM.pair[2]
        mask2 = mask2 .| (type .== i)
    end
    f1 = f[:, mask1]
    f2 = f[:, mask2]
    x1 = y[:, mask1]
    x2 = y[:, mask2]
    vol1 = vol[mask1]
    vol2 = vol[mask2]
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
    f[:, mask2] .= f2 
    return nothing
end


"""
    short_range_repulsion!(y,f,type,RepusionModel)

Updates (inplace) the repulsive acceleration of material points.
"""
function short_range_repulsion!(y, f, type, vol, RM::RepulsionModel11)
    mask1 = false
    for j in RM.type
        mask1 = mask1 .| (type .== j)
    end
    f1 = f[:,mask1]
    x1 = y[:,mask1]
    vol1 = vol[mask1]
    for i in 1:size(RM.neighs,2)
        for k in 1:size(RM.neighs,1)
            j = RM.neighs[k,i]
            if j>0
                force = repulsion_force(x1[:, i] .- x1[:, j], RM)
                f1[:,i] .+= force / vol1[j]
                f1[:,j] .-= force / vol1[i]
            end
            if j==0 break end
        end
    end
    f[:,mask1] .= f1 
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
function collision_box(x1::Array{Float64,2}, x2::Array{Float64,2}, skin::Float64)
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
function update_repulsive_neighs!(y,type,RM::RepulsionModel12)
    _mask1 = false
    for i in RM.pair[1]
        _mask1 = _mask1 .| (type .== i)
    end
    _mask2 = false
    for i in RM.pair[2]
        _mask2 = _mask2 .| (type .== i)
    end
    x11 = y[:, _mask1]
    x22 = y[:, _mask2]
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
        println("Average repulsive neighs: $(sum(RM.neighs .> 0.5)/size(x1, 2))")
    else
        RM.neighs[1:end, 1:end] .= 0
        println("No repulsive neighs.")
    end
end


"""
    update_repulsive_neighs!(y,type,RepulsionModel11)

Update neighbour list for repulsive force calculation (1-1 interaction).

"""
function update_repulsive_neighs!(y, type, RM::RepulsionModel11)
    mask = false
    for j in RM.type
        mask = mask .| (type .== j)
    end
    x = y[:, mask]
    fill!(RM.neighs,0)
    cells, cell_neighs = get_cells(x, RM.distance)
    Threads.@threads for cell_i in 1:length(cells)
        for ca_id in 1:length(cells[cell_i])
            ind = 1
            ca = cells[cell_i][ca_id]
            a1,b1,c1 = x[1,ca],x[2,ca],x[3,ca]
            for neigh_id in 1:length(cell_neighs[cell_i])
                neighs = cell_neighs[cell_i][neigh_id]
                for fa_id in 1:length(cells[neighs])
                    fa = cells[neighs][fa_id]
                    noskip = true
                    for k in 1:size(RM.material.family,1)
                        if fa < RM.material.family[k,ca]
                            break
                        end
                        if fa==RM.material.family[k,ca]
                            if RM.material.intact[k,ca]
                                noskip = false
                            end
                            break
                        end
                    end
                    if noskip
                        a2,b2,c2 = x[1,fa], x[2,fa], x[3,fa]
                        if 1e-4<((a1-a2)*(a1-a2)+(b1-b2)*(b1-b2)+(c1-c2)*(c1-c2))<RM.distance^2
                            RM.neighs[ind,ca] = fa
                            ind += 1
                        end
                    end
                end
            end
        end
    end
    println("Average repulsive neighs: $(sum(RM.neighs .> 0.5)/size(x, 2))")
end


include("./LJ.jl")
include("./shortrange.jl")
include("./nonlinear.jl")
include("./linear.jl")


#
