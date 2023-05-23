export force_density_T, BondBasedMaterial, BondBasedSpecific, PeridynamicsMaterial

"""
Specific bond based material type.
"""
struct BondBasedSpecific{T <: AbstractFloat} <: SpecificMaterial
    bond_stiffness::AbstractArray{T,2}
    critical_stretch::AbstractArray{T,2}
    density::AbstractArray{T,1}
    func::Function
end


function BondBasedSpecific(bond_stiffness::AbstractArray, critical_stretch::AbstractArray, density::AbstractArray; func=nothing)
    if isa(func, Nothing)
        func = bond_force
    end
    BondBasedSpecific(make_matrix(bond_stiffness), make_matrix(critical_stretch), density, func)
end

function BondBasedSpecific(bond_stiffness::Real, critical_stretch::Real, density::Real; kwargs...)
    BondBasedSpecific([bond_stiffness], [critical_stretch], [density]; kwargs...)
end

struct BondBasedMaterial <: PeridynamicsMaterial
    @PeridynamicsMaterial_gf
    specific::BondBasedSpecific
end

function PeridynamicsMaterial(name, type, bid, gen, spc::BondBasedSpecific)
    BondBasedMaterial(name, type, bid, gen, spc)
end


"""
    force_density_T(mat::BondBasedMaterial)

Calculates force density (actually acceleration) for bond based material type.
"""
function force_density_T(f, y::AbstractArray{Float64,2}, limits, mat::BondBasedMaterial; kwargs...)
    force_density_T(f, y, limits, mat, :cpu; kwargs)
end

@inline function bond_force(s, bond_stiffness)
    return bond_stiffness * s
end

function force_density_T(f, y::AbstractArray{Float64,2}, limits, mat::BondBasedMaterial, ::Type{Val{:cpu}}; particles=nothing)
    types = mat.general.type
    x = mat.general.x
    intact = mat.general.intact
    family = mat.general.family
    skip_bb = mat.general.skip_bb
    func = mat.specific.func
    stiffness = mat.specific.bond_stiffness

    if isnothing(particles)
        _N = 1:size(family, 2) + limits[1] - 1
    else
        _N = particles + limits[1] - 1
    end

    M = size(family, 1)
    ARGS = map((i) -> (i, 1:M), _N)

    function with_if_cal_force_ij(i, k)
        if intact[k,i]
            j = family[k,i]
            X = [x[1,j]-x[1,i],x[2,j]-x[2,i],x[3,j]-x[3,i]]
            Y = [y[1,j]-y[1,i],y[2,j]-y[2,i],y[3,j]-y[3,i]]
            _X = get_magnitude(X)
            _Y = get_magnitude(Y)
            ext = _Y - _X
            s = ext/_X
            type1 = types[i]- mat.type.start + 1
            type2 = types[j] - mat.type.start + 1
            if  skip_bb || ((type1!=-1) && (type2!=-1) && (s < mat.specific.critical_stretch[type1, type2]))
                return func(s, stiffness[type1, type2]) * mat.general.volume[j] * (Y/_Y)
            else
                intact[k,i] = false
                return [0.0, 0.0, 0.0]
            end
        else
            return [0.0, 0.0, 0.0]
        end
    end

    inner_map(i, inds) = map_reduce((j)-> with_if_cal_force_ij(i,j), +, inds; init=[0.0, 0.0, 0.0])

    # for x in ARGS
    #     i, inds = x
    #     for j in inds
    #         with_if_cal_force_ij(i,j)
    #     end
    # end

    # # inner_map(i, inds) = mapreduce((j)-> with_if_cal_force_ij(i,j), +, inds)

    outer_map(ARGS) = map((x)->inner_map(x[1], x[2]), ARGS)

    return hcat(outer_map(ARGS)...)
end



"""
    force_density_T(mat::BondBasedMaterial)

Calculates force density (actually acceleration) for bond based material type.
"""
function force_density_T(force::AbstractArray{Float64,2}, y::AbstractArray{Float64,2}, limits, mat::BondBasedMaterial, ::Type{Val{:cuda}}; particles=nothing)

    force = CUDA.CuArray(zeros(eltype(y), size(y)))
    y = CUDA.CuArray(y)
    types = CUDA.CuArray(mat.general.type)
    tstart = mat.type.start
    x = CUDA.CuArray(mat.general.x)
    intact = CUDA.CuArray(mat.general.intact)
    family = CUDA.CuArray(mat.general.family)
    volume = CUDA.CuArray(mat.general.volume)
    stiffness = CUDA.CuArray(mat.specific.bond_stiffness)
    critical_stretch = CUDA.CuArray(mat.specific.critical_stretch)
    skip_bb = CUDA.CuArray([mat.general.skip_bb])
    func = mat.specific.func

    if isnothing(particles)
        _N = 1:size(family, 2)
    else
        _N = particles
    end

    _cuN = CUDA.CuArray(_N)

    function cal_force(y, x, force, family, intact, stiffness, volume, types, _N, critical_stretch, skip_bb, func)
        index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        stride = blockDim().x * gridDim().x

        for ind in index:stride:length(_N)
            i = _N[ind]
            for k in 1:size(family, 1)
                if intact[k,i]
                    j = family[k,i]
                    _X = sqrt((x[1,j]-x[1,i])^2 + (x[2,j]-x[2,i])^2 + (x[3,j]-x[3,i])^2)
                    Y1 = y[1,j]-y[1,i]
                    Y2 = y[2,j]-y[2,i]
                    Y3 = y[3,j]-y[3,i]
                    _Y = sqrt((Y1)^2 + (Y2)^2 + (Y3)^2)
                    ext = _Y - _X
                    s = ext/_X
                    type1 = types[i] - tstart+ 1
                    type2 = types[j] - tstart + 1
                    if skip_bb[1] || ((type1!=-1) && (type2!=-1) && (s < critical_stretch[type1, type2]))
                        t = func(s, stiffness[type1, type2])
                        force[1, i] += t*volume[j] * Y1/_Y
                        force[2, i] += t*volume[j] * Y2/_Y
                        force[3, i] += t*volume[j] * Y3/_Y
                    else
                        intact[k,i] = false
                    end
                end
            end
        end
        return nothing
    end

    kernel = CUDA.@cuda launch=false cal_force(y, x, force, family, intact, stiffness, volume, types, _cuN, critical_stretch, skip_bb, func)
    config = launch_configuration(kernel.fun)
    nthreads = Base.min(length(_N), config.threads)
    nblocks =  cld(length(_N), nthreads)

    CUDA.@sync kernel(y, x, force, family, intact, stiffness, volume, types, _cuN, critical_stretch, skip_bb, func; threads=nthreads, blocks=nblocks)
    mat.general.intact .= Array(intact)
    return  CUDA.Array(force)[:, _N]
end



