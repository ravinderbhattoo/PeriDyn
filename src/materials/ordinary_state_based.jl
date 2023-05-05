export OrdinaryStateBasedMaterial, OrdinaryStateBasedSpecific, force_density_T, PeridynamicsMaterial

mutable struct OrdinaryStateBasedSpecific <: SpecificMaterial
    bulk_modulus::AbstractArray{Float64, 2}
    shear_modulus::AbstractArray{Float64, 2}
    critical_stretch::AbstractArray{Float64, 2}
    density::AbstractArray{Float64, 1}
    kernel::Any
    theta::AbstractArray{Float64, 1}
end

function init(spc::OrdinaryStateBasedSpecific, gen::GeneralMaterial)
    volume = gen.volume
    spc.theta = zeros(eltype(volume), size(volume))
    spc
end


"""
    OrdinaryStateBasedSpecific(bulk_modulus::AbstractArray{Float64, 1}, shear_modulus::AbstractArray{Float64,1}, critical_stretch::AbstractArray{Float64,1}, density::AbstractArray{Float64,1})
Specific ordinary state based materail type.
"""
function OrdinaryStateBasedSpecific(bulk_modulus::AbstractArray{Float64, 1}, shear_modulus::AbstractArray{Float64,1}, critical_stretch::AbstractArray{Float64,1}, density::AbstractArray{Float64,1})
    return OrdinaryStateBasedSpecific(make_matrix(bulk_modulus), make_matrix(shear_modulus), make_matrix(critical_stretch), density, nothing, [])
end

"""
    OrdinaryStateBasedMaterial(type::UnitRange{Int64}, general::GeneralMaterial, specific::OrdinaryStateBasedSpecific)
"""
struct OrdinaryStateBasedMaterial <: PeridynamicsMaterial
    @PeridynamicsMaterial_gf
    specific::OrdinaryStateBasedSpecific
end

function PeridynamicsMaterial(name, type, bid, gen, spc::OrdinaryStateBasedSpecific)
    OrdinaryStateBasedMaterial(name, type, bid, gen, spc)
end



"""
    force_density_T(mat::BondBasedMaterial)

Calculates force density (actually acceleration) for bond based material type.
"""
function force_density_T(f, y::AbstractArray{Float64,2}, limits, mat::OrdinaryStateBasedSpecific; kwargs...)
    force_density_T(f, y, limits, mat, :cpu; kwargs)
end


"""
    force_density_T(f, y, limits, mat::OrdinaryStateBasedMaterial, device::Type{Val{:cpu}}; particles=nothing)

Calculates force density (actually acceleration) for ordinary state based material type.
"""
function force_density_T(f_, y_, limits, mat::OrdinaryStateBasedMaterial, device::Type{Val{:cpu}}; particles=nothing)
    types = mat.general.type
    gen = mat.general
    m = gen.weighted_volume
    intact = gen.intact
    family = gen.family
    volume = gen.volume
    theta = mat.specific.theta
    x = gen.x
    y = y_[:, limits[1]:limits[2]]

    if isnothing(particles)
        N = 1:size(family, 2)
        N_ = N .+ (limits[1] .- 1)
    else
        N_ = particles
        N = N_ .- (limits[1] .- 1)
    end

    M = size(family, 1)
    ARGS = map((i) -> (i, 1:M), N)

    dilatation!(theta, y, x, intact, family, volume, m,
                            mat.general.particle_size,
                            mat.general.horizon,
                            device)

    @inbounds function with_if_cal_force_ij(i, k)
        if intact[k,i]
            j = family[k,i]
            X = get_ij(j,i,x)
            Y = get_ij(j,i,y)
            yij = get_magnitude(Y)
            xij = get_magnitude(X)

            extention = yij - xij
            s = extention/xij

            wij = influence_function(X)
            wji = influence_function(-X)
            type1 = types[i] - mat.type.start + 1
            type2 = types[j] - mat.type.start + 1
            if s < mat.specific.critical_stretch[type1, type2]
                K = mat.specific.bulk_modulus[type1, type2]
                G = mat.specific.shear_modulus[type1, type2]
                a_, b_ = wij/m[i], wji/m[j]
                t =  (3*K-5*G) * (theta[i]*xij*a_ + theta[j]*xij*b_)  + 15*G*extention*(a_ + b_)
                t = t*volume[j]
                return t / yij *Y
            else
                intact[k,i] = false
                return [0.0, 0.0, 0.0]
            end
        else
            return [0.0, 0.0, 0.0]
        end
    end


    inner_map(i, inds) = map_reduce((j)-> with_if_cal_force_ij(i,j), +, inds)

    # inner_map(i, inds) = mapreduce((j)-> with_if_cal_force_ij(i,j), +, inds)
    outer_map(ARGS) = map((x)->inner_map(x[1], x[2]), ARGS)

    f_[:, N_] .+= hcat(outer_map(ARGS)...)
end


"""
    force_density_T(force::AbstractArray{Float64,2}, y::AbstractArray{Float64,2}, limits, mat::OrdinaryStateBasedMaterial, device::Type{Val{:cuda}}; particles=nothing)

Calculates force density (actually acceleration) for bond based material type.
"""
function force_density_T(force::AbstractArray{Float64,2}, y::AbstractArray{Float64,2}, limits, mat::OrdinaryStateBasedMaterial, device::Type{Val{:cuda}}; particles=nothing)

    types = mat.general.type
    tstart = mat.type[1]
    x = mat.general.x
    intact = mat.general.intact
    family = mat.general.family
    volume = mat.general.volume
    m = mat.general.weighted_volume
    critical_stretch = mat.specific.critical_stretch
    bulk_modulus = mat.specific.bulk_modulus
    shear_modulus = mat.specific.shear_modulus

    if isnothing(particles)
        __N = collect(1:size(family, 2))
        _N = __N .+ (limits[1] .- 1)
    else
        _N = particles
        __N = _N .- (limits[1] .- 1)
    end

    theta = mat.specific.theta
    l1 = limits[1] - 1
    dilatation!(theta, l1, y, x, intact, family, volume, m,
                mat.general.particle_size, mat.general.horizon, device)


    __N = CUDA.CuArray(__N)
    _N = CUDA.CuArray(_N)


    if isa(mat.specific.kernel, Nothing)

        function cal_force(force, y, x, family, intact, volume, types, _N, __N,
                            critical_stretch, theta, m, bulk_modulus, shear_modulus)
            index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
            stride = blockDim().x * gridDim().x

            for ind in index:stride:length(_N)
                _i = _N[ind]
                __i = __N[ind]
                for k in 1:size(family, 1)
                    if intact[k, __i]
                        j = family[k, __i]
                        _j = j + _i - __i

                        xij = sqrt((x[1,j]-x[1,__i])^2 + (x[2,j]-x[2,__i])^2 + (x[3,j]-x[3,__i])^2)

                        Y1 = y[1,_j] - y[1,_i]
                        Y2 = y[2,_j] - y[2,_i]
                        Y3 = y[3,_j] - y[3,_i]

                        _Y = sqrt((Y1)^2 + (Y2)^2 + (Y3)^2)
                        extention = _Y - xij
                        s = extention/xij

                        wij = 1/xij
                        wji = 1/xij

                        type1 = types[__i] - tstart + 1
                        type2 = types[j] - tstart + 1
                        if (type1!=-1) && (type2!=-1) && (s < critical_stretch[type1, type2])
                            K = bulk_modulus[type1, type2]
                            G = shear_modulus[type1, type2]
                            a_, b_ = wij/m[__i], wji/m[j]
                            t =  (3*K-5*G) * (theta[__i]*xij*a_ + theta[j]*xij*b_)  + 15*G*extention*(a_ + b_)
                            t = t * volume[j]
                            force[1, _i] += t * Y1/_Y
                            force[2, _i] += t * Y2/_Y
                            force[3, _i] += t * Y3/_Y
                        else
                            intact[k,__i] = false
                        end
                    end
                end
            end
            return nothing
        end

        kernel = CUDA.@cuda launch=false cal_force(force, y, x, family, intact, volume, types, _N, __N,
                                            critical_stretch, theta, m, bulk_modulus, shear_modulus)
        config = launch_configuration(kernel.fun)
        nthreads = Base.min(length(_N), config.threads)
        nblocks =  cld(length(_N), nthreads)

        function fn(args...)
            CUDA.@sync kernel(args...; threads=nthreads, blocks=nblocks)
        end
        mat.specific.kernel = fn
    end


    mat.specific.kernel(force, y, x, family, intact, volume, types, _N, __N,
                        critical_stretch, theta, m, bulk_modulus, shear_modulus)
end



