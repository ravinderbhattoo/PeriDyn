export OrdinaryStateBasedMaterial, OrdinaryStateBasedSpecific, force_density_T, PeridynamicsMaterial

"""
Specific ordinary state based materail type.
"""
struct OrdinaryStateBasedSpecific <: SpecificMaterial
    bulk_modulus::Array{Float64, 2}
    shear_modulus::Array{Float64, 2}
    critical_stretch::Array{Float64, 2}
    density::Array{Float64, 1}
end

"""
    OrdinaryStateBasedSpecific(bulk_modulus::Array{Float64, 1}, shear_modulus::Array{Float64,1}, critical_stretch::Array{Float64,1}, density::Array{Float64,1})
Specific ordinary state based materail type.
"""
function OrdinaryStateBasedSpecific(bulk_modulus::Array{Float64, 1}, shear_modulus::Array{Float64,1}, critical_stretch::Array{Float64,1}, density::Array{Float64,1})
    return OrdinaryStateBasedSpecific(make_matrix(bulk_modulus), make_matrix(shear_modulus), make_matrix(critical_stretch), density)
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
function force_density_T(y::Array{Float64,2}, mat::OrdinaryStateBasedSpecific; kwargs...)
    force_density_T(y, mat, :cpu; kwargs)
end



# """
#     force_density_T(y::Array{Float64,2},mat::OrdinaryStateBasedMaterial)

# Calculates force density (actually acceleration) for ordinary state based material type.
# """
# function PeriDyn.force_density_T(y::Array{Float64,2}, mat::OrdinaryStateBasedMaterial; particles=nothing)

#     force = Zygote.Buffer(y)

#     types = mat.general.type
#     gen = mat.general
#     m = gen.weighted_volume
#     intact = gen.intact

#     intact = Zygote.Buffer(intact)

#     family = gen.family
#     x = gen.x
#     volume = gen.volume

#     if isnothing(particles)
#         _N = 1:size(family, 2)
#     else
#         _N = particles
#     end

#     theta = dilatation(y, gen, m)

#     N = size(x, 2)
#     for i in _N
#         force[1,i] = 0.0
#         force[2,i] = 0.0
#         force[3,i] = 0.0
#         for k in 1:size(family,1)
#             j = family[k,i]
#             if !intact[k,i]
#                 intact[k,i] = false
#             else
#                 X = get_ij(j,i,x)
#                 Y = get_ij(j,i,y)
#                 yij = get_magnitude(Y)
#                 xij = get_magnitude(X)

#                 extention = yij - xij
#                 wij = influence_function(X)
#                 wji = influence_function(-X)
#                 type1 = types[i] - mat.type.start + 1
#                 type2 = types[j] - mat.type.start + 1
#                 if (extention/xij) < mat.specific.critical_stretch[type1, type2]
#                     K = mat.specific.bulk_modulus[type1, type2]
#                     G = mat.specific.shear_modulus[type1, type2]
#                     a_, b_ = wij/m[i], wji/m[j]
#                     t =  (3*K-5*G) * (theta[i]*xij*a_ + theta[j]*xij*b_)  + 15*G*extention*(a_ + b_)
#                     t = t*volume[j] / yij
#                     force[1,i] += t*Y[1]
#                     force[2,i] += t*Y[2]
#                     force[3,i] += t*Y[3]
#                 else
#                     intact[k,i] = false
#                 end
#             end
#         end
#     end
#     return copy(force)
# end



"""
    force_density_T(y::Array{Float64,2},mat::OrdinaryStateBasedMaterial)

Calculates force density (actually acceleration) for ordinary state based material type.
"""
function force_density_T(y::Array{Float64,2}, mat::OrdinaryStateBasedMaterial, device::Type{Val{:cpu}}; particles=nothing)

    types = mat.general.type
    gen = mat.general
    m = gen.weighted_volume
    intact = gen.intact
    family = gen.family
    x = gen.x
    volume = gen.volume

    if isnothing(particles)
        _N = 1:size(family, 2)
    else
        _N = particles
    end

    theta = dilatation(y, x, intact, family, volume, m,
                            mat.general.particle_size,
                            mat.general.horizon,
                            device)

    M = size(family, 1)
    ARGS = map((i) -> (i, 1:M), _N)

    function with_if_cal_force_ij(i, k)
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

    return hcat(outer_map(ARGS)...)


end


"""
    force_density_T(mat::OrdinaryStateBasedMaterial)

Calculates force density (actually acceleration) for bond based material type.
"""
function force_density_T(y::Array{Float64,2}, mat::OrdinaryStateBasedMaterial, device::Type{Val{:cuda}}; particles=nothing)

    force = CUDA.CuArray(zeros(eltype(y), size(y)))
    y = CUDA.CuArray(y)
    types = CUDA.CuArray(mat.general.type)
    tstart = mat.type.start
    x = CUDA.CuArray(mat.general.x)
    intact = CUDA.CuArray(mat.general.intact)
    family = CUDA.CuArray(mat.general.family)
    volume = CUDA.CuArray(mat.general.volume)
    m = CUDA.CuArray(mat.general.weighted_volume)
    critical_stretch = CUDA.CuArray(mat.specific.critical_stretch)
    bulk_modulus = CUDA.CuArray(mat.specific.bulk_modulus)
    shear_modulus = CUDA.CuArray(mat.specific.shear_modulus)

    if isnothing(particles)
        _N = 1:size(family, 2)
    else
        _N = particles
    end

    theta = CUDA.CuArray(zeros(eltype(volume), size(volume)))

    dilatation!(theta, y, x, intact, family, volume, m,
                            mat.general.particle_size,
                            mat.general.horizon,
                            device)

    _N = CUDA.CuArray(_N)

    function cal_force(y, x, force, family, intact, volume, types, _N,
                        critical_stretch, theta, m, bulk_modulus, shear_modulus)
        index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        stride = blockDim().x * gridDim().x

        for ind in index:stride:length(_N)
            i = _N[ind]
            for k in 1:size(family, 1)
                if intact[k,i]
                    j = family[k,i]
                    xij = sqrt((x[1,j]-x[1,i])^2 + (x[2,j]-x[2,i])^2 + (x[3,j]-x[3,i])^2)
                    Y1 = y[1,j]-y[1,i]
                    Y2 = y[2,j]-y[2,i]
                    Y3 = y[3,j]-y[3,i]
                    _Y = sqrt((Y1)^2 + (Y2)^2 + (Y3)^2)
                    extention = _Y - xij
                    s = extention/xij

                    wij = 1/xij
                    wji = 1/xij

                    type1 = types[i] - tstart+ 1
                    type2 = types[j] - tstart + 1
                    if (type1!=-1) && (type2!=-1) && (s < critical_stretch[type1, type2])
                        K = bulk_modulus[type1, type2]
                        G = shear_modulus[type1, type2]
                        a_, b_ = wij/m[i], wji/m[j]
                        t =  (3*K-5*G) * (theta[i]*xij*a_ + theta[j]*xij*b_)  + 15*G*extention*(a_ + b_)
                        t = t * volume[j]
                        force[1, i] += t * Y1/_Y
                        force[2, i] += t * Y2/_Y
                        force[3, i] += t * Y3/_Y
                    else
                        intact[k,i] = false
                    end
                end
            end
        end
        return nothing
    end

    kernel = CUDA.@cuda launch=false cal_force(y, x, force, family, intact, volume, types, _N,
                                        critical_stretch, theta, m, bulk_modulus, shear_modulus)
    config = launch_configuration(kernel.fun)
    nthreads = Base.min(length(_N), config.threads)
    nblocks =  cld(length(_N), nthreads)

    CUDA.@sync kernel(y, x, force, family, intact, volume, types, _N,
                        critical_stretch, theta, m, bulk_modulus, shear_modulus; threads=nthreads, blocks=nblocks)
    mat.general.intact .= Array(intact)
    return  CUDA.Array(force)
end







#
