export force_density_T, PairwiseNNMaterial, PairwiseNNSpecific

"""
Specific bond based material type (NN) .
"""
struct PairwiseNNSpecific <: SpecificMaterial
    NNs::SymMat
    critical_stretch::Array{Float64, 2}
    density::Array{Float64, 1}
end

function PairwiseNNSpecific(layers, critical_stretch, density::Array{Float64, 1}; act=Flux.relu)
    PairwiseNNSpecific(make_NN(layers, length(critical_stretch); act=act), make_matrix(critical_stretch), density)
end

"""
Bond based material type (NN).
"""
struct PairwiseNNMaterial <: PeridynamicsMaterial
    @PeridynamicsMaterial_gf
    specific::PairwiseNNSpecific
end

function PeridynamicsMaterial(name, type, bid, gen, spc::PairwiseNNSpecific)
    PairwiseNNMaterial(name, type, bid, gen, spc)
end


"""
    force_density_T!(mat::PairwiseNNSpecific)

Calculates force density (actually acceleration) for bond based material type.
"""
function force_density_T!(y::Array{Float64,2}, mat::PairwiseNNSpecific; kwargs...)
    force_density_T!(y, mat, :cpu; kwargs)
end

# """
#     force_density_T!(mat::PairwiseNNMaterial)

# Calculates force density (actually acceleration) for bond based material type.
# """
# function force_density_T!(y::Array{Float64,2}, mat::PairwiseNNMaterial; particles=nothing)
#     types = mat.general.type
#     x = mat.general.x
#     intact = mat.general.intact
#     family = mat.general.family
#     if isnothing(particles)
#         _N = 1:size(family, 2)
#     else
#         _N = particles
#     end

#     M = size(family, 1)
#     ARGS = map((i) -> (i, 1:M), _N)

#     function with_if_cal_force_ij(i, k)
#         if intact[k,i]
#             j = family[k,i]
#             X = [x[1,j]-x[1,i],x[2,j]-x[2,i],x[3,j]-x[3,i]]
#             Y = [y[1,j]-y[1,i],y[2,j]-y[2,i],y[3,j]-y[3,i]]
#             _X = get_magnitude(X)
#             _Y = get_magnitude(Y)
#             ext = _Y - _X
#             # s = ext/_X
#             type1 = types[i]- mat.type.start + 1
#             type2 = types[j] - mat.type.start + 1
#             # if s < mat.specific.critical_stretch[type1, type2]
#                 # return cal_force_ij((x) -> first(mat.specific.NNs[type1, type2](x)), Y, [ext, _X])
#                 return first(mat.specific.NNs[type1, type2]([ext])) * Y / _Y
#             # else
#             #     return [0.0, 0.0, 0.0]
#             # end
#         else
#             return [0.0, 0.0, 0.0]
#         end
#     end

#     inner_map(i, inds) = map_reduce((j)-> with_if_cal_force_ij(i,j), +, inds)

#     # inner_map(i, inds) = mapreduce((j)-> with_if_cal_force_ij(i,j), +, inds)

#     outer_map(ARGS) = map((x)->inner_map(x[1], x[2]), ARGS)

#     return hcat(outer_map(ARGS)...)
# end



function force_density_T!(y::Array{Float64,2}, mat::PairwiseNNMaterial, ::Type{Val{:cpu}}; particles=nothing)
    types = mat.general.type
    x = mat.general.x
    intact = mat.general.intact
    family = mat.general.family
    skip_bb = mat.general.skip_bb

    if isnothing(particles)
        _N = 1:size(family, 2)
    else
        _N = particles
    end

    M = size(family, 1)
    ARGS = map((i) -> (i, 1:M), _N)

    function with_if_cal_force_ij(i, k)
        if intact[k,i]
            j = family[k,i]
            X = [x[1,j]-x[1,i],x[2,j]-x[2,i],x[3,j]-x[3,i]]
            Y = [y[1,j]-y[1,i],y[2,j]-y[2,i],y[3,j]-y[3,i]]
            _X =  get_magnitude(X)
            _Y =  get_magnitude(Y)
            ext = _Y - _X
            s = ext/_X
            type1 = types[i]- mat.type.start + 1
            type2 = types[j] - mat.type.start + 1
            if skip_bb || ((type1!=-1) && (type2!=-1) && (s < mat.specific.critical_stretch[type1, type2]))
                # return first(mat.specific.NNs[type1, type2]([s])) * Y / _Y
                return s, mat.general.volume[j] * Y / _Y
            else
                return s, [0.0, 0.0, 0.0]
            end
        else
            return 0.0, [0.0, 0.0, 0.0]
        end
    end

    function getdata()
        ind = 1
        data = []
        for i in _N
            push!(data, map((k)->with_if_cal_force_ij(i, k), 1:M))
            ind += 1
        end
        data
    end

    data =  Zygote.@ignore getdata()


    force = map((data_) -> map_reduce((k) -> first(mat.specific.NNs[1, 1]([data_[k][1]])) * data_[k][2], +, 1:M), data)
    return hcat(force...)

    # inner_map(i, inds) = PeriDyn.map_reduce((j)-> with_if_cal_force_ij(i,j), +, inds)

    # # for x in ARGS
    # #     i, inds = x
    # #     for j in inds
    # #         with_if_cal_force_ij(i,j)
    # #     end
    # # end

    # # # inner_map(i, inds) = mapreduce((j)-> with_if_cal_force_ij(i,j), +, inds)

    # outer_map(ARGS) = map((x)->inner_map(x[1], x[2]), ARGS)

    # return hcat(outer_map(ARGS)...)
end


# """
#     force_density_T!(mat::PairwiseNNMaterial)

# Calculates force density (actually acceleration) for bond based material type.
# """
# function force_density_T!(y::Array{Float64,2}, mat::PairwiseNNMaterial, ::Type{Val{:cpu}}; particles=nothing)
#     force_buf = zeros(eltype(y), size(y)...)
#     # force_buf = Zygote.Buffer(force)
#     types = mat.general.type
#     x = mat.general.x
#     intact = mat.general.intact
#     family = mat.general.family
#     if isnothing(particles)
#         _N = 1:size(family, 2)
#     else
#         _N = particles
#     end

#     M = size(family, 1)

#     function with_if_cal_force_ij(i, k)
#         if intact[k,i]
#             j = family[k,i]
#             X = [x[1,j]-x[1,i],x[2,j]-x[2,i],x[3,j]-x[3,i]]
#             Y = [y[1,j]-y[1,i],y[2,j]-y[2,i],y[3,j]-y[3,i]]
#             _X = get_magnitude(X)
#             _Y = get_magnitude(Y)
#             ext = _Y - _X
#             s = ext/_X
#             type1 = types[i]- mat.type.start + 1
#             type2 = types[j] - mat.type.start + 1
#             if s < mat.specific.critical_stretch[type1, type2]
#                 return first(mat.specific.NNs[type1, type2]([s])) * Y / _Y
#             else
#                 return [0.0, 0.0, 0.0]
#             end
#         else
#             return [0.0, 0.0, 0.0]
#         end
#     end

#     for i in _N
#         for k in 1:M
#             f_ =  with_if_cal_force_ij(i, k)
#             force_buf[1, i] += f_[1]
#             force_buf[2, i] += f_[2]
#             force_buf[3, i] += f_[3]
#         end
#     end
#     return force_buf
# end

function unroll_model(N)
    ex = Expr(
        :block,
        :(B = s)
        )
    for i in 3:3:N
        a = i - 2
        b = i - 1
        c = i
        tmp = quote
            W = args[$a]
            b = args[$b]
            act = args[$c]
            for i in 1:size(W, 1)
                for j in 1:size(W, 2)
                    b[i] = b[i] + W[i, j]*B[j]
                end
                b[i] = act(b[i])
            end
            B = b
        end
        (x -> push!(ex.args, x)).(tmp.args)
    end
    push!(ex.args, :(return B))

    exf = :((s, args...) -> nothing)
    (x -> push!(exf.args[2].args, x)).(ex.args)
    return ex
end

"""
    force_density_T!(mat::PairwiseNNMaterial)

Calculates force density (actually acceleration) for bond based material type.
"""
function force_density_T!(y::Array{Float64,2}, mat::PairwiseNNMaterial, ::Type{Val{:cuda}}; particles=nothing)

    force = CUDA.CuArray(zeros(eltype(y), size(y)))
    y = CUDA.CuArray(y)
    types = CUDA.CuArray(mat.general.type)
    tstart = mat.type.start
    x = CUDA.CuArray(mat.general.x)
    intact = CUDA.CuArray(mat.general.intact)
    family = CUDA.CuArray(mat.general.family)
    volume = CUDA.CuArray(mat.general.volume)
    critical_stretch = CUDA.CuArray(mat.specific.critical_stretch)
    skip_bb = CUDA.CuArray([mat.general.skip_bb])
    model = mat.specific.NNs[1,1]

    if isnothing(particles)
        _N = 1:size(family, 2)
    else
        _N = particles
    end

    _cuN = CUDA.CuArray(_N)

    wbas = []
    for i in model
        if isa(i, Dense)
            push!(wbas, CUDA.CuArray(i.weight))
            push!(wbas, CUDA.CuArray(i.bias))
            push!(wbas, i.σ)
        else
            push!(wbas, CUDA.CuArray([1.0;;]))
            push!(wbas, CUDA.CuArray([1.0]))
            push!(wbas, i)
        end
    end

    function func(s, args...)
        f = eval(unroll_model(length(size_wbas)))
        return f(s, args...)
    end

    # function func(s, args...)
    #     s * 33009.91412276346
    # end

    function cal_force(y, x, force, family, intact, volume, types, _N, critical_stretch, skip_bb, func, args...)
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
                        t = func(s, args...)
                        t = t * volume[j] / _Y

                        force[1, i] += t * Y1
                        force[2, i] += t * Y2
                        force[3, i] += t * Y3
                    else
                        intact[k,i] = false
                    end
                end
            end
        end
        return nothing
    end

    kernel = CUDA.@cuda launch=false cal_force(y, x, force, family, intact, volume, types, _cuN, critical_stretch, skip_bb, func, wbas...)
    config = launch_configuration(kernel.fun)
    nthreads = Base.min(length(_N), config.threads)
    nblocks =  cld(length(_N), nthreads)

    CUDA.@sync kernel(y, x, force, family, intact, volume, types, _cuN, critical_stretch, skip_bb, func, wbas...; threads=nthreads, blocks=nblocks)
    mat.general.intact .= Array(intact)
    return  CUDA.Array(force)[:, _N]
end




