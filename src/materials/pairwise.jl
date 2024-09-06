export PairWiseMaterial, force_density_T


PairWiseMaterial = Union{BondBasedMaterial,PairwiseNNMaterial}


function force_density_T!(y::Array{Float64,2}, mat::T, ::Type{Val{:cpu}}; particles=nothing) where T <: PairWiseMaterial
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
            _X = get_magnitude(X)
            _Y = get_magnitude(Y)
            ext = _Y - _X
            s = ext/_X
            type1 = types[i]- mat.type.start + 1
            type2 = types[j] - mat.type.start + 1
            if  skip_bb || ((type1!=-1) && (type2!=-1) && (s < mat.specific.critical_stretch[type1, type2]))
                return bond_force(s, mat.specific) * mat.general.volume[j] * (Y/_Y)
            else
                intact[k,i] = false
                return [0.0, 0.0, 0.0]
            end
        else
            return [0.0, 0.0, 0.0]
        end
    end

    inner_map(i, inds) = map_reduce((j)-> with_if_cal_force_ij(i,j), +, inds)

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
