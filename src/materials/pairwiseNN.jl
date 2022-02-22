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
    name::String
    type::UnitRange{Int64}
    general::GeneralMaterial
    specific::PairwiseNNSpecific
end

function PeridynamicsMaterial(name, type, gen, spc::PairwiseNNSpecific)
    PairwiseNNMaterial(name, type, gen, spc) 
end

function cal_force_ij(f, Y, inp)
    M = Y./get_magnitude(Y)
    t = f(inp)
    return t.*M    
end

"""
    force_density_T(mat::PairwiseNNMaterial)

Calculates force density (actually acceleration) for bond based material type.
"""
function force_density_T(y::Array{Float64,2}, mat::PairwiseNNMaterial; particles=nothing)
    types = mat.general.type
    x = mat.general.x 
    intact = mat.general.intact
    family = mat.general.family
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
            if s < mat.specific.critical_stretch[type1, type2]
                return cal_force_ij((x) -> first(mat.specific.NNs[type1, type2](x)), Y, [ext, _X])
            else
                return [0.0, 0.0, 0.0]
            end            
        else
            return [0.0, 0.0, 0.0]
        end
    end

    # inner_map(i, inds) = Folds.mapreduce((j)-> with_if_cal_force_ij(i,j), +, inds)

    inner_map(i, inds) = sum(map((j)-> with_if_cal_force_ij(i,j), inds))
    outer_map(ARGS) = map((x)->inner_map(x[1], x[2]), ARGS)

    return hcat(outer_map(ARGS)...)
end

#
