export force_density_T, BondBasedMaterial, BondBasedSpecific

"""
Specific bond based material type.
"""
struct BondBasedSpecific <: SpecificMaterial
    bond_stiffness::Array{Float64,2}
    critical_stretch::Array{Float64, 2}
    density::Array{Float64, 1}
end

function BondBasedSpecific(S, critical_stretch, density::Array{Float64, 1})
    BondBasedSpecific(make_matrix(S), make_matrix(critical_stretch), density)
end

"""
Bond based material type.
"""
struct BondBasedMaterial <: PeridynamicsMaterial
    name::String
    type::UnitRange{Int64}
    general::GeneralMaterial
    specific::BondBasedSpecific
end

function PeridynamicsMaterial(gen, spc::BondBasedSpecific; name="Default")
    type = minimum(gen.type):maximum(gen.type)
    BondBasedMaterial(name, type, gen, spc)
end


"""
    force_density_T(mat::BondBasedMaterial)

Calculates force density (actually acceleration) for bond based material type.
"""
function force_density_T(y::Array{Float64,2}, mat::BondBasedMaterial)
    types = mat.general.type
    x = mat.general.x 
    intact = mat.general.intact
    family = mat.general.family
    N = size(family, 2)
    M = size(family, 1)
    ARGS = map((i) -> (i, 1:M), 1:N)    

    function cal_force_ij(f, Y, inp)
        M = Y./@_magnitude(Y)
        t = f(inp)
        return t.*M    
    end    
    
    function with_if_cal_force_ij(i, k)
        if intact[k,i]
            j = family[k,i]
            X = [x[1,j]-x[1,i],x[2,j]-x[2,i],x[3,j]-x[3,i]]
            Y = [y[1,j]-y[1,i],y[2,j]-y[2,i],y[3,j]-y[3,i]]
            _X = @_magnitude(X)
            _Y = @_magnitude(Y)
            ext = _Y - _X
            s = ext/_X
            type1 = types[i]- mat.type.start + 1
            type2 = types[j] - mat.type.start + 1
            if s < mat.specific.critical_stretch[type1, type2]
                return cal_force_ij((s) -> mat.specific.bond_stiffness[type1, type2]*s, Y, s)
            else
                intact[k,i] = false
                return [0.0, 0.0, 0.0]
            end            
        else
            return [0.0, 0.0, 0.0]
        end
    end

    # inner_map(i, inds) = pmapreduce((j)-> with_if_cal_force_ij(i,j), +, inds)
    
    inner_map(i, inds) = sum(map((j)-> with_if_cal_force_ij(i,j), inds)) / mat.general.volume[i]
    outer_map(ARGS) = map((x)->inner_map(x[1], x[2]), ARGS)

    return hcat(outer_map(ARGS)...)
end



