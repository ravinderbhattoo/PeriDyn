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
    S = mat.general
    types = mat.general.type
    X = Float64[1.0,0.0,0.0]
    Y = Float64[1.0,0.0,0.0]
    force = zeros(size(S.x)...)
    for i in 1:size(S.x,2)
        for k in 1:size(S.family,1)
            if !S.intact[k,i]
                force[1,i] += 0.0
                force[2,i] += 0.0
                force[3,i] += 0.0
            else
                j = S.family[k,i]::Int64
                X[1],X[2],X[3] = S.x[1,j]-S.x[1,i],S.x[2,j]-S.x[2,i],S.x[3,j]-S.x[3,i]
                Y[1],Y[2],Y[3] = y[1,j]-y[1,i],y[2,j]-y[2,i],y[3,j]-y[3,i]
                s = (_magnitude(Y) - _magnitude(X))/_magnitude(X)
                type1 = types[i]- mat.type.start + 1
                type2 = types[j] - mat.type.start + 1
                if s < mat.specific.critical_stretch[type1, type2]
                    M = Y./_magnitude(Y)
                    t = mat.specific.bond_stiffness[type1, type2]*s
                    force[:,i] .+= t.*M
                else
                    force[1,i] += 0.0
                    force[2,i] += 0.0
                    force[3,i] += 0.0
                    S.intact[k,i] = false
                end
            end
        end
    end
    return force
end



#
