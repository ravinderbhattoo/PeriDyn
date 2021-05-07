export force_density_T, BondBasedMaterial, BondBasedSpecific

"""
Specific bond based material type.
"""
struct BondBasedSpecific <: SpecificMaterial
    bond_stiffness::Array{Float64,2}
    critical_stretch::Array{Float64, 2}
    density::Array{Float64, 1}
end

function BondBasedSpecific(S::Array{Float64,1}, critical_stretch, density)
    bs = make_matrix(S)
    BondBasedSpecific(bs, make_matrix(critical_stretch), density)
end

"""
Bond based material type.
"""
struct BondBasedMaterial<:PeridynamicsMaterial
    type::UnitRange{Int64}
    general::GeneralMaterial
    specific::BondBasedSpecific
end

function PeridynamicsMaterial(gen, spc::BondBasedSpecific)
    type = minimum(gen.type):maximum(gen.type)
    BondBasedMaterial(type, gen, spc)
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
        for j in 1:size(S.family,1)
            if !S.intact[j,i]
                force[1,i] += 0.0
                force[2,i] += 0.0
                force[3,i] += 0.0
            else
                jj = S.family[j,i]
                X[1],X[2],X[3] = S.x[1,jj]-S.x[1,i],S.x[2,jj]-S.x[2,i],S.x[3,jj]-S.x[3,i]
                Y[1],Y[2],Y[3] = y[1,jj]-y[1,i],y[2,jj]-y[2,i],y[3,jj]-y[3,i]
                e = (magnitude(Y) - magnitude(X))/magnitude(X)
                type1 = types[i]- mat.type.start + 1
                type2 = types[jj] - mat.type.start + 1
                if e<mat.specific.critical_stretch[type1, type2]
                    M = Y./magnitude(Y)
                    t = mat.specific.bond_stiffness[type1, type2]*e
                    force[:,i] .+= t.*M
                else
                    force[1,i] += 0.0
                    force[2,i] += 0.0
                    force[3,i] += 0.0
                    S.intact[j,i] = false
                end
            end
        end
    end
    return force
end



#
