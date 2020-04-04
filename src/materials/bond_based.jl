struct BondBasedSpecific
    bond_stiffness::Float64
end

struct BondBasedMaterial<:PeridynamicsMaterial
    general::GeneralMaterial
    specific::BondBasedSpecific
end

function s_force_density_T(mat::BondBasedMaterial)
    S = mat.general
    force = zeros(size(S.x)...)
    for i in 1:size(S.x,2)
        for j in 1:size(S.family,1)
            if !S.intact[j,i]
                force[1,i] .+= 0.0
                force[2,i] .+= 0.0
                force[3,i] .+= 0.0
            else
            E = [S.family[i,j],i]
            Y = s_X(E...,S.y)
            X = s_X(E...,S.x)
            e = (s_magnitude(Y) - s_magnitude(X))/s_magnitude(X)
                if e<S.critical_stretch
                    M = Y./s_magnitude(Y)
                    t = mat.specific.bond_stiffness*e
                    force[:,i] .+= t.*M
                else
                    force[1,i] .+= 0.0
                    force[2,i] .+= 0.0
                    force[3,i] .+= 0.0
                    S.intact[j,i] = false
                end
            end
        end
    end
    return force
end



#
