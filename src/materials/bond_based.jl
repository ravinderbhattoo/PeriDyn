struct BondBasedSpecific
    bond_stiffness::Float64
end

struct BondBasedMaterial<:PeridynamicsMaterial
    general::Any
    specific::Any
end

function s_force_density_T(mat::BondBasedMaterial)
    S = mat.general
    force = zeros(size(S.x)...)
    for i in 1:size(S.x,1)
        for j in 1:size(S.family,2)
            if !S.damage[i,j]
                force[i,:] .+= [0,0,0]
            else
            E = [i,S.family[i,j]]
            Y = s_Y(E,S.y)
            X = s_X(E,S.x)
            e = (s_magnitude(Y) - s_magnitude(X))/s_magnitude(X)
                if e<S.critical_stretch
                    M = Y./s_magnitude(Y)
                    t = mat.specific.bond_stiffness*e
                    force[i,:] .+= t.*M
                else
                    force[i,:] .+= [0,0,0]
                    S.damage[i,j] = true
                end
            end
        end
    end
    return force
end



#
