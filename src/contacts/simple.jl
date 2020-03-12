struct SimpleRepulsionModel<:RepulsionModel
    pair::Array{Int64,1}
    densities::Array{Float64,1}
    alpha::Float64
    epsilon::Float64
    equi_dist::Float64
    neighs::Array{Int64,2}
    distance::Float64
    max_neighs::Int64
end


function SimpleRepulsionModel(alpha::Float64,epsilon::Float64, mat1::PeridynamicsMaterial,mat2::PeridynamicsMaterial;distanceX=5,max_neighs=50)
    p_size = (mat1.general.particle_size+mat2.general.particle_size)/2
    neighs = zeros(size(mat1.general.x,1),max_neighs)
    SimpleRepulsionModel([mat1.type,mat2.type],[mat1.general.density,mat2.general.density],alpha,epsilon,p_size,neighs,p_size*distanceX,max_neighs)
end


function repulsion_force(dr,den_i,RepMod::SimpleRepulsionModel)
    del_x = RepMod.equi_dist - s_magnitude(dr)
    if del_x<0
        return zeros(size(dr)...)
    else
        return -(RepMod.epsilon*del_x^(RepMod.alpha-1)).*dr/s_magnitude(dr)/RepMod.densities[den_i]
    end
end
