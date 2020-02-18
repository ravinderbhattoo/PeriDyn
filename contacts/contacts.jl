abstract type RepulsionModel end

function short_range_repulsion!(y,f,type,RM::RepulsionModel)
    mask1 = type.==RM.pair[1]
    mask2 = type.==RM.pair[2]
    f1 = f[mask1,:]
    f2 = f[mask2,:]
    x1 = y[mask1,:]
    x2 = y[mask2,:]
    for i in 1:size(RM.neighs,1)
        for k in 1:size(RM.neighs,2)
            j = RM.neighs[i,k]
            if j>0
                f1[i,:] .+= -repulsion_force(x1[i,:].-x2[j,:],1,RM)
                f2[j,:] .+= repulsion_force(x1[i,:].-x2[j,:],2,RM)
            end
        end
    end
    f[mask1,:] = f1
    f[mask2,:] = f2
    return nothing
end


function update_repulsive_neighs!(y,type,RM::RepulsionModel)
    x1 = y[type.==RM.pair[1],:]
    x2 = y[type.==RM.pair[2],:]
    family = zeros(Float64,size(x1,1),size(x2,1))
    max_neighs = min(size(x2,1),RM.max_neighs)
    for i in 1:size(x1,1)
        for j in 1:size(x2,1)
            if 1.0e-3<s_magnitude(x2[j,:].-x1[i,:])<RM.distance
                family[i,j] = j
            end
        end
    end
    for i in 1:size(x1,1)
        family[i,:] = sort(family[i,:])
    end
    RM.neighs[:,1:max_neighs] = family[:,end:-1:end-max_neighs+1]
end

#
