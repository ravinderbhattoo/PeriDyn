abstract type RepulsionModel11 end
abstract type RepulsionModel12 end

function short_range_repulsion!(y,f,type,RM::RepulsionModel12)
    mask1 = type.==RM.pair[1]
    mask2 = type.==RM.pair[2]
    f1 = f[:,mask1]
    f2 = f[:,mask2]
    x1 = y[:,mask1]
    x2 = y[:,mask2]
    for i in 1:size(RM.neighs,2)
        for k in 1:size(RM.neighs,1)
            j = RM.neighs[k,i]
            if j>0
                f1[:,i] .+= -repulsion_acc(x1[:,i].-x2[:,j],1,RM)
                f2[:,j] .+= repulsion_acc(x1[:,i].-x2[:,j],2,RM)
            end
            if j==0 break end
        end
    end
    f[:,mask1] = f1
    f[:,mask2] = f2
    return nothing
end


function short_range_repulsion!(y,f,type,RM::RepulsionModel11)
    mask1 = type.==RM.type
    f1 = f[:,mask1]
    x1 = y[:,mask1]
    for i in 1:size(RM.neighs,2)
        for k in 1:size(RM.neighs,1)
            j = RM.neighs[k,i]
            if j>0
                f1[:,i] .+= -repulsion_acc(x1[:,i].-x1[:,j],RM)
                f1[:,j] .+= repulsion_acc(x1[:,i].-x1[:,j],RM)
            end
            if j==0 break end
        end
    end
    f[:,mask1] = f1
    return nothing
end


function update_repulsive_neighs!(y,type,RM::RepulsionModel12)
    x1 = y[:,type.==RM.pair[1]]
    x2 = y[:,type.==RM.pair[2]]
    family = zeros(Float64,size(x2,2),size(x1,2))
    for i in 1:size(x1,2)
        a1,b1,c1 = x1[1,i],x1[2,i],x1[3,i]
        for j in 1:size(x2,2)
            a2,b2,c2 = x2[1,j],x2[2,j],x2[3,j]
            if (a1-a2)*(a1-a2)+(b1-b2)*(b1-b2)+(c1-c2)*(c1-c2)>RM.distance^2
            else
                family[j,i] = j
            end
        end
    end
    family = sort(family,dims=1)
    RM.neighs[1:end,1:end] = family[end:-1:end+1-RM.max_neighs,1:end]
end

function update_repulsive_neighs!(y,type,RM::RepulsionModel11)
    x = y[:,type.==RM.type]
    fill!(RM.neighs,0)
    cells, cell_neighs = get_cells(x,RM.distance)
    for cell_i in 1:length(cells)
        for ca_id in 1:length(cells[cell_i])
            ind = 1
            ca = cells[cell_i][ca_id]
            a1,b1,c1 = x[1,ca],x[2,ca],x[3,ca]
            for neigh_id in 1:length(cell_neighs[cell_i])
                neighs = cell_neighs[cell_i][neigh_id]
                for fa_id in 1:length(cells[neighs])
                    fa = cells[neighs][fa_id]
                    noskip = true
                    for k in 1:size(RM.material.family,1)
                        if fa < RM.material.family[k,ca]
                            break
                        end
                        if fa==RM.material.family[k,ca]
                            if RM.material.intact[k,ca]
                                noskip = false
                            end
                            break
                        end
                    end
                    if noskip
                        a2,b2,c2 = x[1,fa], x[2,fa], x[3,fa]
                        if 1e-4<((a1-a2)*(a1-a2)+(b1-b2)*(b1-b2)+(c1-c2)*(c1-c2))<RM.distance^2
                            RM.neighs[ind,ca] = fa
                            ind += 1
                        end
                    end
                end
            end
        end
    end
end

#
