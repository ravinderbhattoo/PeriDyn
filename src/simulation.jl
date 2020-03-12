struct AbstractEnv
    type::Array{Int64,1}
    y::Array{Float64,2}
    v::Array{Float64,2}
    f::Array{Float64,2}
    time::Float64
    dt::Float64
    N::Int64
    neighs::Array{Int64,2}
    short_range_repulsion::Any
    material_blocks::Any
end


function Env(materials,short_range_repulsion,time,dt)
    type = materials[1].type*ones(Int64,size(materials[1].general.y,1))
    for i in 2:size(materials,1)
        type = vcat(type,materials[i].type*ones(Int64,size(materials[i].general.y,1)))
    end
    y = materials[1].general.y
    for mat in materials[2:end]
        y = vcat(y,mat.general.y)
    end
    v = materials[1].general.velocity
    for mat in materials[2:end]
        v = vcat(v,mat.general.velocity)
    end
    N = Int64(time/dt)
    return AbstractEnv(type,y,v,0v,time,dt,N,zeros(2,2),short_range_repulsion,materials)
end


#
