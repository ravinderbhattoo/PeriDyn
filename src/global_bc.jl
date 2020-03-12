abstract type BoundaryCondition end

struct Move_BC<:BoundaryCondition
    block_ID::Int64
    bool::Array{Bool,1}
    rate::Array{Float64,1}
    start::Array{Float64,2}
end

function apply_bc!(env,BC::Move_BC)
    y1 = env.y[:,env.type.==BC.block_ID]
    y1[1,BC.bool] = BC.start[1,:] .+ env.time_step*env.dt*BC.rate[1]
    y1[2,BC.bool] = BC.start[2,:] .+ env.time_step*env.dt*BC.rate[2]
    y1[3,BC.bool] = BC.start[3,:] .+ env.time_step*env.dt*BC.rate[3]
    env.y[:,env.type.==BC.block_ID] = y1
end
