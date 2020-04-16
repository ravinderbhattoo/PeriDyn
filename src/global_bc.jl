"""
This modules contains the boundary conditions definitions.
"""

"""
Abstract class for boundary conditions.
"""
abstract type BoundaryCondition end


struct FixBC<:BoundaryCondition
    blockID::Int64
    bool::Array{Bool,1}
    start::Array{Float64,2}
end

"""
    apply_bc!(env,BC::FixBC)

Applies fixed boundary condition to a given material block.
"""
function apply_bc!(env,BC::FixBC)
    y1 = env.y[:,env.type.==BC.blockID]
    y1[1,BC.bool] = BC.start[1,:]
    y1[2,BC.bool] = BC.start[2,:]
    y1[3,BC.bool] = BC.start[3,:]
    env.y[:,env.type.==BC.blockID] = y1
end

struct MoveBC<:BoundaryCondition
    blockID::Int64
    bool::Array{Bool,1}
    rate::Array{Float64,1}
    start::Array{Float64,2}
end

"""
    apply_bc!(env,BC::MoveBC)

Applies moving boundary condition to a given material block. Material will move with the given rate (constant velocity).
"""
function apply_bc!(env,BC::MoveBC)
    y1 = env.y[:,env.type.==BC.blockID]
    y1[1,BC.bool] = BC.start[1,:] .+ env.time_step*env.dt*BC.rate[1]
    y1[2,BC.bool] = BC.start[2,:] .+ env.time_step*env.dt*BC.rate[2]
    y1[3,BC.bool] = BC.start[3,:] .+ env.time_step*env.dt*BC.rate[3]
    env.y[:,env.type.==BC.blockID] = y1
end

struct ScaleBC<:BoundaryCondition
    blockID::Int64
    bool::Array{Bool,1}
    rate::Array{Float64,1}
    FixPoint::Array{Float64,1}
    LastTimeStep::Array{Int64,1}
end

function ScaleBC(blockID,bool,rate,FixPoint)
    LastTimeStep = [-1]
    ScaleBC(blockID,bool,rate,FixPoint,LastTimeStep)
end

"""
    apply_bc!(env,BC::ScaleBC)

It scales a given material block about a given fixed point.
"""
function apply_bc!(env,BC::ScaleBC)
    if BC.LastTimeStep[1]==env.time_step
    else
        BC.LastTimeStep[1]=env.time_step
        y1 = env.y[:,env.type.==BC.blockID]
        y1[1,BC.bool] = (y1[1,BC.bool] .- BC.FixPoint[1]) * env.dt*BC.rate[1] .+ BC.FixPoint[1]
        y1[2,BC.bool] = (y1[2,BC.bool] .- BC.FixPoint[2]) * env.dt*BC.rate[2] .+ BC.FixPoint[2]
        y1[3,BC.bool] = (y1[3,BC.bool] .- BC.FixPoint[3]) * env.dt*BC.rate[3] .+ BC.FixPoint[3]
        env.y[:,env.type.==BC.blockID] = y1
    end
end

struct JustScaleBC<:BoundaryCondition
    blockID::Int64
    bool::Array{Bool,1}
    rate::Array{Float64,1}
    start::Array{Float64,2}
    FixPoint::Array{Float64,1}
    LastTimeStep::Array{Int64,1}
end

function JustScaleBC(blockID,bool,rate,start,FixPoint)
    LastTimeStep = [-1]
    JustScaleBC(blockID,bool,rate,start,FixPoint,LastTimeStep)
end


"""
    apply_bc!(env,BC::JustScaleBC)

It scales a given material block about a given fixed point.
"""
function apply_bc!(env,BC::JustScaleBC)
    if BC.LastTimeStep[1]==env.time_step
    else
        BC.LastTimeStep[1]=env.time_step
        y1 = env.y[:,env.type.==BC.blockID]
        y1[1,BC.bool] = (BC.start[1,:] .- BC.FixPoint[1]) .* env.dt*BC.rate[1] .+ BC.FixPoint[1]
        y1[2,BC.bool] = (BC.start[2,:] .- BC.FixPoint[2]) .* env.dt*BC.rate[2] .+ BC.FixPoint[2]
        y1[3,BC.bool] = (BC.start[3,:] .- BC.FixPoint[3]) .* env.dt*BC.rate[3] .+ BC.FixPoint[3]
        env.y[:,env.type.==BC.blockID] = y1
    end
end

struct ScaleFixBC<:BoundaryCondition
    blockID::Int64
    bool::Array{Bool,1}
    m::Array{Float64,1}
    Fixbool::Array{Bool,1}
    Fix::Array{Float64,2}
    LastTimeStep::Array{Int64,1}
end

function ScaleFixBC(blockID,bool,m,Fixbool)
    Fix = zeros(Float64,3,sum(Fixbool))
    LastTimeStep = [-1]
    ScaleFixBC(blockID,bool,m,Fixbool,Fix,LastTimeStep)
end

"""
    apply_bc!(env,BC::ScaleFixBC)
    
It scales as well as fix a given material block. (like scale and fix edges, helpful in tensile quasi-static simulations.)
"""
function apply_bc!(env,BC::ScaleFixBC)
    if BC.LastTimeStep[1]==env.time_step
    else
        BC.LastTimeStep[1]=env.time_step
        m_ = (env.dt*BC.m) ./ ( 1 .+ env.dt*BC.m .* (env.time_step .- 1))
        y1 = env.y[:,env.type.==BC.blockID]
        y1[1,BC.bool] = (y1[1,BC.bool]) .* (1+m_[1])
        y1[2,BC.bool] = (y1[2,BC.bool]) .* (1+m_[2])
        y1[3,BC.bool] = (y1[3,BC.bool]) .* (1+m_[3])
        env.y[:,env.type.==BC.blockID] = y1
        BC.Fix[:,:] = y1[:,BC.Fixbool]
    end
    y1 = env.y[:,env.type.==BC.blockID]
    y1[1,BC.Fixbool] = BC.Fix[1,:]
    y1[2,BC.Fixbool] = BC.Fix[2,:]
    y1[3,BC.Fixbool] = BC.Fix[3,:]
    env.y[:,env.type.==BC.blockID] = y1
end
