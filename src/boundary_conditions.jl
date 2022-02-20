"""
This modules contains the boundary conditions definitions.
"""

export FixBC, MoveBC, ToFroBC, ScaleBC, ScaleFixBC, ScaleMoveBC, JustScaleBC, apply_bc!

"""
Abstract class for boundary conditions.
"""
abstract type BoundaryCondition end
abstract type BoundaryConditionat0 end

function Base.show(io::IO, i::BoundaryCondition)  
    println(io, typeof(i))
    for j in fieldnames(typeof(i))
        if j in [:bool, :start]
        else
        println(io, j, ": ", getproperty(i, j))
        end
    end
end

function check!(BC::T, env) where T <: BoundaryCondition
    # error("check! method Not implemented for $(T)")
end

struct FixBC<:BoundaryCondition
    bool::Array{Bool,1}
    start::Array{Float64,2}
    onlyatstart::Bool
end

function FixBC(bool; onlyatstart=false)
    FixBC(bool, zeros(Float64, 3, sum(bool)), onlyatstart)
end

function apply_bc_at0!(env, BC::FixBC)
    y1 = env.y[:, BC.bool]
    BC.start[:, :] .= y1
end

"""
    apply_bc!(env,BC::FixBC)

Applies fixed boundary condition to a given material block.
"""
function apply_bc!(env, BC::FixBC)
    y1 = env.y[:, BC.bool]
    v1 = env.v[:, BC.bool]
    y1[:, :] .= BC.start
    v1[:, :] .= 0.0
    env.y[:, BC.bool] = y1
    env.v[:, BC.bool] = v1
end

struct MoveBC<:BoundaryCondition
    bool::Array{Bool,1}
    rate::Array{Float64,1}
    start::Array{Float64,2}
    onlyatstart::Bool
end

function MoveBC(bool, rate; onlyatstart=false)
    MoveBC(bool, rate, zeros(Float64, 3, sum(bool)), onlyatstart)
end

function apply_bc_at0!(env, BC::MoveBC)
    y1 = env.y[:, BC.bool]
    BC.start[:, :] .= y1
end

"""
    apply_bc!(env,BC::MoveBC)

Applies moving boundary condition to a given material block. Material will move with the given rate (constant velocity).
"""
function apply_bc!(env,BC::MoveBC)
    env.y[:, BC.bool] .= BC.start .+ vec(env.time_step*env.dt*BC.rate)
    env.v[:, BC.bool] .= vec(BC.rate)
end

struct ToFroBC<:BoundaryCondition
    bool::Array{Bool,1}
    rate::Array{Float64,1}
    direction::Ref{Int64}
    freq::Int64
    applyafter::Int64
    onlyatstart::Bool
end

function ToFroBC(bool, rate, freq; applyafter=0, onlyatstart=false)
    ToFroBC(bool, rate, Ref(1), freq, applyafter, onlyatstart)
end

function apply_bc_at0!(env, BC::ToFroBC)
end

"""
    apply_bc!(env,BC::ToFroBC)

Applies to-fro boundary condition to a given material block. Material will move with the given rate (constant velocity).
"""
function apply_bc!(env, BC::ToFroBC)
    y1 = env.y[:, BC.bool]
    v1 = env.v[:, BC.bool]
    y1[:, :] .+=  vec(env.dt*BC.direction[]*BC.rate) .- env.dt*v1
    v1[:, :] .= vec(BC.direction[]*BC.rate)
    env.y[:, BC.bool] = y1
    env.v[:, BC.bool] = v1
end


function check!(BC::ToFroBC, env)
    if ( env.time_step > BC.applyafter ) & (env.time_step % BC.freq == 0)
        BC.direction[] = -1*BC.direction[]
    end    
end


struct ScaleBC<:BoundaryCondition
    bool::Array{Bool,1}
    rate::Array{Float64,1}
    FixPoint::Array{Float64,1}
    start::Array{Float64,2}
    onlyatstart::Bool
end

function ScaleBC(bool, rate, FixPoint; onlyatstart=false)
    ScaleBC(bool,rate,FixPoint,zeros(Float64, 3, sum(bool)), onlyatstart)
end

function ScaleBC(bool, rate; onlyatstart=false)
    ScaleBC(bool,rate,[0.0, 0, 0],zeros(Float64, 3, sum(bool)), onlyatstart)
end

function apply_bc_at0!(env, BC::ScaleBC)
    y1 = env.y[:, BC.bool]
    BC.start[:, :] .= y1
end


"""
    apply_bc!(env,BC::ScaleBC)

It scales a given material block about a given fixed point.
"""
function apply_bc!(env,BC::ScaleBC)
    y1 = env.y[:, BC.bool]
    v1 = env.v[:, BC.bool]
    y1[:, :] .= (BC.start[:, :] .- vec(BC.FixPoint)) .* (1 .+ vec(env.time_step*env.dt*BC.rate)) .+ vec(BC.FixPoint)
    v1[:, :] .= (BC.start[:, :] .- vec(BC.FixPoint)) .* (vec(BC.rate))
    env.y[:, BC.bool] = y1
    env.v[:, BC.bool] = v1
end

struct ScaleMoveBC<:BoundaryCondition
    bool::Array{Bool,1}
    rates::Array{Float64,1}
    FixPoint::Array{Float64,1}
    ratem::Array{Float64,1}
    start::Array{Float64,2}
    onlyatstart::Bool
end

function ScaleMoveBC(bool, rates, FixPoint, ratem; onlyatstart=false )
    ScaleMoveBC(bool,rates,FixPoint,ratem,zeros(Float64, 3, sum(bool)), onlyatstart)
end

function ScaleMoveBC(bool, rates, ratem; onlyatstart=false)
    ScaleMoveBC(bool,rates,[0.0, 0, 0],ratem,zeros(Float64, 3, sum(bool)), onlyatstart)
end

function apply_bc_at0!(env, BC::ScaleMoveBC)
    y1 = env.y[:, BC.bool]
    BC.start[:, :] .= y1
end


"""
    apply_bc!(env,BC::ScaleMoveBC)

It scales a given material block about a given fixed point.
"""
function apply_bc!(env, BC::ScaleMoveBC)
    y1 = env.y[:, BC.bool]
    v1 = env.v[:, BC.bool]
    y1[:, :] .= (BC.start[:, :] .- vec(BC.FixPoint)) .* (1 .+ vec(env.time_step*env.dt*BC.rates)) .+ vec(BC.FixPoint) .+ vec(env.time_step*env.dt*BC.ratem)
    v1[:, :] .= (BC.start[:, :] .- vec(BC.FixPoint)) .* (vec(BC.rates)) .+ vec(BC.ratem)
    env.y[:, BC.bool] = y1
    env.v[:, BC.bool] = v1
end

# struct ScaleFixBC<:BoundaryCondition
#     blockID::Int64
#     bool::Array{Bool,1}
#     m::Array{Float64,1}
#     Fixbool::Array{Bool,1}
#     Fix::Array{Float64,2}
#     LastTimeStep::Array{Int64,1}
# end

# """
#     ScaleFixBC(blockID, bool, m, Fixbool)

# """
# function ScaleFixBC(blockID, bool, m, Fixbool)
#     Fix = zeros(Float64,3,sum(Fixbool))
#     LastTimeStep = [-1]
#     ScaleFixBC(blockID,bool,m,Fixbool,Fix,LastTimeStep)
# end

# """
#     apply_bc!(env,BC::ScaleFixBC)

# It scales as well as fix a given material block. (like scale and fix edges, helpful in tensile quasi-static simulations.)
# """
# function apply_bc!(env,BC::ScaleFixBC)
#     if BC.LastTimeStep[1]==env.time_step
#     else
#         BC.LastTimeStep[1]=env.time_step
#         m_ = (env.dt*BC.m) ./ ( 1 .+ env.dt*BC.m .* (env.time_step .- 1))
#         y1 = env.y[:,env.type.==BC.blockID]
#         y1[1,BC.bool] = (y1[1,BC.bool]) .* (1+m_[1])
#         y1[2,BC.bool] = (y1[2,BC.bool]) .* (1+m_[2])
#         y1[3,BC.bool] = (y1[3,BC.bool]) .* (1+m_[3])
#         env.y[:,env.type.==BC.blockID] = y1
#         BC.Fix[:,:] = y1[:,BC.Fixbool]
#     end
#     y1 = env.y[:,env.type.==BC.blockID]
#     y1[1,BC.Fixbool] = BC.Fix[1,:]
#     y1[2,BC.Fixbool] = BC.Fix[2,:]
#     y1[3,BC.Fixbool] = BC.Fix[3,:]
#     env.y[:,env.type.==BC.blockID] = y1
# end
