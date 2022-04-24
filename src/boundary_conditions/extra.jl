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
