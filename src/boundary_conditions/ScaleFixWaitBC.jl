"""
Module containing definitions for ScaleFixWaitBC boundary conditions.
"""

export ScaleFixWaitBC

"""
Structure representing a ScaleFixWaitBC boundary condition.
"""
struct ScaleFixWaitBC <: BoundaryCondition
    @general_bc_p
    checkF::Function
end

"""
    ScaleFixWaitBC(bool, scale, fixpoint, wait, scalebool; applyafter=0, onlyatstart=false)

Construct a ScaleFixWaitBC boundary condition.

# Arguments
- `bool`: Boolean array specifying the affected elements.
- `scale`: Scale factor for the elements.
- `fixpoint`: Fix point for the elements.
- `wait`: Number of time steps to wait before applying the condition.
- `scalebool`: Boolean array specifying the elements to be scaled.
- `applyafter`: Number of time steps after which the condition is applied (default: 0).
- `onlyatstart`: Boolean indicating whether the condition is applied only at the start (default: false).

# Returns
- An instance of ScaleFixWaitBC boundary condition.
"""
function ScaleFixWaitBC(bool, scale, fixpoint, wait, scalebool; applyafter=0, onlyatstart=false)
    fixpoint = deviceconvert(fixpoint)
    scale = deviceconvert(scale)
    last = zeros(Float64, SPATIAL_DIMENSIONS_REF[], sum(bool))
    xF = (env, BC) -> (BC.last, BC.last)
    vF = (env, BC) -> (0.0, BC.last)
    checkF = (env, BC) -> begin
        if (env.time_step % wait == 0)
            y = env.y[:, scalebool]
            y .= (y .- reshape(fixpoint, :)) .* sqrt.(1 .+ 2*reshape(scale, :)*env.dt) .+ reshape(fixpoint, :)
            env.y[:, scalebool] .= y
            BC.last .= env.y[:, BC.bool]
        end
    end
    deviceconvert(ScaleFixWaitBC(bool, last, onlyatstart, xF, vF, checkF))
end

"""
    apply_bc_at0!(env, BC::ScaleFixWaitBC)

Apply the ScaleFixWaitBC boundary condition at time 0.

# Arguments
- `env`: Environment in which the condition is applied.
- `BC`: Instance of ScaleFixWaitBC boundary condition.
"""
function apply_bc_at0!(env, BC::ScaleFixWaitBC)
    BC.last .= env.y[:, BC.bool]
end

"""
    check!(BC::ScaleFixWaitBC, env)

Check if the ScaleFixWaitBC boundary condition needs to be applied.

# Arguments
- `BC`: Instance of ScaleFixWaitBC boundary condition.
- `env`: Environment in which the condition is applied.
"""
function check!(BC::ScaleFixWaitBC, env)
    BC.checkF(env, BC)
end

