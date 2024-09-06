"""
Module containing definitions for ScaleFixWaitBC boundary conditions.
"""

export ScaleFixWaitBC

"""
    ScaleFixWaitBC

Structure representing a ScaleFixWaitBC boundary condition.

# Fields
- `bool`: Boolean array specifying the affected elements.
- `dims`: Boolean array specifying the affected dimensions.
- `last`: Last position of the affected elements.
- `onlyatstart`: Flag indicating if the boundary condition is applied only at the start.
- `xF`: function for updating the position.
- `vF`: function for updating the velocity.
- `checkF`: Function for checking if the boundary condition needs to be applied.
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

# Keyword Arguments
- `dims`: Boolean array specifying the affected dimensions (default: `[true, true, true]`).
- `applyafter`: Number of time steps after which the condition is applied (default: 0).
- `onlyatstart`: Boolean indicating whether the condition is applied only at the start (default: false).

# Returns
- An instance of ScaleFixWaitBC boundary condition.
"""
function ScaleFixWaitBC(bool, scale, fixpoint, wait, scalebool; dims=[true, true, true], applyafter=0, onlyatstart=false)
    dims = dims[1:SPATIAL_DIMENSIONS_REF[]]
    fixpoint = fixpoint[1:SPATIAL_DIMENSIONS_REF[]]
    scale = scale[1:SPATIAL_DIMENSIONS_REF[]]
    fixpoint = deviceconvert(fixpoint)
    scale = deviceconvert(scale)
    last = 1*zeros(Float64, SPATIAL_DIMENSIONS_REF[], sum(bool))
    xF = (env, BC) -> (BC.last, BC.last)
    vF = (env, BC) -> (0.0, BC.last)
    checkF = (env, BC) -> begin
        if (env.time_step % wait == 0)
            y = env.y[BC.dims, scalebool]
            y .= (y .- reshape(fixpoint, :)) .* sqrt.(1 .+ 2*reshape(scale, :)*env.dt) .+ reshape(fixpoint, :)
            env.y[BC.dims, scalebool] .= y
            BC.last .= env.y[BC.dims, BC.bool]
        end
    end
    deviceconvert(ScaleFixWaitBC(bool, dims, last, onlyatstart, xF, vF, checkF))
end

"""
    apply_bc_at0!(env, BC::ScaleFixWaitBC)

Apply the ScaleFixWaitBC boundary condition at time 0.

# Arguments
- `env`: Environment in which the condition is applied.
- `BC`: Instance of ScaleFixWaitBC boundary condition.
"""
function apply_bc_at0!(env, BC::ScaleFixWaitBC)
    BC.last .= env.y[BC.dims, BC.bool]
end

"""
    check!(BC::ScaleFixWaitBC, env)

Check if the ScaleFixWaitBC boundary condition needs to be applied.

# Arguments
- `BC`: Instance of ScaleFixWaitBC boundary condition.
- `env`: Environment in which the condition is applied.
"""
function check!(env, BC::ScaleFixWaitBC)
    BC.checkF(env, BC)
end

