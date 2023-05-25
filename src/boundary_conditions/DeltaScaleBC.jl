"""
This module contains definitions for DeltaScaleBC of boundary conditions.
"""

export DeltaScaleBC

"""
Struct representing the DeltaScaleBC boundary condition.
"""
struct DeltaScaleBC <: BoundaryCondition
    @general_bc_p
    # bool::AbstractArray{Bool, 1}
    # last::AbstractArray{Float64, 2}
    # onlyatstart::Bool
    # xF::Function
    # vF::Function
end

"""
    DeltaScaleBC(bool, scale, fixpoint; onlyatstart=false)

Construct a DeltaScaleBC boundary condition.

# Arguments
- `bool`: Boolean array specifying the affected elements.
- `scale`: Scale factor applied to the elements.
- `fixpoint`: Reference point used for scaling.
- `onlyatstart`: Flag indicating if the boundary condition is applied only at the start (default: `false`).

# Returns
A DeltaScaleBC object representing the boundary condition.
"""
function DeltaScaleBC(bool, scale, fixpoint; onlyatstart=false)
    scale = deviceconvert(scale)
    fixpoint = deviceconvert(fixpoint)
    last = zeros(Float64, SPATIAL_DIMENSIONS_REF[], sum(bool))
    xF = (env, BC) -> begin
                y = (BC.last .- reshape(fixpoint, :)) .* (1 .+ reshape(scale, :)*env.dt*env.time_step) .+ reshape(fixpoint, :)
                (y, BC.last)
            end
    vF = (env, BC) -> ( 0.0, BC.last )
    deviceconvert(DeltaScaleBC(bool, last, onlyatstart, xF, vF))
end


"""
    apply_bc_at0!(env, BC::DeltaScaleBC)

Apply the DeltaScaleBC boundary condition at time 0 to the given environment `env`.

# Arguments
- `env`: The environment to which the boundary condition is applied.
- `BC`: The DeltaScaleBC boundary condition to apply.
"""
function apply_bc_at0!(env, BC::DeltaScaleBC)
    BC.last .= env.y[:, BC.bool]
end

"""
    check!(BC::DeltaScaleBC, env)

Perform a check on the DeltaScaleBC boundary condition.

# Arguments
- `BC`: The DeltaScaleBC boundary condition to check.
- `env`: The environment associated with the boundary condition.
"""
function check!(BC::DeltaScaleBC, env)
    # pass
end