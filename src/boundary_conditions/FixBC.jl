"""
This module contains definitions for FixBC boundary conditions.
"""

export FixBC

"""
Struct representing the FixBC boundary condition.
"""
struct FixBC <: BoundaryCondition
    @general_bc_p
    # bool::AbstractArray{Bool, 1}
    # last::AbstractArray{Float64, 2}
    # onlyatstart::Bool
    # xF::Function
    # vF::Function
end


"""
    FixBC(bool; onlyatstart=false)

Construct a FixBC boundary condition.

# Arguments
- `bool`: Boolean array specifying the affected elements.
- `onlyatstart`: Flag indicating if the boundary condition is applied only at the start (default: `false`).

# Returns
A FixBC object representing the boundary condition.
"""
function FixBC(bool; onlyatstart=false)
    last = zeros(Float64, SPATIAL_DIMENSIONS_REF[], sum(bool))
    xF = (env, BC) -> (BC.last, BC.last)
    vF = (env, BC) -> (0.0, BC.last)
    deviceconvert(FixBC(bool, last, onlyatstart, xF, vF))
end

"""
    apply_bc_at0!(env, BC::FixBC)

Apply the FixBC boundary condition at time 0 to the given environment `env`.

# Arguments
- `env`: The environment to which the boundary condition is applied.
- `BC`: The FixBC boundary condition to apply.
"""
function apply_bc_at0!(env, BC::FixBC)
    BC.last .= env.y[:, BC.bool]
end

