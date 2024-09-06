"""
This module contains definitions for DeltaScaleBC of boundary conditions.
"""

export DeltaScaleBC

"""
    DeltaScaleBC

Struct representing the DeltaScaleBC boundary condition.

# Fields
- `bool`: Boolean array specifying the affected elements.
- `dims`: Boolean array specifying the affected dimensions.
- `last`: Last position of the affected elements.
- `onlyatstart`: Flag indicating if the boundary condition is applied only at the start.
- `xF`: function for updating the position.
- `vF`: function for updating the velocity.
"""
struct DeltaScaleBC <: BoundaryCondition
    @general_bc_p
    # bool::AbstractArray{Bool, 1}
    # dims::AbstractArray{Bool, 1}
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

## Keyword Arguments
- `dims`: Boolean array specifying the affected dimensions (default: `[true, true, true]`).
- `fixpoint`: Reference point used for scaling.
- `onlyatstart`: Flag indicating if the boundary condition is applied only at the start (default: `false`).

# Returns
A DeltaScaleBC object representing the boundary condition.
"""
function DeltaScaleBC(bool, scale, fixpoint; dims=[true, true, true], onlyatstart=false)
    dims = dims[1:SPATIAL_DIMENSIONS_REF[]]
    scale = scale[1:SPATIAL_DIMENSIONS_REF[]]
    fixpoint = fixpoint[1:SPATIAL_DIMENSIONS_REF[]]
    scale = deviceconvert(scale)
    fixpoint = deviceconvert(fixpoint)
    last = 1*zeros(eltype(fixpoint), SPATIAL_DIMENSIONS_REF[], sum(bool))
    xF = (env, BC) -> begin
                y = (BC.last .- reshape(fixpoint, :)) .* (1 .+ reshape(scale, :)*env.dt*env.time_step) .+ reshape(fixpoint, :)
                (y, BC.last)
            end
    vF = (env, BC) -> ( 0.0*unit(eltype(env.v)), BC.last )
    deviceconvert(DeltaScaleBC(bool, dims, last, onlyatstart, xF, vF))
end


"""
    apply_bc_at0!(env, BC::DeltaScaleBC)

Apply the DeltaScaleBC boundary condition at time 0 to the given environment `env`.

# Arguments
- `env`: The environment to which the boundary condition is applied.
- `BC`: The DeltaScaleBC boundary condition to apply.
"""
function apply_bc_at0!(env, BC::DeltaScaleBC)
    BC.last .= env.y[BC.dims, BC.bool]
end

"""
    check!(BC::DeltaScaleBC, env)

Perform a check on the DeltaScaleBC boundary condition.

# Arguments
- `BC`: The DeltaScaleBC boundary condition to check.
- `env`: The environment associated with the boundary condition.
"""
function check!(env, BC::DeltaScaleBC)
    # pass
end