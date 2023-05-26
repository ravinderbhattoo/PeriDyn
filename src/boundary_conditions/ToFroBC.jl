"""
This module contains definitions for ToFroBC and MoveBC of boundary conditions.
"""

export ToFroBC, MoveBC

"""
    ToFroBC

Struct representing the ToFroBC boundary condition.

# Fields
- `bool`: Boolean array specifying the affected elements.
- `last`: Last position of the affected elements.
- `onlyatstart`: Flag indicating if the boundary condition is applied only at the start.
- `xF`: Function for updating the velocity.
- `vF`: Function for updating the position.
- `direction`: Direction of movement.
- `freq`: Frequency at which the direction of movement changes.
- `applyafter`: Number of steps after which the frequency is applied.
"""
struct ToFroBC <: BoundaryCondition
    @general_bc_p
    direction::Ref{Int64}
    freq::Union{Int64,Float64}
    applyafter::Int64  # apply freq after this many steps
end


"""
    ToFroBC(bool, rate, freq; applyafter=0, onlyatstart=false)

Construct a ToFroBC boundary condition.

# Arguments
- `bool`: Boolean array specifying the affected elements.
- `rate`: Rate at which the elements move.
- `freq`: Frequency at which the direction of movement changes.
- `applyafter`: Number of steps after which the frequency is applied (default: 0).
- `onlyatstart`: Flag indicating if the boundary condition is applied only at the start (default: `false`).

# Returns
A ToFroBC object representing the boundary condition.
"""
function ToFroBC(bool, rate, freq; applyafter=0, onlyatstart=false)
    rate = deviceconvert(rate)
    last = zeros(Float64, SPATIAL_DIMENSIONS_REF[], sum(bool))
    xF = (env, BC) -> begin
        y = BC.last .+ reshape(env.dt * BC.direction[] * rate, :)
        (y, y)
    end
    vF = (env, BC) -> ( reshape(BC.direction[] * rate, :), BC.last)
    deviceconvert(ToFroBC(bool, last, onlyatstart, xF, vF, Ref(1), freq, applyafter))
end


"""
    MoveBC(bool, rate; kwargs...)

Create a MoveBC boundary condition.

# Arguments
- `bool`: Boolean array specifying the affected elements.
- `rate`: Rate at which the elements move.
- `kwargs`: Additional keyword arguments passed to `ToFroBC`.

# Returns
A MoveBC object representing the boundary condition.
All MoveBC objects are ToFroBC objects with frequency Inf.
"""
function MoveBC(bool, rate; kwargs...)
    ToFroBC(bool, rate, Inf; kwargs...)
end



"""
    apply_bc_at0!(env, BC::ToFroBC)

Apply the ToFroBC boundary condition at time 0 to the given environment `env`.

# Arguments
- `env`: The environment to which the boundary condition is applied.
- `BC`: The ToFroBC boundary condition to apply.
"""
function apply_bc_at0!(env, BC::ToFroBC)
    BC.last .= env.y[:, BC.bool]
end

"""
    check!(BC::ToFroBC, env)

Perform a check on the ToFroBC boundary condition.

# Arguments
- `BC`: The ToFroBC boundary condition to check.
- `env`: The environment associated with the boundary condition.
"""
function check!(env, BC::ToFroBC)
    if ( env.time_step > BC.applyafter ) & (env.time_step % BC.freq == 0)
        BC.direction[] = -1*BC.direction[]
    end
end

