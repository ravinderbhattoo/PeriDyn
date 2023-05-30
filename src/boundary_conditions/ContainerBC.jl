"""
This module contains definitions for ContainerBC boundary conditions.
"""

export ContainerBC

"""
    ContainerBC

Struct representing the ContainerBC boundary condition.

# Fields
- `bool`: Boolean array specifying the affected elements.
- `last`: Last position of the affected elements.
- `onlyatstart`: Flag indicating if the boundary condition is applied only at the start.
- `xF`: function for updating the position.
- `vF`: function for updating the velocity.
"""
struct ContainerBC <: BoundaryCondition
    @general_bc_p
    # bool::AbstractArray{Bool, 1}
    # last::AbstractArray{Float64, 2}
    # onlyatstart::Bool
    # xF::Function
    # vF::Function
end


"""
    ContainerBC(bool; limits=nothing, onlyatstart=false)

Construct a ContainerBC boundary condition.

# Arguments
- `bool`: Boolean array specifying the affected elements.
- `limits`: Limits of the container (default: `nothing`).
- `onlyatstart`: Flag indicating if the boundary condition is applied only at the start (default: `false`).

# Returns
A ContainerBC object representing the boundary condition.
"""
function ContainerBC(bool; limits=nothing, onlyatstart=false)
    last = [nothing;;] # Abstract Matrix of last positions
    if isa(limit, Nothing)
        xF = (env, BC) -> begin
                        nothing
                    end
        vF = (env, BC) -> begin
                        nothing
                    end
        return deviceconvert(ContainerBC(bool, last, onlyatstart, xF, vF))
    else
        xF = (env, BC) -> begin
                        y = @view env.y[:, BC.bool]
                        y, BC.last
                    end
        vF = (env, BC) -> begin
                        v = @view env.v[:, BC.bool]
                        y = @view env.y[:, BC.bool]
                        mask = reshape(any(y .< limits[:, 1] .|| y .> limits[:, 2], dims=1), :)
                        v[:, mask] .= 0.0
                        v, BC.last
                    end
        return deviceconvert(ContainerBC(bool, last, onlyatstart, xF, vF))
    end
end

"""
    apply_bc_at0!(env, BC::ContainerBC)

Apply the ContainerBC boundary condition at time 0 to the given environment `env`.

# Arguments
- `env`: The environment to which the boundary condition is applied.
- `BC`: The ContainerBC boundary condition to apply.
"""
function apply_bc_at0!(env, BC::ContainerBC)
    nothing
end

