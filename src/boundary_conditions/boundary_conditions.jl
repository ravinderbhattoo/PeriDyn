"""
This module contains definitions for different types of boundary conditions.
"""

export BoundaryCondition, check!, apply_bc!, apply_bc_at0!

"""
    BoundaryCondition

Abstract type for boundary conditions.
"""
abstract type BoundaryCondition end

function _cudaconvert(x::Vector{T}) where T <: BoundaryCondition
    _cudaconvert.(x)
end

function _cudaconvert(x::T) where T <: BoundaryCondition
    T((_cudaconvert(getfield(x, k)) for k in fieldnames(T))...)
end

function Base.show(io::IO, i::BoundaryCondition)
    println(io, typeof(i))
    for j in fieldnames(typeof(i))
        item = getproperty(i, j)
        if isa(item, AbstractArray)
            println(io, j, ": $(typeof(item)) $(size(item))")
        else
            println(io, j, ": ", item)
        end
    end
end

"""
    check!(env, BC::BoundaryCondition)

Check if the boundary condition `BC` has changed and updates the `BC`. Used for dynamic boundary conditions.
"""
function check!(env, BC::T) where T <: BoundaryCondition
    # error("check! method Not implemented for $(T)")
    log_detail("check! method not implemented for $(T)")
end

"""
    apply_bc_at0!(env, BC::BoundaryCondition)

Apply the boundary condition `BC` to the given environment `env` at the start of the simulation (t=0).
"""
function apply_bc_at0!(env, BC::T) where T <: BoundaryCondition
    # error("check! method Not implemented for $(T)")
    log_detail("apply_bc_at0! method not implemented for $(T)")
end

"""
    general_bc_p()

Macro for defining the common fields of a general boundary condition.

The following fields are defined:
- `bool`: Boolean array specifying the affected material points.
- `last`: Array of the same type as the state vector specifying the last state of the affected material points.
- `onlyatstart`: Flag indicating if the boundary condition is applied only at the start (default: `false`).
- `xF`: Function for updating the velocity.
- `vF`: Function for updating the position.
"""
@def general_bc_p begin
    bool::AbstractArray{T, 1} where T
    last::AbstractArray{T, 2} where T
    onlyatstart::Bool
    xF::Function
    vF::Function
end


"""
    apply_bc!(env, BC::BoundaryCondition, on::Symbol)

Apply the specified boundary condition `BC` to the given environment `env` on the specified aspect `on`.

# Arguments
- `env`: The environment to which the boundary condition is applied.
- `BC`: The boundary condition to apply.
- `on::Symbol`: The aspect on which the boundary condition is applied (`:position` or `:velocity`).
"""
function apply_bc!(env, BC::T, on::Symbol) where T<:BoundaryCondition
    apply_bc!(env, BC, Val{on})
end

"""
    apply_bc!(env, BC::BoundaryCondition, ::Type{Val{:position}})

Apply the general boundary condition `BC` to the position aspect of the given environment `env`.

# Arguments
- `env`: The environment to which the boundary condition is applied.
- `BC`: The general boundary condition to apply.
"""
function apply_bc!(env, BC::T, ::Type{Val{:position}}) where T <: BoundaryCondition
    a, b = BC.xF(env, BC)
    env.y[:, BC.bool] .= a
    BC.last .= b
end

"""
    apply_bc!(env, BC::BoundaryCondition, ::Type{Val{:velocity}})

Apply the general boundary condition `BC` to the velocity aspect of the given environment `env`.

# Arguments
- `env`: The environment to which the boundary condition is applied.
- `BC`: The general boundary condition to apply.
"""
function apply_bc!(env, BC::T, ::Type{Val{:velocity}}) where T <: BoundaryCondition
    a, b = BC.vF(env, BC)
    env.v[:, BC.bool] .= a
    BC.last .= b
end


include("./FixBC.jl")
include("./ToFroBC.jl")
include("./DeltaScaleBC.jl")
include("./ScaleFixWaitBC.jl")


