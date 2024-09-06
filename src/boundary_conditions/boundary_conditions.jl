"""
This module contains definitions for different types of boundary conditions.
"""

export BoundaryCondition, check!, apply_bc!, apply_bc_at0!

"""
    BoundaryCondition

Abstract type for boundary conditions.
"""
abstract type BoundaryCondition end

function Base.show(io::IO, i::BoundaryCondition)
    print(io, getPanel(i))
end

function getPanel(i::BoundaryCondition; ptype=SPanel, width=Term.default_width())
    txt = ""
    for j in fieldnames(typeof(i))
        item = getproperty(i, j)
        # check if iterable
        txt = txt * "$(variable_color(j)): "
        if isa(item, AbstractArray)
            txt = txt * array_repr(item)
        else
            txt = txt * "$(item)\n"
        end
    end
    txt = txt[1:end-1]
    return ptype(txt,
        title="$(typeof(i))",
        justify=:left,
        width=width)
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
- `dims`: Boolean array specifying the affected dimensions.
- `last`: Array of the same type as the state vector specifying the last state of the affected material points.
- `onlyatstart`: Flag indicating if the boundary condition is applied only at the start (default: `false`).
- `xF`: function for updating the position.
- `vF`: function for updating the velocity.
"""
@def general_bc_p begin
    bool::AbstractArray{T, 1} where T
    dims::AbstractArray{Bool, 1}
    last::AbstractArray{T, 2} where T
    onlyatstart::Bool
    xF::Function
    vF::Function
end


###############################################
# apply BCs
###############################################

"""
    apply_bc!(env, BC::BoundaryCondition, on::Symbol)

Apply the specified boundary condition `BC` to the given environment `env` on the specified aspect `on`.

# Arguments
- `env`: The environment to which the boundary condition is applied.
- `BC`: The boundary condition to apply.
- `on::Symbol`: The aspect on which the boundary condition is applied (`:position` or `:velocity`).
"""
function apply_bc!(env, BC::T, on::Symbol) where T<:BoundaryCondition
    @timeit timings "BC $(on)" apply_bc!(env, BC, Val{on})
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
    env.y[BC.dims, BC.bool] .= a
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
    env.v[BC.dims, BC.bool] .= a
    BC.last .= b
end

###############################################
# Unit stripping
###############################################

function uconvert_to(DIMS, BC::T) where T <: BoundaryCondition
    args = (uconvert_to(DIMS, getfield(BC, k)) for k in fieldnames(T))
    return T(args...)
end

function ustrip(BC::T) where T <: BoundaryCondition
    args = (ustrip(getfield(BC, k)) for k in fieldnames(T))
    return T(args...)
end

###############################################
# cuda
###############################################
function _cudaconvert(x::Vector{T}) where T <: BoundaryCondition
    _cudaconvert.(x)
end

function _cudaconvert(x::T) where T <: BoundaryCondition
    T((_cudaconvert(getfield(x, k)) for k in fieldnames(T))...)
end


###############################################
# add BCs
###############################################

include("./FixBC.jl")
include("./ToFroBC.jl")
include("./DeltaScaleBC.jl")
include("./ScaleFixWaitBC.jl")
include("./ContainerBC.jl")


