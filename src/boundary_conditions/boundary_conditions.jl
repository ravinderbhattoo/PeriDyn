"""
This modules contains the boundary conditions definitions.
"""

export apply_bc!

"""
Abstract class for boundary conditions.
"""
abstract type BoundaryCondition end
abstract type BoundaryConditionat0 end

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

function check!(BC::T, env) where T <: BoundaryCondition
    # error("check! method Not implemented for $(T)")
end

@def general_bc_p begin
    bool::AbstractArray{T, 1} where T
    last::AbstractArray{T, 2} where T
    onlyatstart::Bool
    xF::Function
    vF::Function
end


function apply_bc!(env, BC::T, on::Symbol) where T<:BoundaryCondition
    apply_bc!(env, BC, Val{on})
end

"""
    apply_bc!(env,BC::BoundaryCondition, :position)

Applies general boundary condition to a given material block.
"""
function apply_bc!(env, BC::T, ::Type{Val{:position}}) where T <: BoundaryCondition
    a, b = BC.xF(env, BC)
    env.y[:, BC.bool] .= a
    BC.last .= b
end

"""
    apply_bc!(env,BC::BoundaryCondition, :velocity)

Applies general boundary condition to a given material block.
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


