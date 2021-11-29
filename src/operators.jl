"""
This module is supposed to contains standard peridynamic state definitions (We will keep it only if it is useful).
"""

export _O, _I, _magnitude, _ij, _unit

function _O(x::AbstractArray)
    return 0*x
end

function _I(x::AbstractArray)
    return x
end

function _magnitude(x::Array{T,1}) where T
    return norm(x)
end

function _unit(x::Array{T,1}) where T
    return x / norm(x)
end

function _ij(j::Int64,i::Int64,x::Array{T,2}) where T
    return @. x[:,j] - x[:,i]
end

#
