"""
This module is supposed to contains standard peridynamic state definitions (We will keep it only if it is useful).
"""

"""
    s_null(x::Array{Float64,1})

`Peridynamics state: Null operator`
"""
function s_null(x::Array{Float64,1})
    return zeros(Float64,size(x)...)
end


"""
    s_identity(x::Array{Float64,1})

`Peridynamics state: Identity operator`
"""
function s_identity(x::Array{Float64,1})
    return x
end


"""
    s_magnitude(x::Array{Float64,1})

`Peridynamics state: Modulus operator`
"""
function s_magnitude(x::Array{Float64,1})
    return sqrt(x[1]*x[1]+x[2]*x[2]+x[3]*x[3])
end


"""
    s_X(x::Array{Float64,1})

`Peridynamics state: Vector {i,j} operator`
"""
function s_X(j::Int64,i::Int64,x::Array{Float64,2})
    return x[:,j] .- x[:,i]
end

#
