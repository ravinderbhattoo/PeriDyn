function s_null(x::Array{Float64})
    return zeros(Float64,size(x)...)
end


function s_identity(x::Array{Float64})
    return x
end


function s_magnitude(x::Array{Float64})
    return sqrt(sum(x.*x))
end


function s_X(E,x)
    return x[E[2],:] .- x[E[1],:]
end

function s_Y(E,y)
    return y[E[2],:] .- y[E[1],:]
end
#
