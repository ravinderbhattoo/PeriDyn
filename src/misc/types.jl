struct VeryBigArray{T, N} <: AbstractArray{T, N}
    x :: Array{T, N}
end 

size(m::VeryBigArray, args...) = size(m.x, args...)
getindex(m::VeryBigArray, args...) = getindex(m.x, args...)
setindex!(m::VeryBigArray, args...) = setindex(m.x, args...)
iterate(m::VeryBigArray, args...) = iterate(m.x, args...)
length(m::VeryBigArray, args...) = length(m.x, args...)
similar(m::VeryBigArray, args...) = VeryBigArray(similar(m.x, args...))

