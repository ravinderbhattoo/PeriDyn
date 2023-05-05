_cudaconvert(item) = if item!=[] CuArray(item) else [] end
_cudaconvert(item::Vector{Any}) = _cudaconvert.(item)
_cudaconvert(item::Union{Number,String,Bool,Nothing,Function}) = item
_cudaconvert(item::Ref{T}) where T = item
_cudaconvert(item::Tuple) = Tuple(_cudaconvert(i) for i in item)
_cudaconvert(item::UnitRange) = item
_cudaconvert(item::BitVector) = _cudaconvert(Vector(item))

function _cudaconvert(item::Dict)
    (; (Symbol(k) => _cudaconvert(v) for (k,v) in item)...)
end

function _cudaconvert2(x::T) where T
    function fn(x, k)
        println(k)
        _cudaconvert(getfield(x, k))
    end
    T(( fn(getfield(x, k)) for k in fieldnames(T))...)
end


function makecuda!(x)
    if DEVICE[]==:cuda
        return _cudaconvert(x)
    end
    x
end

function deviceconvert(x)
    if DEVICE[]==:cuda
        x = _cudaconvert(x)
    end
    x
end

function apply_kernel(fn, args...; Nthreads=256)
    cuargs = CuArray.(args)
    kernel = CUDA.@cuda launch=false fn(cuargs...)
    config = launch_configuration(kernel.fun)
    nthreads = Base.min(Nthreads, config.threads)
    nblocks =  cld(Nthreads, nthreads)
    CUDA.@sync kernel(cuargs...; threads=nthreads, blocks=nblocks)
    cuargs
end

