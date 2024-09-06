function set_device(device)
    device = get_valid_device(device)
    DEVICE[] = device
    log_impinfo("PeriDyn: DEVICE set to $(DEVICE[])")
end

function get_valid_device(x)
    out = x
    if x==:cuda
        if !CUDA.functional()
            out = :cpu
            log_impinfo("PeriDyn: CUDA is not available.")
        end
    else
        out = :cpu
        log_impinfo("PeriDyn: Number of threads = $(Threads.nthreads()).")
    end
    return out
end

function reset_cuda()
    set_device(:cuda)
end

_cudaconvert(item) = if item!=[] CuArray(item) else [] end
_cudaconvert(item::Vector{Any}) = _cudaconvert.(item)
_cudaconvert(item::Union{Number,String,Bool,Nothing,Function}) = item
_cudaconvert(item::Ref{T}) where T = item
_cudaconvert(item::Tuple) = Tuple(_cudaconvert(i) for i in item)
_cudaconvert(item::BitVector) = _cudaconvert(Vector(item))
_cudaconvert(item::UnitRange) = item
_cudaconvert(item::VeryBigArray) = item

function _cudaconvert(x::Vector{Bool})
    range_ = 1:length(x)
    _cudaconvert(range_[x])
end

function _cudaconvert(item::Dict)
    Dict((Symbol(k) => _cudaconvert(v) for (k,v) in item)...)
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

