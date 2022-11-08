function _makecuda!(env::T) where T
    for name in fieldnames(T)
        item = getfield(env, name)
        if typeof(item) <: AbstractArray
            try
                setfield!(env, name, CuArray(item))
            catch
            end
        end
    end
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
