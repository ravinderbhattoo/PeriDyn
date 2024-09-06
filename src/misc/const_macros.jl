export TIMEIT_REF, CHECK_NAN, SPATIAL_DIMENSIONS_REF, MULTI_THREAD_REF
export check_nan, time_code, sinfo, swarning, sdetail, simpinfo
export _magnitude, _ij
export refresh, set_multi_threading

# Constants and variables
const timings = TimerOutput(); disable_timer!(timings)
SOLVERS =  Dict()
const PDBlockID = Ref{Int64}(1)
const DEVICE = Ref{Symbol}(:cpu)
const TIMEIT_REF = Ref(false)
const CHECK_NAN = Ref(true)
const MULTI_THREAD_REF = Ref(true)
const SPATIAL_DIMENSIONS_REF = Ref(3)

"""
    def(name, defination)

Define a macro.

# Arguments
- `name`: Symbol, the name of the macro.
- `defination`: Any, the defination of the macro.

"""
macro def(name, definition)
    return quote
        macro $(esc(name))()
            esc($(Expr(:quote, definition)))
        end
    end
end

"""
    set_multi_threading(x)

Set the multi threading.

# Arguments
- `x`: Bool, the multi threading.

"""
function set_multi_threading(x)
    MULTI_THREAD_REF[] = x
    refresh()
    println("Multi threading: ", x)
end


"""
    check_nan(var, name)

Check if there is nan or inf in the variable.

# Arguments
- `var`: Any, the variable.
- `name`: String, the name of the variable.

"""
macro check_nan(var, name)
    if CHECK_NAN[]
        quote
            local var = $(esc(var))
            local name = $(esc(name))
            if any(isnan, var)
                error("nan in $name")
            end
            if any(isinf, var)
                error("inf in $name")
            end
        end
    else
    end
end


"""
    time_code(ex, name)

Time the execution of the expression.

# Arguments
- `ex`: Any, the expression.
- `name`: String, the name of the expression.

"""
macro time_code(ex, name)
    if TIMEIT_REF[]
        quote
            local name = $(esc(name))
            local val = @time $(esc(ex))
            println("Time taken for $name: ")
            println("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
            println("#############################################################")
            val
        end
    else
        quote
            $(esc(ex))
        end
    end
end

"""
    _magnitude(a)

Get the magnitude of a vector depending on the spatial dimensions.

# Arguments
- `a`: Vector, the vector.

"""
macro _magnitude(a)
    if PeriDyn.SPATIAL_DIMENSIONS_REF[]==3
        quote
            sqrt($(esc(a))[1]*$(esc(a))[1] + $(esc(a))[2]*$(esc(a))[2] + $(esc(a))[3]*$(esc(a))[3])
        end
    else
        quote
            sqrt($(esc(a))[1]*$(esc(a))[1] + $(esc(a))[2]*$(esc(a))[2])
        end
    end
end

"""
    get_magnitude(a)

Get the magnitude of a vector depending on the spatial dimensions.

# Arguments
- `a`: Vector, the vector.

"""
get_magnitude(a) = @_magnitude(a)


"""
    _ij(j, i, x)

Get the difference of two vectors depending on the spatial dimensions.

``X_j - X_i``

# Arguments
- `j`: Int64, the index of the first vector.
- `i`: Int64, the index of the second vector.
- `x`: Matrix, the matrix of the vectors.

"""
macro _ij(j, i, x)
    if PeriDyn.SPATIAL_DIMENSIONS_REF[]==3
        quote
            ($(esc(x))[1,$(esc(j))] - $(esc(x))[1,$(esc(i))],
                    $(esc(x))[2,$(esc(j))] - $(esc(x))[2,$(esc(i))],
                    $(esc(x))[3,$(esc(j))] - $(esc(x))[3,$(esc(i))])
        end
    else
        quote
            ($(esc(x))[1,$(esc(j))] - $(esc(x))[1,$(esc(i))],
                    $(esc(x))[2,$(esc(j))] - $(esc(x))[2,$(esc(i))])
        end
    end
end

"""
    get_ij(j, i, x)

Get the difference of two vectors depending on the spatial dimensions.

``X_j - X_i``

# Arguments
- `j`: Int64, the index of the first vector.
- `i`: Int64, the index of the second vector.
- `x`: Matrix, the matrix of the vectors.

"""
get_ij(j, i, a) = @_ij(j, i, a)

"""
    applyops(x)

Apply the operations to the vectors depending on the spatial dimensions.

# Arguments
- `x`: String, the operations.

"""
macro applyops(x)
    local a,b,c = replace(x, "#"=>1), replace(x, "#"=>2), replace(x, "#"=>3)
    if PeriDyn.SPATIAL_DIMENSIONS_REF[]==3
        return quote
            Meta.parse($a)
            Meta.parse($b)
            Meta.parse($c)
        end
    else
        return quote
            Meta.parse($a)
            Meta.parse($b)
        end
    end
end

"""
    gatherops(x)

Gather the operations to the vectors depending on the spatial dimensions.

# Arguments
- `x`: String, the operations.

"""
macro gatherops(x)
    local a,b,c = replace(x, "#"=>1), replace(x, "#"=>2), replace(x, "#"=>3)
    if PeriDyn.SPATIAL_DIMENSIONS_REF[]==3
        return quote
            (
                Meta.parse($a),
                Meta.parse($b),
                Meta.parse($c)
            )
        end
    else
        return quote
            (
                Meta.parse($a),
                Meta.parse($b)
            )
        end
    end
end

"""
    map_reduce(f, op, iter, init)

Map and reduce the operations to the vectors depending on the spatial dimensions.

# Arguments
- `f`: Function, the function.
- `op`: Function, the operation.
- `iter`: Any, the iterator.
- `init`: Any, the initial value.

"""
macro map_reduce(f, op, iter, init)
    if PeriDyn.MULTI_THREAD_REF[]
        quote
            Folds.mapreduce($(esc(f)), $(esc(op)), $(esc(iter)); init=$(esc(init)))
        end
    else
        quote
            mapreduce($(esc(f)), $(esc(op)), $(esc(iter)); init=$(esc(init)))
        end
    end
end

"""
    map_reduce(f, op, iter; init=0.0)

Map and reduce the operations to the vectors depending on the spatial dimensions.

# Arguments
- `f`: Function, the function.
- `op`: Function, the operation.
- `iter`: Any, the iterator.

# Keyword Arguments
- `init`: Any, the initial value.

"""
map_reduce(f, op, iter; init = 0.0) = @map_reduce(f, op, iter, init)

"""
    refresh()

Refresh the functions which depends on macros.

"""
function refresh()
    @eval(
        begin
            PeriDyn.get_magnitude(a) = PeriDyn.@_magnitude(a)
            PeriDyn.get_ij(j, i, a) = PeriDyn.@_ij(j, i, a)
            PeriDyn.map_reduce(f, op, iter; init=0.0) = PeriDyn.@map_reduce(f, op, iter, init)
            PeriDyn.log_detail(x) = PeriDyn.@sdetail(x)
            PeriDyn.log_info(x) = PeriDyn.@sinfo(x)
            PeriDyn.log_impinfo(x) = PeriDyn.@simpinfo(x)
            PeriDyn.log_warning(x) = PeriDyn.@swarning(x)
        end)
end
