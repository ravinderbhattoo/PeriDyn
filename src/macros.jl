export TIMEIT_REF, CHECK_NAN, SPATIAL_DIMENSIONS_REF, MULTI_THREAD_REF
export check_nan, timeit, sinfo, swarning, sdetail, simpinfo, printsize, whatis
export _magnitude, _ij
export refresh, set_multi_threading

macro def(name, definition)
    return quote
        macro $(esc(name))()
            esc($(Expr(:quote, definition)))
        end
    end
end

TIMEIT_REF = Ref(false)
CHECK_NAN = Ref(true)
MULTI_THREAD_REF = Ref(true)
SPATIAL_DIMENSIONS_REF = Ref(3)
LOGLEVEL=Ref(2)

macro swarning(x)
    if LOGLEVEL[] >= 1
        quote
            printstyled("WARNING: ", $(esc(x)), "\n"; color = :red)
        end
    else
    end
end

log_warning(x) = @swarning(x)

macro simpinfo(x)
    if LOGLEVEL[] >= 2
        quote
            printstyled("INFO: ", $(esc(x)), "\n"; color = :yellow, bold = true)
        end
    else
    end
end

log_impinfo(x) = @simpinfo(x)

macro sinfo(x)
    if LOGLEVEL[] >= 3
        quote
            printstyled("INFO: ", $(esc(x)), "\n"; color = :green)
        end
    else
    end
end

log_info(x) = @sinfo(x)

macro sdetail(x)
    if LOGLEVEL[] >= 4
        quote
            printstyled("INFO: ", $(esc(x)), "\n"; color = :blue)
        end
    else
    end
end

log_detail(x) = @sdetail(x)

function set_loglevel(x::Int64)
    LOGLEVEL[] = x
    @simpinfo "Log level: $(LOGLEVEL[])"
    refresh()
end


function set_multi_threading(x)
    MULTI_THREAD_REF[] = x
    refresh()
    println("Multi threading: ", x)
end

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

function printsize(var...)
    for var1 in var
        print(size(var1), ", ")
    end
    var
end

macro timeit(ex, name)
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

get_magnitude(a) = @_magnitude(a)

macro _ij(j, i, x)
    if PeriDyn.SPATIAL_DIMENSIONS_REF[]==3
        quote
            [$(esc(x))[1,$(esc(j))] - $(esc(x))[1,$(esc(i))],
                    $(esc(x))[2,$(esc(j))] - $(esc(x))[2,$(esc(i))],
                    $(esc(x))[3,$(esc(j))] - $(esc(x))[3,$(esc(i))]]
        end
    else
        quote
            [$(esc(x))[1,$(esc(j))] - $(esc(x))[1,$(esc(i))],
                    $(esc(x))[2,$(esc(j))] - $(esc(x))[2,$(esc(i))]]
        end
    end
end

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



get_ij(j, i, a) = @_ij(j, i, a)

function whatis(a, b, c)
    println(a)
    println(b)
    println(c)
end

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

map_reduce(f, op, iter; init = 0.0) = @map_reduce(f, op, iter, init)

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
