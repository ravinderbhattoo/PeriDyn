export TIMEIT_REF, CHECK_NAN, SPATIAL_DIMENSIONS_REF
export check_nan, timeit
export _magnitude, _ij
export refresh

TIMEIT_REF = Ref(false)
CHECK_NAN = Ref(true)
SPATIAL_DIMENSIONS_REF = Ref(3)

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
            SVector($(esc(x))[1,$(esc(j))] - $(esc(x))[1,$(esc(i))], 
                    $(esc(x))[2,$(esc(j))] - $(esc(x))[2,$(esc(i))], 
                    $(esc(x))[3,$(esc(j))] - $(esc(x))[3,$(esc(i))])
        end
    else
        quote
            SVector($(esc(x))[1,$(esc(j))] - $(esc(x))[1,$(esc(i))], 
                    $(esc(x))[2,$(esc(j))] - $(esc(x))[2,$(esc(i))])
        end
    end
end

get_ij(j, i, a) = @_ij(j, i, a) 

function refresh()
    @eval(
        begin
            PeriDyn.get_magnitude(a) = PeriDyn.@_magnitude(a)
            PeriDyn.get_ij(j, i, a) = PeriDyn.@_ij(j, i, a)
        end)
end
