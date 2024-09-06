# define some types
QF = Union{Float64, Unitful.Quantity{Float64, D, U} where {D, U}}

function Base.convert(::Type{Union{Float64, Quantity{Float64}}}, x::Union{Int, Quantity{Int}})
    Float64(x)
end

function set_default_units(;ulength = u"mm", umass = u"mg", utime = u"ms")

    uarea = ulength*ulength
    uvolume = uarea*ulength

    uweighted_vol = ulength * uvolume
    udensity = umass/uvolume

    uforce = umass*ulength/utime/utime
    uvelocity = ulength/utime
    uacceleration = uvelocity/utime
    umomentum = umass*uvelocity

    uforce_density = uforce/uvolume/uvolume # tij
    ubulk_modulus = uforce/uarea

    uplastic_yield_value = ubulk_modulus^2/ulength^5

    PeriDyn.DEFAULT_UNITS = Dict(
        Unitful.dimension(i)=>i for i in
        [u"m/m", ulength, umass, utime,
        uarea, uvolume, uweighted_vol, udensity,
        uvelocity, uacceleration, umomentum,
        uforce, uforce_density,
        ubulk_modulus,
        uplastic_yield_value
        ]
    )
end


DEFAULT_UNITS = Dict()

# set SI_UNITS
SI_UNITS = set_default_units(;
                ulength = u"m",
                umass = u"kg",
                utime = u"s")

# set DEFAULT_UNITS
set_default_units(;ulength=u"mm", umass=u"mg", utime=u"μs")

############################################
# uconvert
############################################

function change_unit!(DIMS, x)
    DIMS[dimension(x)] = x
end

function uconvert_to(DIMS, x)
    s = DIMS[Unitful.dimension(x)]
    return Unitful.uconvert(s, x)
end

function uconvert_to(DIMS, x::AbstractArray)
    u = unit(eltype(x))
    s = DIMS[Unitful.dimension(eltype(x))]
    return Unitful.uconvert(s, 1u) * ustrip(x)
end

function uconvert_to(DIMS, x::Dict)
    return Dict(k => uconvert_to(DIMS, v) for (k,v) in x)
end

function uconvert_to(DIMS, x::T) where T<:Union{Function,
    Base.RefValue,
    String,
    Pair{AbstractVector{Int64}, AbstractVector{Int64}},
    BitVector,
    UnitRange{Int64},
    Nothing
    }
    return x
end

function uconvert_to_default(x)
    uconvert_to(DEFAULT_UNITS, x)
end

function uconvert_to_si(x)
    uconvert_to(SI_UNITS, x)
end

############################################
# ustrip
############################################
function ustrip(x)
    return 1Unitful.ustrip(x)
end

function ustrip(x::Dict)
    return Dict(k => ustrip(v) for (k,v) in x)
end

function ustrip(x::T) where T<:Union{Function,
    Base.RefValue,
    String,
    Pair{AbstractVector{Int64}, AbstractVector{Int64}},
    BitVector,
    UnitRange{Int64},
    Nothing
    }
    return x
end

function ustrip_to_default(x)
    ustrip(uconvert_to(DEFAULT_UNITS, x))
end

function ustrip_to_si(x)
    ustrip(uconvert_to(SI_UNITS, x))
end





