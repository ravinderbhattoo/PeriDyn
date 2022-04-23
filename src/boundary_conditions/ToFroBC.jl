export ToFroBC, MoveBC

struct ToFroBC <: BoundaryCondition
    @general_bc_p
    direction::Ref{Int64}
    freq::Union{Int64,Float64}
    applyafter::Int64  # apply freq after this many steps
end


# Definitions
function ToFroBC(bool, rate, freq; applyafter=0, onlyatstart=false)
    last = zeros(Float64, SPATIAL_DIMENSIONS_REF[], sum(bool))
    xF = (env, BC) -> begin
        y = BC.last .+ vec(env.dt * BC.direction[] * rate)
        (y, y)
    end
    vF = (env, BC) -> ( vec(BC.direction[] * rate), BC.last)
    ToFroBC(bool, last, onlyatstart, xF, vF, Ref(1), freq, applyafter)
end

function MoveBC(bool, rate; kwargs...)
    ToFroBC(bool, rate, Inf; kwargs...)
end



# Implementation
function apply_bc_at0!(env, BC::ToFroBC)
    BC.last .= env.y[:, BC.bool]
end

# Check
function check!(BC::ToFroBC, env)
    if ( env.time_step > BC.applyafter ) & (env.time_step % BC.freq == 0)
        BC.direction[] = -1*BC.direction[]
    end    
end

