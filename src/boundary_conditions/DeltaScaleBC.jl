export DeltaScaleBC

struct DeltaScaleBC <: BoundaryCondition
    @general_bc_p
end

# Definitions
function DeltaScaleBC(bool, scale, fixpoint; onlyatstart=false)
    scale = deviceconvert(scale)
    fixpoint = deviceconvert(fixpoint)
    last = zeros(Float64, SPATIAL_DIMENSIONS_REF[], sum(bool))
    xF = (env, BC) -> begin
                y = (BC.last .- reshape(fixpoint, :)) .* (1 .+ reshape(scale, :)*env.dt*env.time_step) .+ reshape(fixpoint, :)
                (y, BC.last)
            end
    vF = (env, BC) -> ( 0.0, BC.last )
    DeltaScaleBC(bool, last, onlyatstart, xF, vF)
end


# Implementation
function apply_bc_at0!(env, BC::DeltaScaleBC)
    BC.last .= env.y[:, BC.bool]
end

# Check
function check!(BC::DeltaScaleBC, env)
    # pass
end