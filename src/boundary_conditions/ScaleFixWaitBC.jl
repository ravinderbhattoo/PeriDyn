export ScaleFixWaitBC, MoveBC

struct ScaleFixWaitBC <: BoundaryCondition
    @general_bc_p
    checkF::Function
end


# Definitions
function ScaleFixWaitBC(bool, scale, fixpoint, wait, scalebool; applyafter=0, onlyatstart=false)
    last = zeros(Float64, SPATIAL_DIMENSIONS_REF[], sum(bool))
    xF = (env, BC) -> (BC.last, BC.last)
    vF = (env, BC) -> (0.0, BC.last)
    checkF = (env, BC) -> begin
        if (env.time_step % wait == 0)
            y = env.y[:, scalebool]
            y = (y .- vec(fixpoint)) .* sqrt.(1 .+ 2*vec(scale)*env.dt) .+ vec(fixpoint)
            env.y[:, scalebool] .= y
            BC.last .= env.y[:, BC.bool]
        end   
    end
    ScaleFixWaitBC(bool, last, onlyatstart, xF, vF, checkF)
end

# Implementation
function apply_bc_at0!(env, BC::ScaleFixWaitBC)
    BC.last .= env.y[:, BC.bool]
end

# Check
function check!(BC::ScaleFixWaitBC, env)
    BC.checkF(env, BC)
end

