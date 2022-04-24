export FixBC

struct FixBC <: BoundaryCondition
    @general_bc_p
    # bool::Array{Bool, 1}
    # last::Array{Float64, 2}
    # onlyatstart::Bool
    # xF::Function
    # vF::Function
end


"""
    FixBC(bool; onlyatstart=false)

"""
function FixBC(bool; onlyatstart=false)
    last = zeros(Float64, SPATIAL_DIMENSIONS_REF[], sum(bool))
    xF = (env, BC) -> (BC.last, BC.last)
    vF = (env, BC) -> (0.0, BC.last)
    FixBC(bool, last, onlyatstart, xF, vF)
end

function apply_bc_at0!(env, BC::FixBC)
    BC.last .= env.y[:, BC.bool]
end

