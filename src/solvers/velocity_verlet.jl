export velocity_verlet_step!, apply_solver!, DSVelocityVerlet


struct DSVelocityVerlet <: DynamicSolver
end

function apply_solver!(env, solver::DSVelocityVerlet)
    velocity_verlet_step!(env, solver)
end

"""
velocity_verlet_step!(env::GeneralEnv)

Implement a single step of velocity verlet algorithm.
"""
function velocity_verlet_step!(env::GeneralEnv, solver::DSVelocityVerlet)

    c = env.dt/2
    env.v .+= c.*env.f
    env.y .+= env.dt.*env.v

    for bc in env.boundary_conditions
        apply_bc!(env, bc, :position)
    end

    update_acc!(env)
    env.v .+= c.*env.f

    for bc in env.boundary_conditions
        apply_bc!(env, bc, :velocity)
    end

    for bc in env.boundary_conditions
        check!(bc, env)
    end

    env.time_step += 1

    check_boundaries!(env)

end

SOLVERS[:vv] = DSVelocityVerlet

