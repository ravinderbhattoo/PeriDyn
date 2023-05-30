export velocity_verlet_step!, apply_solver!, DSVelocityVerlet


"""
    DSVelocityVerlet

The velocity verlet algorithm.
"""
struct DSVelocityVerlet <: DynamicSolver
end

"""
    apply_solver!(env, solver::DSVelocityVerlet)

Apply the velocity verlet algorithm to the environment.

# Arguments
- `env`: GeneralEnv, the environment.
- `solver`: DSVelocityVerlet, the velocity verlet algorithm.

# Example
```
solver = DSVelocityVerlet()
apply_solver!(env, solver)
```
"""
function apply_solver!(env, solver::DSVelocityVerlet)
    velocity_verlet_step!(env, solver)
end

"""
    velocity_verlet_step!(env, solver::DSVelocityVerlet)

Apply one step of the velocity verlet algorithm to the environment.

# Arguments
- `env`: GeneralEnv, the environment.
- `solver`: DSVelocityVerlet, the velocity verlet algorithm.

# Example
```
solver = DSVelocityVerlet()
velocity_verlet_step!(env, solver)
```
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
        check!(env, bc)
    end

    env.time_step += 1

    check_boundaries!(env)

end

SOLVERS[:vv] = DSVelocityVerlet

