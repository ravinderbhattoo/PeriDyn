export velocity_verlet_step!, apply_solver!, DSVelocityVerlet


"""
    DSVelocityVerlet

The velocity verlet algorithm.
"""
struct DSVelocityVerlet <: DynamicSolver
end

"""
    apply_solver!(env, solver::DSVelocityVerlet)

Apply the velocity verlet algorithm to the simulation environment.

# Arguments
- `env`: the simulation environment.
- `solver`: DSVelocityVerlet, the velocity verlet algorithm.

# Example
```
solver = DSVelocityVerlet()
apply_solver!(env, solver)
```
"""
function apply_solver!(env, solver::DSVelocityVerlet)
    @timeit timings "VelocityVerlet step" velocity_verlet_step!(env, solver)
end

"""
    velocity_verlet_step!(env, solver::DSVelocityVerlet)

Apply one step of the velocity verlet algorithm to the simulation environment.

### Steps:

1. Update the velocity.
2. Update the position.
3. Apply the boundary conditions to the position.
4. Update the acceleration.
5. Update the velocity.
6. Apply the boundary conditions to the velocity.
7. Check the boundary conditions.
8. Check the boundaries.

# Arguments
- `env`: the simulation environment.
- `solver`: DSVelocityVerlet, the velocity verlet algorithm.

# Example
```
solver = DSVelocityVerlet()
velocity_verlet_step!(env, solver)
```
"""
function velocity_verlet_step!(env, solver::DSVelocityVerlet)

    c = env.dt/2
    dt = env.dt

    @timeit timings "state update 1" begin
        @. env.v += $c*env.f
        @. env.y += $dt*env.v
    end

    for bc in env.boundary_conditions
        apply_bc!(env, bc, :position)
    end

    update_acc!(env)

    @timeit timings "state update 2" begin
        @. env.v += $c.*env.f
    end

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

