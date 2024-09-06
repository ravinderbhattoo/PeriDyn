# Solvers

Here's a breakdown of the solvers implemented in the package:

## Solvers in PeriDyn

The package implements two primary types of solvers ([`Solver`](@ref)): quasi-static solvers and dynamic solvers. 

- **Quasi-static solvers** aim to find the equilibrium state of a system by gradually minimizing forces and energy. 
- **Dynamic solvers**, on the other hand, explicitly compute the time evolution of a system based on forces, velocities, and accelerations. 

The choice of solver depends on the specific phenomenon you're simulating and the desired level of accuracy versus computational expense.

### 1. Dynamic Solvers

#### **[`DSVelocityVerlet`](@ref):**

This solver implements the velocity Verlet algorithm, a numerical integration method commonly used in molecular dynamics and peridynamics simulations. It offers a good balance between accuracy and computational efficiency.

**Implementation Details:**

The [`DSVelocityVerlet`](@ref) solver's implementation follows these steps within each time step:

1.  **Half-Step Velocity Update:**  Calculate the velocity at the half-time step (t + Δt/2) using the acceleration at the current time step (t):  
    `v(t + Δt/2) = v(t) + (Δt/2) * a(t)`
2.  **Position Update:** Update the particle positions using the calculated half-step velocities:  
    `x(t + Δt) = x(t) + Δt * v(t + Δt/2)`
3.  **Apply Boundary Conditions to Positions:** Enforce any defined boundary conditions that constrain particle positions.
4.  **Acceleration Update:** Calculate the acceleration at the new time step (t + Δt) based on the updated positions and the peridynamic force calculations:  
    `a(t + Δt) = F(x(t + Δt)) / m` 
5.  **Velocity Update:**  Calculate the velocity at the new time step (t + Δt) using the acceleration at the new time step (t + Δt):
    `v(t + Δt) = v(t + Δt/2) + (Δt/2) * a(t + Δt)` 
6.  **Apply Boundary Conditions to Velocities:** Enforce any defined boundary conditions that affect particle velocities.
7.  **Check and Update Boundary Conditions:** Perform checks and updates for boundary conditions that change dynamically during the simulation (e.g., `ToFroBC`).
8.  **Time Step Increment:** Advance the simulation time: `t = t + Δt`

### 2. Quasi-static Solvers

#### **[`QSDrag`](@ref)**: 

This solver is designed for quasi-static simulations, where the goal is to reach an equilibrium state rather than explicitly simulating dynamic behaviour. It introduces a drag force to dampen the system's motion, allowing it to settle into a stable configuration.

**Implementation Details:**

The [`QSDrag`](@ref) solver uses a drag force to gradually reduce kinetic energy and drive the system towards equilibrium. 

1.  **Apply Boundary Conditions:** Initially, it applies all specified boundary conditions to both position and velocity.
2.  **Iterative Minimization:** The core of the `QSDrag` solver lies in its iterative minimization process, which continues until:
    - The maximum number of iterations (specified by `max_iter`) is reached, or
    - Both the maximum particle displacement and the maximum force fall below their respective tolerances (`x_tol` and `f_tol`).
3.  **Drag Force Calculation:** In each iteration, the solver calculates a drag force for each particle based on its velocity.  The drag force opposes the direction of motion, effectively slowing down the particles.
4.  **Position and Velocity Update:** Using the calculated drag forces, the solver updates particle positions and velocities. Importantly, it limits both displacement and velocity within each step to maintain stability during the minimization process.
5.  **Boundary Condition Checks and Updates:**  After updating the system's state, the solver checks for boundary condition violations and updates the system accordingly. This step is crucial to ensure that the system converges to a state that satisfies all imposed constraints.
6.  **Time Step Increment:** Once the iteration criteria are met, the solver increments the time step.

Implementations related to solvers:

- [`Solver`](@ref)
- [`QuasiStaticSolver`](@ref)
    - [`QSDrag`](@ref)
    - [`apply_solver!(::Any, ::QSDrag)`](@ref)
    - [`minimize!`](@ref)
- [`DynamicSolver`](@ref)
    - [`DSVelocityVerlet`](@ref)
    - [`apply_solver!(::Any, ::DSVelocityVerlet)`](@ref)
    - [`velocity_verlet_step!`](@ref)


## Running the Simulation

The simulation can be run using the [`simulate!`](@ref) or [`run!`](@ref) functions.

- **[`simulate!`](@ref)** is a higher-level function that takes care of setting up the output directory and calling the [`run!`](@ref) function.
- **[`run!`](@ref)** is the core function that executes the simulation loop.

### The [`simulate!`](@ref) Function

The [`simulate!`](@ref) function has the following signature:

```julia
function simulate!(args...; out_dir="datafile", append_date=true, kwargs...)
```

**Arguments:**

- `args...`: Arguments passed to the [`run!`](@ref) function.
- `out_dir="datafile"`: Directory where the data files are saved.
- `append_date=true`: If `true`, the date is appended to the `out_dir` path.
- `kwargs...`: Keyword arguments passed to the [`run!`](@ref) function.

**Usage:**

```julia
simulate!(envs, steps, solver; kwargs...)
```

where:

- `envs` is an array of `Environment` ([`GeneralEnvironment`](@ref)) objects representing the simulation environments.
- `steps` is the number of simulation steps to run.
- `solver` is a [`Solver`](@ref) object specifying the numerical integration scheme to use.

### The [`run!`](@ref) Function

The [`run!`](@ref) function has the following signature:

```julia
function run!(envs, N::Int64, solver::Solver; filewrite_freq::Int64=10,
    neigh_update_freq::Int64=1, average_prop_freq::Int64=1,
    out_dir::String="datafile",
    start_at::Int64=0, write_from::Int=0, ext::Symbol=:jld,
    max_part=30)
```

**Arguments:**

- `envs`: Array of `Environment` ([`GeneralEnvironment`](@ref)) objects representing the simulation environments.
- `N`: Number of simulation steps to run.
- `solver`: [`Solver`](@ref) object specifying the numerical integration scheme to use.

**Keyword Arguments:**

- `filewrite_freq::Int64=10`: Frequency of writing data files to disk.
- `neigh_update_freq::Int64=1`: Frequency of updating neighbour lists.
- `average_prop_freq::Int64=1`: Frequency of calculating average properties.
- `out_dir::String="datafile"`: Directory where the data files are saved.
- `start_at::Int64=0`: Starting step.
- `write_from::Int=0`: Starting index of the data files.
- `ext::Symbol=:jld`: Extension of the data files.
- `max_part=30`: kwargs for neighbour update function.

**Simulation Loop:**

The [`run!`](@ref) function performs the following steps in a loop:

1. **Update Neighbour Lists:** Update the neighbour lists for each `Environment` ([`GeneralEnvironment`](@ref)) if `i % neigh_update_freq == 0`, where `i` is the current simulation step.
2. **Apply Solver:** Apply the specified `solver` to each `Environment` ([`GeneralEnvironment`](@ref)).
3. **Print Data:** Write the simulation data to disk if `i % filewrite_freq == 0`.

This loop continues until the specified number of `steps` have been completed.


## Updating the Simulation State

The simulation state is updated by the following functions:

- **[`update_acc!](@ref)`**: This function updates the acceleration of all material points in a simulation environment. It calculates the acceleration due to material deformation ([]`update_mat_acc!`](#ref)) and contact forces ([`update_contact_acc!`](#ref)) and updates the momentum ([`update_misc!`](#ref)).
- **[`update_neighs!`](#ref)**: Updates the neighbour lists of material points for contact force calculations. This is important for efficient computation of contact forces between particles that are close to each other.
- **[`apply_bc_at0`](#ref)**: This function applies the boundary conditions at the beginning of the simulation (`start_at` = 0). This sets the initial conditions for the simulation, such as fixing certain particles in place or assigning initial velocities.

Within the main simulation loop, these functions are called at specific frequencies determined by the `neigh_update_freq` and `filewrite_freq` parameters. For example, the neighbour lists are updated every `neigh_update_freq` steps. The acceleration is updated at each step of the simulation.

