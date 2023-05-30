# Boundary conditions

Boundary conditions in PeriDyn are used to specify the behavior of material particles at the selected material location in a simulation. The page lists different types of boundary conditions that are predefined in the package.

## Predefined boundary conditions

- [`FixBC`](bc.html/#PeriDyn.FixBC) (Fix Boundary Condition): This boundary condition fixes the position and velocity of particles at specific boundary regions. It is typically used to simulate rigid walls or immovable boundaries. For example, it can be used to fix the position of particles at the start of the simulation to represent a fixed boundary.

```@docs
FixBC
FixBC(bool; onlyatstart=false)
```

- [`ToFroBC`](bc.html/#PeriDyn.ToFroBC) (To and Fro Boundary Condition): This boundary condition imparts a prescribed motion to particles at specific boundary regions. It allows particles to move back and forth within a specified direction and frequency. It is useful for simulating boundaries with oscillatory motion. For instance, it can be used to simulate a boundary that moves back and forth periodically.

```@docs
ToFroBC
ToFroBC(bool, rate, freq; applyafter=0, onlyatstart=false)
```

- [`MoveBC`](bc.html/#PeriDyn.MoveBC) (Move Boundary Condition): This boundary condition imparts a constant velocity to particles at specific boundary regions. It is commonly used to simulate boundaries subjected to external forces or displacements. For example, it can be applied to particles at the boundary of a bar to represent a constant velocity movement.

```@docs
MoveBC
```

- [`DeltaScaleBC`](bc.html/#PeriDyn.DeltaScaleBC) (Delta Scale Boundary Condition): This boundary condition applies a scaling factor to particles at specific boundary regions. It can be used to simulate deformations or changes in the size of the material. For instance, it can be used to scale particles near a boundary to represent a stretching or compression effect.

```@docs
DeltaScaleBC
DeltaScaleBC(bool, scale, fixpoint; onlyatstart=false)
```

- [`ScaleFixWaitBC`](bc.html/#PeriDyn.ScaleFixWaitBC) (Scale Fix Wait Boundary Condition): This boundary condition applies a scaling factor to particles after a specified number of time steps. It allows particles to be scaled or deformed over time. It can be useful for simulating time-dependent deformations. For example, it can be used to gradually apply a scaling effect to particles at the boundary after a certain waiting period.

```@docs
ScaleFixWaitBC
ScaleFixWaitBC(bool, scale, fixpoint, wait, scalebool; applyafter=0, onlyatstart=false)
```

# Custom boundary conditions
A custom boundary condition can also be defined by the user which should be a subtype of [`BoundaryCondition`](bc.html/#PeriDyn.BoundaryCondition) abstract type. These boundary conditions are applied to the position and velocity aspects of the particles in the simulation defined by `xF` and `vF` functions respectively. They define the behavior of the particles at the boundaries and can be customized based on the specific requirements of the simulation. [`apply_bc!`](bc.html/#PeriDyn.apply_bc!) and [`apply_bc_at0!`](bc.html/#PeriDyn.apply_bc_at0!) functions are used to apply the boundary conditions to the particles in the simulation. The [`check!`](bc.html/#PeriDyn.check!) function is used to check if the boundary conditions are applied correctly to the particles in the simulation. [`apply_bc!`](bc.html/#PeriDyn.apply_bc!), [`apply_bc_at0!`](bc.html/#PeriDyn.apply_bc_at0!) and [`check!`](bc.html/#PeriDyn.check!) functions are defined for [`BoundaryCondition`](bc.html/#PeriDyn.BoundaryCondition) abstract type and can be overloaded by the user to define custom boundary conditions. The default implementation of these functions is as follows.

```julia
"""
    apply_bc!(env, BC::BoundaryCondition, ::Type{Val{:position}})

Apply the general boundary condition `BC` to the position aspect
of the given environment `env`.

# Arguments
- `env`: The environment to which the boundary condition is applied.
- `BC`: The general boundary condition to apply.
"""
function apply_bc!(env, BC::T, ::Type{Val{:position}}) where
                                        T <: BoundaryCondition
    a, b = BC.xF(env, BC)
    env.y[:, BC.bool] .= a
    BC.last .= b
end

"""
    apply_bc!(env, BC::BoundaryCondition, ::Type{Val{:velocity}})

Apply the general boundary condition `BC` to the velocity aspect of the
given environment `env`.

# Arguments
- `env`: The environment to which the boundary condition is applied.
- `BC`: The general boundary condition to apply.
"""
function apply_bc!(env, BC::T, ::Type{Val{:velocity}}) where
                                        T <: BoundaryCondition
    a, b = BC.vF(env, BC)
    env.v[:, BC.bool] .= a
    BC.last .= b
end

"""
    apply_bc_at0!(env, BC::BoundaryCondition)

Apply the boundary condition `BC` to the given environment `env` at the
start of the simulation (t=0).
"""
function apply_bc_at0!(env, BC::T) where T <: BoundaryCondition
    # error("check! method Not implemented for $(T)")
    log_detail("apply_bc_at0! method not implemented for $(T)")
end

"""
    check!(env, BC::BoundaryCondition)

Check if the boundary condition `BC` has changed and updates the `BC`.
Used for dynamic boundary conditions.
"""
function check!(env, BC::T) where T <: BoundaryCondition
    # error("check! method Not implemented for $(T)")
    log_detail("check! method not implemented for $(T)")
end

```

The following code shows how one can define the `CustomNameBC` boundary condition.

```julia
"""
This module contains definitions for CustomNameBC boundary conditions.
"""

export CustomNameBC

"""
    CustomNameBC

Struct representing the CustomNameBC boundary condition.

# Fields
- `bool`: Boolean array specifying the affected elements.
- `last`: Last position of the affected elements.
- `onlyatstart`: Flag indicating if the boundary condition is applied
    only at the start.
- `xF`: function for updating the position.
- `vF`: function for updating the velocity.
"""
struct CustomNameBC <: BoundaryCondition
    @general_bc_p
    # macro will expand to:
    # bool::AbstractArray{Bool, 1}
    # last::AbstractArray{Float64, 2}
    # onlyatstart::Bool
    # xF::Function
    # vF::Function
end


"""
    CustomNameBC(args...; onlyatstart=false)

Construct a CustomNameBC boundary condition.

# Arguments
- `args`: Arguments for the boundary condition.
- `onlyatstart`: Flag indicating if the boundary condition is applied only
    at the start (default: `false`).

# Returns
A CustomNameBC object representing the boundary condition.
"""
function CustomNameBC(args...; onlyatstart=false)
    last = zeros(Float64, SPATIAL_DIMENSIONS_REF[], sum(bool))
    bool = # generate bool array depending on the args passed
    xF = (env, BC) -> begin
                        # write logic here
                        (a, b)  # return a and b as the new position
                                # and BC.last position
                    end
    vF = (env, BC) -> begin
                        # write logic here
                        (a, b)  # return a and b as the new velocity
                                # and BC.last position
    deviceconvert(CustomNameBC(bool, last, onlyatstart, xF, vF))
end

"""
    apply_bc_at0!(env, BC::CustomNameBC)

Apply the CustomNameBC boundary condition at time 0 to the given
environment `env`.

# Arguments
- `env`: The environment to which the boundary condition is applied.
- `BC`: The CustomNameBC boundary condition to apply.
"""
function apply_bc_at0!(env, BC::CustomNameBC)
    # write logic here for e.g. BC.last .= env.y[:, BC.bool]
end

"""
    check!(env, BC::CustomNameBC)

Check if the CustomNameBC boundary condition has changed and updates
the `CustomNameBC`. Used for dynamic boundary conditions.

# Arguments
- `env`: The environment to which the boundary condition is applied.
- `BC`: The CustomNameBC boundary condition to apply.
"""
function check!(env, BC::CustomNameBC)
    # write logic here for e.g. BC.bool .= env.y[:, 1] .> 0.5
end
```

```@docs
BoundaryCondition
```

```@docs
apply_bc!
apply_bc_at0!
check!
```
