# Solvers

```@docs
Solver
QuasiStaticSolver
DynamicSolver
apply_solver!
```

## QuasiStaticSolver

```@docs
QSDrag
apply_solver!(env, solver::QSDrag)
minimize!
```

## DynamicSolver

```@docs
DSVelocityVerlet
apply_solver!(env, solver::DSVelocityVerlet)
velocity_verlet_step!
```

## Running the simulation

```@docs
simulate!
run!
```

```@docs
update_acc!
update_neighs!
apply_bc_at0
```
## Saving the simulation state
```@docs
filepath_
save_state!
save_state_ovito_bc!
print_data_file!
```
