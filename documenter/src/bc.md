# Boundary conditions

```@docs
BoundaryCondition
```

```@docs
apply_bc!
apply_bc_at0!
check!
```

```@docs
FixBC
FixBC(bool; onlyatstart=false)
```

```@docs
ToFroBC
ToFroBC(bool, rate, freq; applyafter=0, onlyatstart=false)
MoveBC
```

```@docs
DeltaScaleBC
DeltaScaleBC(bool, scale, fixpoint; onlyatstart=false)
```

```@docs
ScaleFixWaitBC
ScaleFixWaitBC(bool, scale, fixpoint, wait, scalebool; applyafter=0, onlyatstart=false)
```
