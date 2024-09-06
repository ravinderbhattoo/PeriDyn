# Automatic documentation

__Table of contents on this page.__
```@contents
Pages = [
            "autodocs.md"
        ]
Depth = 4
```

## Peridynamics functions

```@autodocs
Modules = [PeriDyn]
Pages   = ["peri.jl"]
```

## Simulation environment

```@autodocs
Modules = [PeriDyn]
Pages   = ["environment.jl"]
```

## Material models

### General material models and functions
```@autodocs
Modules = [PeriDyn]
Pages = ["material.jl", "general_material.jl"]
```

### Specific material models and functions
```@autodocs
Modules = [PeriDyn]
Pages = ["skip_specific.jl", "bond_based.jl", "ordinary_state_based.jl", "elasto_plastic.jl"]
```

## Contact models and functions

```@autodocs
Modules = [PeriDyn]
Pages   = ["contacts.jl", "nonlinearspring.jl", "linearspring.jl", "shortrange.jl"]
```

## Boundary conditions

```@autodocs
Modules = [PeriDyn]
Pages   = ["boundary_conditions.jl", "FixBC.jl", "ToFroBC.jl", "DeltaScaleBC.jl", "ScaleFixWaitBC.jl", "ContainerBC.jl"]
```

## Solvers and integrators

### General solver functions
```@autodocs
Modules = [PeriDyn]
Pages   = ["solvers.jl"]
```

### Quasi-static solver
```@autodocs
Modules = [PeriDyn]
Pages   = ["minimize.jl"]
```

### Dynamic solver
```@autodocs
Modules = [PeriDyn]
Pages   = ["velocity_verlet.jl"]
```

## Input and output functions

```@autodocs
Modules = [PeriDyn]
Pages   = ["io.jl", "ff_ovito.jl"]
```


## Miscellaneous functions

### Logging
```@autodocs
Modules = [PeriDyn]
Pages   = ["log.jl"]
```

### Macros
```@autodocs
Modules = [PeriDyn]
Pages   = ["const_macros.jl"]
```

### Representation
```@autodocs
Modules = [PeriDyn]
Pages   = ["repr.jl"]
```

### Utilities
```@autodocs
Modules = [PeriDyn]
Pages   = ["utils.jl"]
```