```@meta
CurrentModule = PeriDyn
```

# PeriDyn.jl Documentation

## Simulation Environment
```@docs
    Env(id::Int64,materials,short_range_repulsion,boundary_conds,dt;state=2)
```
## Mesh
```@docs
    create_block(lattice::Array{Float64,1}, N::Array{Int64,1})
```

## Boundary Conditions
```@docs
    apply_bc!(env,BC::FixBC)
    apply_bc!(env,BC::MoveBC)
    apply_bc!(env,BC::ScaleBC)
    apply_bc!(env,BC::JustScaleBC)
    apply_bc!(env,BC::ScaleFixBC)
```

## Solvers
```@docs
velocity_verlet!(envs::Any,N::Int64;freq1=10,freq2=50,file_prefix="datafile",start_at::Int64=0)
minimize!(env::GeneralEnv,step_size::Float64; max_iter::Int64=500,min_step_tol_per::Float64=5.0)
quasi_static!(envs::Any,N::Int64,step_size::Float64; max_iter::Int64=100, min_step_tol_per::Float64=0.5, freq1::Int64=10, freq2::Int64=50, file_prefix::String="datafile",start_at::Int64=0)
```

## Util
```@docs
    write_data(filename,every,type,pos)
    write_data(filename::String,type::Array{Int64,1},y::Array{Float64,2})
    write_data(filename::String,type::Array{Int64,1},y::Array{Float64,2}, v::Array{Float64,2},f::Array{Float64,2},)
```
