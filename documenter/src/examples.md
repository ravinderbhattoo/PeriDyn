# Examples

## Basic examples

### Tensile simulation of a bar
To simulate the tensile behavior of a bar using the PeriDyn package, you can follow these steps:

Activate the PeriDyn environment by running the following code from the PeriDyn package directory. This will activate the environment and install the required packages or install package first.
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```
Import the PeriDyn and PDMaterialPoints packages.
```julia
using PeriDyn
using PDMaterialPoints
```

Create the geometry of the material block using a Cuboid shape and a specified resolution. Then create a material generator using the created geometry and resolution. In this example, we create a bar with dimensions 20x5x5 mm and a resolution of 0.5 mm.

!!! info
    **Units are abstracted in the PeriDyn package, so you can use any unit system as long as it is consistent.**

```julia
resolution = 0.5

x1, v1, y1, vol1, type1 = unpack(create(Cuboid([0 20; 0 5; 0 5]), resolution=resolution))  # mm
mat_gen1 = GeneralMaterial(y1, v1, x1, vol1, type1, 3.0*resolution; max_neigh=200)
```

Define the material parameters such as Young's modulus (Es), Poisson's ratio (nu), density (den), and critical stretch (cstretch). Then create a bond-based material block using the specified material parameters:
```julia
Es = 70 # GPa
nu = 0.15
K = Es/3/(1-2nu)
G = Es/2/(1+nu)
den = 2.2 * 1000.0 # Kg/m3
cstretch = 0.15

mat_spec1 = BondBasedSpecific([K], [cstretch], [den])
# mat_spec1 = OrdinaryStateBasedSpecific([K], [G], [cstretch], [den])

block1 = PeridynamicsMaterial(mat_gen1, mat_spec1)
```

Choose a solver for the simulation. In this example, we use the DSVelocityVerlet solver and set the number of steps (Steps), file write frequency (fwf), and neighbor update frequency (nuf):
```julia
solver = DSVelocityVerlet()
Steps, fwf, nuf = 2000, 100, 10
```

Define the boundary conditions for the simulation. Here, we fix the left part of the bar (y1[1, :] .< 4.0) and move the right part of the bar (y1[1, :] .> 16.0) with a constant velocity:
```julia
BC1 = FixBC(y1[1, :] .< 4.0)
velocity = [0.05*20/Steps/dt, 0.0, 0.0]
BC2 = MoveBC(y1[1, :] .> 16.0, velocity)
```

Define a repulsive contact model. In this example, we use a LinearRepulsionModel with a spring constant (k) and the material block defined earlier (block1):
```julia
k = 1.0
RM1 = LinearRepulsionModel(k, block1; distanceX=3, max_neighs=200)
RMs = [RM1]
```

Create a PeriDyn environment by specifying the material blocks, repulsion models, boundary conditions, and time step (dt):
```julia
dt = 1e-3
env = PeriDyn.Env(1, [block1], RMs, [BC1, BC2], dt)
```

Customize the PeriDyn environment by setting parameters and defining a collection function. In this example, we set the "left" parameter to select the left part of the bar and collect the force values in the "Out" dictionary:
```julia
env.Params = Dict("left" => (y

1[1,:] .< 4))

env.Out = Dict("Force" => zeros(3, Steps))

env.Collect! = function (env, step)
    env.Out["Force"][:, step] = sum(env.f[:, env.Params["left"]], dims=2)
end
```

Run the simulation using the `run!` function and specify the number of steps, solver, file write frequency, output directory, and start index. Additionally, write the simulation output to a JLD file or Ovito files for visualization:
```julia
out_dir = "./output/tensile_sim_BB/$(typeof(solver))"

run!([env], Steps, solver;
    filewrite_freq=fwf, neigh_update_freq=nuf, out_dir=out_dir, start_at=0, ext=:jld)

PeriDyn.write_data("$(out_dir)/env_Out.jld"; Out=env.Out)

PeriDyn.jld2ovito("$(out_dir)/env_1_step_*.jld", Steps; start=0, step=fwf)
```

By following these steps, you can perform a tensile simulation of a bar using the PeriDyn package and visualize the results using Ovito.