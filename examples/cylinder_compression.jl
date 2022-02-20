using Pkg

Pkg.activate(".")

using PeriDyn
using PDMesh

using Random
Random.seed!(42)

resolution = 1.0
x1, v1, y1, vol1, type1 = unpack(create(Cylinder(20.0, 5.0, 100.0), resolution=resolution, type=1, rand_=0.02))

hor1 = 3.0*resolution
mat_gen1 = GeneralMaterial(y1, v1, x1, vol1, type1, hor1, max_neigh=200)

K, G = 1.0, 1.0
den1 = 1000.0
cstretch1 = 0.2
mat_spec1 = OrdinaryStateBasedSpecific([K], [G], [cstretch1], [den1])
# mat_spec1 = SkipSpecific([den1])

block1 = PeridynamicsMaterial(mat_gen1, mat_spec1)

k = 10.0
RM1 = LinearRepulsionModel(k, block1, distanceX=2, max_neighs=200)

bc1 = FixBC(y1[3, :] .< 5.0)
bc2 = MoveBC(y1[3, :] .> 95.0, [0.0, 0.0, -0.01])

env = Env(1, [block1], [RM1], [bc1, bc2], 0.5)

velocity_verlet!([env], 2000, filewrite_freq=100, neigh_update_freq=10, out_dir="./output/cylinder", start_at=0)

#
