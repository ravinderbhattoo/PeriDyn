using Pkg
Pkg.activate(".")

using Revise
using PeriDyn
using PDMesh
using CUDA

device = :cuda

PeriDyn.set_device(device)
PDMesh.set_device(device)

resolution = 2

x1, v1, y1, vol1, type1 = unpack(create(Cuboid([-5 5; -5 5; -5 5]); resolution=resolution, type=1))

hor1 = 3.0*resolution
mat_gen1 = GeneralMaterial(y1, v1, x1, vol1, type1, hor1, max_neigh=200)

K, G = 100.0, 1.0
den1 = 1000.0
cstretch1 = 0.5

# mat_spec1 = BondBasedSpecific([K], [cstretch1], [den1])
mat_spec1 = SkipSpecific([den1])
# mat_spec1 = OrdinaryStateBasedSpecific([K], [G], [cstretch1], [den1])

block1 = PeridynamicsMaterial(mat_gen1, mat_spec1; name="block 1")

BLOCKS = [block1]

x2, v2, y2, vol2, type2 = unpack(create(Cuboid([-5 5; -5 5; 15 25]); resolution=resolution, type=2))

# v2 = Array(v2)
# v2[3,:] .-= 0.01


hor2 = 3.0*resolution
mat_gen2 = GeneralMaterial(y2, v2, x2, vol2, type2, hor2, max_neigh=200)



K2, G2 = 10.0, 10.0
den2 = 1000.0
cstretch2 = 0.5

# mat_spec2 = BondBasedSpecific([K2], [cstretch2], [den2])
mat_spec2 = SkipSpecific([den2])
# mat_spec2 = OrdinaryStateBasedSpecific([K2], [G2], [cstretch2], [den2])

type_ = 2
block2 = PeridynamicsMaterial(mat_gen2, mat_spec2; name="block $type_")

k = 1.0
max_neighs = 400

RM = LinearRepulsionModel(k, block1, block2, distanceX=2, max_neighs=max_neighs)
RM1 = LinearRepulsionModel(k/1000, block1, distanceX=2, max_neighs=max_neighs)
RM2 = LinearRepulsionModel(k/1000, block2, distanceX=2, max_neighs=max_neighs)

env = Env(1, [block1, block2], [RM], [], 0.01)

PeriDyn.@check_nan env.v "v"


solver = DSVelocityVerlet()
Steps, fwf, nuf = 20000, 100, 10

out_dir="./output/$(length(BLOCKS))blocks"


run!([env], Steps, solver;
    filewrite_freq=fwf, neigh_update_freq=100, out_dir=out_dir, start_at=0, ext=:data)

# PeriDyn.jld2ovito("$(out_dir)/env_1_step_*.jld", steps1; start=0, step=20)
