using Pkg

Pkg.activate(".")

using Revise

using PeriDyn
using PDMesh

using Random


Random.seed!(42)

resolution = 2

x1, v1, y1, vol1, type1 = unpack(create(Cuboid([-5 5; -5 5; 0 10]), resolution=resolution, type=1))

hor1 = 3.0*resolution
mat_gen1 = GeneralMaterial(y1, v1, x1, vol1, type1, hor1, max_neigh=200)

K, G = 1.0, 1.0
den1 = 1000.0
cstretch1 = 0.5
mat_spec1 = OrdinaryStateBasedSpecific([K], [G], [cstretch1], [den1])
# mat_spec1 = SkipSpecific([den1])

block1 = PeridynamicsMaterial(mat_gen1, mat_spec1; name="block 1")

BLOCKS = [block1]
RMS = []
    
x2, v2, y2, vol2, type2 = unpack(create(Cuboid([-5 5; -5 5; 15 25]), resolution=resolution, type=2))

v2[3,:] .-= 0.01

hor2 = 3.0*resolution
mat_gen2 = GeneralMaterial(y2, v2, x2, vol2, type2, hor2, max_neigh=200)

K2, G2 = 10.0, 10.0
den2 = 1000.0
cstretch2 = 0.5
mat_spec2 = OrdinaryStateBasedSpecific([K2], [G2], [cstretch2], [den2])
# mat_spec2 = SkipSpecific([den2])

type_ = 2
block2 = PeridynamicsMaterial(mat_gen2, mat_spec2; name="block $type_")

push!(BLOCKS, block2)

k = 1000.0
RM = LinearRepulsionModel(k, block1, block2, distanceX=2, max_neighs=200)
RM1 = LinearRepulsionModel(k/1000, block1, distanceX=2, max_neighs=200)
RM2 = LinearRepulsionModel(k/1000, block2, distanceX=2, max_neighs=200)

push!(RMS, RM)
push!(RMS, RM1)
push!(RMS, RM2)

env = Env(1, BLOCKS, RMS, [], 0.5)

PeriDyn.TIMEIT_REF[] = false

steps1 = 5000
velocity_verlet!([env], steps1; filewrite_freq=200, neigh_update_freq=100, out_dir="./output/2blocks")
#
