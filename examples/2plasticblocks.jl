using Pkg

Pkg.activate(".")

using PeriDyn
using PDMesh

using Random
Random.seed!(42)

x1, v1, y1, vol1, type1 = create(Cuboid([0 100; 0 100; 0 5]), resolution=1, type=1)

hor1 = 3.0
mat_gen1 = GeneralMaterial(y1, v1, x1, vol1, type1, hor1, max_neigh=200)

K, G = 1.0, 1.0
den1 = 1000.0
cstretch1 = 0.5
σy = 1.0
mat_spec1 = ElastoPlasticSolidSpecific([K], [G], [cstretch1], [den1], [σy])


block1 = PeridynamicsMaterial(mat_gen1, mat_spec1)

x2, v2, y2, vol2, type2 = create(move(Cuboid([0 10; 0 10; 0 10]), by=[45.0, 45, 7]), resolution=1, type=2)

v2[3,:] .-= 0.01

hor2 = 3.0
mat_gen2 = GeneralMaterial(y2, v2, x2, vol2, type2, hor2, max_neigh=200)

K2, G2 = 1.0, 1.0
den2 = 1000.0
cstretch2 = 0.5
mat_spec2 = OrdinaryStateBasedSpecific([K2], [G2], [cstretch2], [den2])

block2 = PeridynamicsMaterial(mat_gen2, mat_spec2)

k = 1000.0
RM = LinearRepulsionModel(k, block1, block2, distanceX=2, max_neighs=200)
RM1 = LinearRepulsionModel(k, block1, distanceX=2, max_neighs=200)
RM2 = LinearRepulsionModel(k, block2, distanceX=2, max_neighs=200)



env = Env(1, [block1, block2], [RM RM1 RM2], [], 0.5)

velocity_verlet!([env], 5000, filewrite_freq=10, neigh_update_freq=10, file_prefix="2blocks", start_at=0)

#
