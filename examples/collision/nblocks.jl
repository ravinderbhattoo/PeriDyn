using Revise
using PeriDyn
using PDMaterialPoints
using CUDA

OUTDIR = joinpath(homedir(), "Downloads", "PeriDyn", "Collision", "2Blocks")

device = :cpu
PeriDyn.set_device(device)
PDMaterialPoints.set_device(device)

##############################################################################
# Block 1
##############################################################################
resolution = 0.5
x1, v1, y1, vol1, type1 = unpack(create(Cuboid([-5 5; -5 5; -5 5]); 
                                    resolution=resolution, type=1))
horizon1 = 3.0*resolution
mat_gen1 = GeneralMaterial(y1, v1, x1, vol1, type1, horizon1, max_neigh=200)

K, G = 1.0, 1.0
density1 = 1.0
cstretch1 = 0.5
mat_spec1 = BondBasedSpecific([K], [cstretch1], [density1])

block1 = PeridynamicsMaterial(mat_gen1, mat_spec1; name="block 1")


##############################################################################
# Block 2
##############################################################################
resolution = 1.0
x2, v2, y2, vol2, type2 = unpack(create(Cuboid([-5 5; -17 -7; -5 5]); 
                                    resolution=resolution, type=2))
v2 = Array(v2)
v2[2, :] .+= 0.0
horizon2 = 3.0*resolution
mat_gen2 = GeneralMaterial(y2, v2, x2, vol2, type2, horizon2, max_neigh=200)

K2, G2 = 1.0, 1.0
density2 = 1.0
cstretch2 = 0.5

mat_spec2 = BondBasedSpecific([K2], [cstretch2], [density2])

type_ = 2
block2 = PeridynamicsMaterial(mat_gen2, mat_spec2; name="block $type_")


##############################################################################
# Block 3
##############################################################################
resolution = 0.5
x3, v3, y3, vol3, type3 = unpack(create(Cuboid([-5 5; -29 -19; -5 5]); 
                                    resolution=resolution, type=3))
v3 = Array(v3)
v3[2, :] .+= 0.01
horizon3 = 3.0*resolution
mat_gen3 = GeneralMaterial(y3, v3, x3, vol3, type3, horizon3, max_neigh=200)

K3, G3 = 1.0, 1.0
density3 = 1.0
cstretch3 = 0.5

mat_spec3 = BondBasedSpecific([K3], [cstretch3], [density3])

type_ = 3
block3 = PeridynamicsMaterial(mat_gen3, mat_spec3; name="block $type_")

##############################################################################
# Contact Model
##############################################################################
k = 1.0e1 # spring stiffness
max_neighs = 400

RM1 = LinearSpringContactModel(k, block1, block2, distanceX=2, max_neighs=max_neighs)
RM2 = LinearSpringContactModel(k, block2, block3, distanceX=2, max_neighs=max_neighs)
# RM1 = LinearSpringContactModel(k/1000, block1, distanceX=2, max_neighs=max_neighs)
# RM2 = LinearSpringContactModel(k/1000, block2, distanceX=2, max_neighs=max_neighs)

##############################################################################
# Simulation Env
##############################################################################
env = Env(1, [block1, block2, block3], [RM1, RM2], [], 0.1)

##############################################################################
# Solver
##############################################################################
solver = DSVelocityVerlet()
Steps, fwf, nuf = 8000, 100, 10


run!([env], Steps, solver;
        filewrite_freq=fwf, neigh_update_freq=nuf, 
        out_dir=OUTDIR, start_at=0, ext=:data)

println("Simulation Finished. :)")
