# Activate PeriDyn environment from root dir
using Pkg
Pkg.activate(".")

# imports 
using Revise
using PeriDyn
using PDMesh
using Random

set_multi_threading(true)

# Create block
Random.seed!(42)
resolution = 1
x1, v1, y1, vol1, type1 = unpack(create(Cuboid([-5 5; -5 5; 0 50]), resolution=resolution, type=1))

hor1 = 3.0*resolution
mat_gen1 = GeneralMaterial(y1, v1, x1, vol1, type1, hor1, max_neigh=200)

K, G = 1000.0, 1000.0
den1 = 1000.0
cstretch1 = 1.0 
mat_spec1 = OrdinaryStateBasedSpecific([K], [G], [cstretch1], [den1])

block1 = PeridynamicsMaterial(mat_gen1, mat_spec1; name="block 1")

BLOCKS = [block1]

RMs = [] # no repulsive models
dt = 0.5 # time step
env_id = 1

# boundary conditions
# 1. Fix on one side
# 2. Move on other side

velocity = 0.05*resolution

bool_left = x1[3, :] .< 5.0
bool_right = x1[3, :] .> 45.0
bool_mid = (x1[3, :] .> 22.0) .* (x1[3, :] .< 28.0) 
move_speed = [0.0, 0.0, velocity] 

bool_all = ones(Bool, length(bool_left))

BCs = [
    # MoveBC(bool_left, 0.1*move_speed),
    # ToFroBC(bool_right, move_speed, 40),
    # DeltaScaleBC(bool_mid, [0.0, 0.0, 0.001], [0.0, 0.0, 25.0])
    ScaleFixWaitBC(bool_left .+ bool_right, [0.0, 0.0, 0.1*velocity], 
                    [0.0, 0.0, 25.0], 1, bool_all)
    ]

env = Env(env_id, BLOCKS, RMs, BCs, dt)

steps1 = 500
velocity_verlet!([env], steps1; filewrite_freq=10, neigh_update_freq=10, out_dir="./output/boundary_conditions/fix_move_vv/")


steps1 = 20
lr = 0.001
quasi_static!([env], steps1, lr; filewrite_freq=1, neigh_update_freq=1, out_dir="./output/boundary_conditions/fix_move_qs/")



#
