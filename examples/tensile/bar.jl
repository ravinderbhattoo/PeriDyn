using Revise
using PeriDyn
using PDMaterialPoints
using CUDA

device = :cpu
PeriDyn.set_device(device)
PDMaterialPoints.set_device(device)

##############################################################################
# Block 1
##############################################################################
resolution = 0.5

obj = Cuboid([-10 10; -5 5; -5 5])

function func(out)
    x = out[:x]
    mask = @. (x[1, :] < 0.5) & (x[1, :] > -0.5) & (x[3, :] > 3.0)
    mask
end
obj = delete(obj, func)

x1, v1, y1, vol1, type1 = unpack(create(obj; 
                                    resolution=resolution, type=1))
horizon1 = 3.0*resolution
mat_gen1 = GeneralMaterial(y1, v1, x1, vol1, type1, horizon1, max_neigh=200)

E = 3.0
nu = 0.15
K = E / 3 / (1-2*nu) 
G = E / 2 / (1+nu)
density1 = 1.0
cstretch1 = 0.2
sigma_y = cstretch1*E / 2


mat_spec1 = assign_mat(K, G, sigma_y, cstretch1, density1)

block1 = PeridynamicsMaterial(mat_gen1, mat_spec1; name="block 1")

OUTDIR = joinpath(homedir(), "Downloads", "PeriDyn", "Tensile", "$(typeof(mat_spec1))")

##############################################################################
# Boundary condition
##############################################################################
speed = 0.01
BC1 = FixBC(y1[1, :] .< -7.0)
BC2 = ToFroBC(y1[1, :] .> 7.0, [1.0, 0.0, 0.0]*speed, Inf)
BCs = [BC1, BC2]

##############################################################################
# Contact Model
##############################################################################
k = 1.0e1 # spring stiffness
max_neighs = 400

RM1 = LinearSpringContactModel(k, block1, distanceX=2, max_neighs=max_neighs)

##############################################################################
# Simulation Env
##############################################################################
env = Env(1, [block1], [RM1], BCs, 0.1)

##############################################################################
# Solver
##############################################################################
solver = DSVelocityVerlet()
Steps, fwf, nuf, apf = 4000, 100, 10, 100

run!([env], Steps, solver;
        average_prop_freq=apf,
        filewrite_freq=fwf, 
        neigh_update_freq=nuf, 
        out_dir=OUTDIR, 
        start_at=0, ext=:data)

println("Simulation Finished. :)")
