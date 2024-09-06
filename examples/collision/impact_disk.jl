using Revise
using PeriDyn
using PDMaterialPoints
using CUDA

OUTDIR = joinpath(homedir(), "Downloads", "PeriDyn", "Collision", "DiskImpact", "ElastoPlastic")

device = :cpu
PeriDyn.set_device(device)
PDMaterialPoints.set_device(device)

##############################################################################
# Block 1
##############################################################################
resolution = 1.0
x1, v1, y1, vol1, type1 = unpack(create(move(Sphere(5.0); by=[0.0, 0.0, 15.0]); 
                                    resolution=resolution, type=1))

v1 = Array(v1)
velocity = 1.0
v1[3, :] .-= velocity

horizon1 = 3.0*resolution
mat_gen1 = GeneralMaterial(y1, v1, x1, vol1, type1, horizon1, max_neigh=200)


Es = 25 # GPa
nu = 0.15
K1 = Es/3/(1-2nu)
G1 = Es/2/(1+nu)

density1 = 10.0
cstretch1 = 1.0

mat_spec1 = BondBasedSpecific([K1], [cstretch1], [density1])

block1 = PeridynamicsMaterial(mat_gen1, mat_spec1; name="block 1")

##############################################################################
# Block 2
##############################################################################
resolution = 1.0
x2, v2, y2, vol2, type2 = unpack(create(Disk(50, 5); 
                                    resolution=resolution, type=2))
horizon2 = 3.0*resolution
mat_gen2 = GeneralMaterial(y2, v2, x2, vol2, type2, horizon2, max_neigh=200)

Es = 25.0 # GPa
nu = 0.15
K2 = Es/3/(1-2nu)
G2 = Es/2/(1+nu)

density2 = 1.0

# cstretch2 = 0.01
# mat_spec2 = BondBasedSpecific([K2], [cstretch2], [density2])

cstretch2 = 0.2
sigma = 0.1*Es
mat_spec2 = ElastoPlasticSolidSpecific([K2], [G2], [cstretch2], [density2], [sigma])

type_ = 2
block2 = PeridynamicsMaterial(mat_gen2, mat_spec2; name="block $type_")


##############################################################################
# Contact Model
##############################################################################
k = 1.0e3 # spring stiffness
max_neighs = 400

RM1 = LinearSpringContactModel(k, block1, block2, distanceX=2, max_neighs=max_neighs)
# RM1 = LinearSpringContactModel(k/1000, block1, distanceX=2, max_neighs=max_neighs)
# RM2 = LinearSpringContactModel(k/1000, block2, distanceX=2, max_neighs=max_neighs)

##############################################################################
# Simulation Env
##############################################################################
dt = velocity / 100
env = Env(1, [block1, block2], [RM1], [], 0.01)

##############################################################################
# Solver
##############################################################################
solver = DSVelocityVerlet()
Steps, fwf, nuf = 4000, 100, 10


run!([env], Steps, solver;
        filewrite_freq=fwf, neigh_update_freq=nuf, 
        out_dir=OUTDIR, start_at=0, ext=:data)

println("Simulation Finished. :)")

















# const Es = 2000.0
# const nu = 0.25 # for bond based peridynamics

# K = Es/3/(1-2nu)
# G = Es/2/(1+nu)

# const rho = 2000.0

# wv = sqrt(Es/rho)
# println("Wave velocity: ", wv)

# const cs = 0.01
# const reso = 0.5
# const horizon = 3*reso

# time_step = 1.0 * reso / wv

# C = 18K/(pi*horizon^4)

# gen_mat = PDBenchmark.NameParam(:GeneralMaterial, (horizon), Dict(:max_neigh=>200))
# spc_mat = [
#             PDBenchmark.NameParam(:BondBasedSpecific, ([C], [cs], [rho], ), Dict()),
#             PDBenchmark.NameParam(:BondBasedSpecific, ([C], [10*cs], [10*rho], ), Dict())
#           ]

# solver = :vv

# st_st = Dict(
#     :qs => Dict(
#             :printevery => 1,
#             :fwf => 1,
#             :Steps => 20,
#             ),
#     :vv => Dict(
#             :printevery => 1,
#             :fwf => 10,
#             :nuf => 10,
#             :Steps => 4000,
#             )
#     )

# out_dir = "ImpactDiskBB_thin4" * string(solver)

# try
#     foreach(rm, filter(endswith(".data"), readdir("./output/"*out_dir, join=true)))
# catch
#     nothing
# end

# test = PDBenchmark.ImpactDisk(;gen_mat=gen_mat, dt=time_step, st_st[solver]..., projectile_velocity=[0.0, 0.0, -0.05],
#                                 spc_mat=spc_mat, resolution=reso, solver_=solver, out_dir=out_dir, makeplot=true, trueE=Es)

# env, solvef! = PDBenchmark.stage!(test)

# solvef!(env)


# # using JLD
# # using Plots
# # using Statistics

# # damage1 = Float64[]
# # damage2 = Float64[]
# # time = Float64[]
# # for i in 0:10:4000
# #     data = load("./output/ImpactDiskBB_vv/env_1_step_$(i).jld2")
# #     type_ = data["type"]
# #     damage_ = data["damage"]
# #     mask1 = type_ .== 1
# #     mask2 = type_ .== 2

# #     push!(damage1, mean(data["damage"][mask1]))
# #     push!(damage2, mean(data["damage"][mask2]))
# #     push!(time, i)
# # end

# # fig = plot(time, damage1; label="Disk")
# # plot!(time, damage2; label="Projectile")
# # xlabel!("Time step")
# # ylabel!("Average damage")
# # savefig("./output/damage_disk.png")

# # using PeriDyn
# # PeriDyn.jld2ovito("./output/$(out_dir)/env_1_step_*.jld2", st_st[solver][:Steps]; start=0, step=Int64(st_st[solver][:Steps]/100))

