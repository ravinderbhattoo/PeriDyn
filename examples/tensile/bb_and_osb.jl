using Pkg
Pkg.activate(".")

using PeriDyn
using PDMaterialPoints

PeriDyn.set_cuda(true)
# PDMaterialPoints.set_cuda(true)

resolution = 0.5

x1, v1, y1, vol1, type1 = unpack(create(Cuboid([0 20; 0 5; 0 5]), resolution=resolution))  # mm
mat_gen1 = GeneralMaterial(y1, v1, x1, vol1, type1, 3.0*resolution; max_neigh=200)

Es = 70 # GPa
nu = 0.15
K = Es/3/(1-2nu)
G = Es/2/(1+nu)
den = 2.2 * 1000.0 # Kg/m3
cstretch = 0.15

# mat_spec1 = BondBasedSpecific([K], [cstretch], [den])
mat_spec1 = OrdinaryStateBasedSpecific([K], [G], [cstretch], [den])

block1 = PeridynamicsMaterial(mat_gen1, mat_spec1)

# solver = QSDrag(0.1, 0.1; max_iter=2000, x_tol=0.001, f_tol=0.000001)
# Steps, fwf, nuf = 20, 1, 1

solver = DSVelocityVerlet()
Steps, fwf, nuf = 2000, 100, 10

dt = 1.0

BC1 = FixBC(y1[1, :] .< 4.0)
BC2 = MoveBC(y1[1, :] .> 16.0, [0.05*20/Steps/dt, 0.0, 0.0])

# BC2 = DeltaScaleBC(y1[1, :] .> 16.0, [1.01, 1.0, 1.0], [4.0, 0.0, 0.0])
# BC2 = ScaleFixWaitBC(y1[1, :] .> 16.0, [0.01, 0.0, 0.0], [0.0, 0.0, 0.0], 100, y1[1, :] .> -1)


RMs = []
k = 1.0
RM1 = LinearRepulsionModel(k, block1; distanceX=3, max_neighs=200)
RMs = [RM1]

env = PeriDyn.Env(1, [block1], RMs, [BC1, BC2], dt)

env.Params = Dict("left" => (y1[1,:] .< 4))

env.Out = Dict("Force" => zeros(3,Steps))

env.Collect! = function (env, step)
    env.Out["Force"][:, step] = sum(env.f[:, env.Params["left"]], dims=2)
end


out_dir="./output/tensile_sim_OSB/$(typeof(solver))"

run!([env], Steps, solver;
    filewrite_freq=fwf, neigh_update_freq=nuf, out_dir=out_dir, start_at=0, ext=:jld)

PeriDyn.write_data("$(out_dir)/env_Out.jld"; Out=env.Out)

PeriDyn.jld2ovito("$(out_dir)/env_1_step_*.jld", Steps; start=0, step=fwf)



