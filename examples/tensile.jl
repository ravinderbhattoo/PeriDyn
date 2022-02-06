using Pkg

Pkg.activate(".")

using Revise
using PeriDyn
using PDMesh

x1, v1, y1, vol1, type1 = unpack(create(Cuboid([0 20; 0 5; 0 5]), resolution=1))
mat_gen1 = GeneralMaterial(y1, v1, x1, vol1, type1, 3.0; max_neigh=200)

Es = 20
nu = 0.2
K = Es/3/(1-2nu)
G = Es/2/(1+nu)
den = 1000.0
cstretch = 0.2

mat_spec1 = OrdinaryStateBasedSpecific([K], [G], [cstretch], [den])

block1 = PeridynamicsMaterial(mat_gen1, mat_spec1)


BC1 = FixBC(y1[1, :] .< 4.0)
BC2 = MoveBC(y1[1, :] .> 16.0, [0.01, 0.0, 0.0])

k = 1.0
RM1 = LinearRepulsionModel(k, block1; distanceX=3, max_neighs=200)

env = PeriDyn.Env(1, [block1], [RM1], [BC1, BC2], 1.0)

Steps = 100

env.Params = Dict("left" => (y1[1,:] .< 4))

env.Out = Dict("Force" => zeros(3,Steps))

env.Collect! = function (env, step)
    env.Out["Force"][:, step] = sum(env.f[:, env.Params["left"]], dims=2)
end

quasi_static!([env], Steps, 0.1; max_iter=1000, 
    filewrite_freq=1, neigh_update_freq=1, out_dir="./output/tensile/minimize", start_at=0)

 