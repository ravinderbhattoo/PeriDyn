using Revise
using PeriDyn
using PDMesh

x1, v1, y1, vol1, type1 = create(Cuboid([0 20; 0 5; 0 5]), resolution=1)
mat_gen1 = GeneralMaterial(y1, v1, x1, vol1, type1, 3.0; max_neigh=200)

Es = 20
nu = 0.2
K = Es/3/(1-2nu)
G = Es/2/(1+nu)
den = 1000.0
cstretch = 0.15

mat_spec1 = OrdinaryStateBasedSpecific([K], [G], [cstretch], [den])

block1 = PeridynamicsMaterial(mat_gen1, mat_spec1)


BC1 = FixBC(y1[1, :] .< 4.0)
BC2 = MoveBC(y1[1, :] .> 16.0, [0.01, 0.0, 0.0])

k = 1.0
RM1 = LinearRepulsionModel(k, block1; distanceX=3, max_neighs=200)

env =  PD.Env(1, [block1], [RM1], [BC1, BC2], 1.0)

Steps = 20

env.Params = Dict("left" => (y1[1,:] .< 4))

env.Out = Dict("Force" => zeros(3,Steps))

env.Collect! = function (Params, Out, step)
    Out["Force"][:, step] = sum(env.f[:, Params["left"]], dims=2)
end

quasi_static!([env],Steps,10.0,filewrite_freq=1,neigh_update_freq=1,file_prefix="minimize",start_at=0)


# using Plots
# dl = 8.1
# plot([1:10]*dl/20/10,10*env.Out["Force"][1,:],marker=4,linewidth=6,label=raw"\sigma_x")
# plot!([1:10]*dl/20/10,10*env.Out["Force"][2,:],marker=4,linewidth=6,label=raw"\sigma_y")
# plot!([1:10]*dl/20/10,10*env.Out["Force"][3,:],marker=4,linewidth=6,label=raw"\sigma_z")
# xlabel!("Strain")
# ylabel!("Stress")


#
# using DelimitedFiles
# mask = y1[1,:].<6
# F = zeros(10)
# for i in 1:10
#     df = readdlm("./output/minimize_env_1_step_$i.data", ',', Float64, '\n', skipstart=2)
#     F[i] = sum(df[mask,end-2])
# end
#
# x = [i for i in 1:10]/3*8.1/20
# y = 10*F
# plot(x, y, marker=4, linewidth=6, label=raw"\sigma_x")
# xlabel!("Strain")
# ylabel!("Stress")
#
# E_ = (y[3]-y[2])/(x[3]-x[2])
