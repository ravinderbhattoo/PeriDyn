

using Revise
using PeriDyn
using PDMesh

x1, v1, y1, vol1, type1 = create(Cuboid([0 20; 0 5; 0 5]), resolution=0.5)
mat_gen1 = GeneralMaterial(y1, v1, x1, vol1, type1, 1.5; max_neigh=200)

Es = 200
nu = 0.2
K = Es/3/(1-2nu)
G = Es/2/(1+nu)
den = 10000.0
cstretch = 0.15

mat_spec1 = BondBasedSpecific([K], [cstretch], [den])

block1 = PeridynamicsMaterial(mat_gen1, mat_spec1)


BC1 = FixBC(y1[1, :] .< 4.0)
BC2 = MoveBC(y1[1, :] .> 16.0, [0.01, 0.0, 0.0])

k = 1.0
RM1 = LinearRepulsionModel(k, block1; distanceX=3, max_neighs=200)

dt = 0.1
env = PeriDyn.Env(1, [block1], [RM1], [BC1, BC2], dt)

udf = Int(1/dt)
Steps = 20*udf

env.Params = Dict("left" => (y1[1,:] .< 4))

env.Out = Dict("y" => zeros(size(env.y)...,Steps))

env.Collect! = function (env, step)
    env.Out["y"][:, :, step] .= env.y
end

velocity_verlet!([env], Steps, filewrite_freq=udf,neigh_update_freq=udf,file_prefix="vv",start_at=0)

mat_spec2 = PairwiseNNSpecific([2, 10, 10, 1], [cstretch], [den])

block2 = PeridynamicsMaterial(mat_gen1, mat_spec2)

# env2 = PeriDyn.Env(1, [block2], [RM1], [BC1, BC2], dt)

FF = force_density_T(env.y, block1)

function loss(y)
    sum((FF .- force_density_T(y, block2)).^2)
end

# function big_loss(ys)
#     l = 0.0
#     for i in 1:10:size(ys, 3) 
#         l += loss(ys[:, :, i])
#     end
#     return l
# end

using Flux

theta = Flux.params(block2.specific.NNs[1])

opt = Flux.ADAM(0.001)

#for i in 1:10
#    println(i)
#    gs = Flux.gradient(()->loss(env.y), theta)
#    Flux.update!(opt, theta, gs)
#    println(loss(env.y))
#end

# velocity_verlet!([env2], Steps, filewrite_freq=udf,neigh_update_freq=udf,file_prefix="ml_vv",start_at=0)



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
