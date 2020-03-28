using Profile
using BenchmarkTools

include("./mesh.jl")
include("./simulation.jl")
include("./operators.jl")
include("./materials/material.jl")
include("./materials/bond_based.jl")
include("./materials/ordinary_state_based.jl")
include("./contacts/contacts.jl")
include("./contacts/simple.jl")
include("./util.jl")
include("./global_bc.jl")
include("./solvers.jl")



x1 = create_block([1.0,1,1],[20,10,10])
x1 .+= -1
v1 = x1*0
y1 = 1x1



mat_gen1 = GeneralMaterial(y1,v1,x1,x1[1,:]*0 .+ 1, 1000.0, 3.0, 0.5, max_neigh=200)


Es = 20
nu = 0.2
K = Es/3/(1-2nu)
G = Es/2/(1+nu)

mat_spec1 = OrdinaryStateBasedSpecific(K,G,mat_gen1)
block1 = OrdinaryStateBasedMaterial(1,mat_gen1,mat_spec1)


bool4 = y1[1,:].>-10000
bool5 = (y1[1,:].<6) .| (y1[1,:].>15)

BC4 = ScaleFixBC(1,bool4,[0.02,0,0],bool5)

RM1 = SimpleRepulsionModel(2.0,1.0,block1, distanceX=3,max_neighs=200,)

env =  Env(1,[block1],[RM1],[BC4],1)


Steps = 30

env.Params = Dict("left" => (y1[1,:].<6))
env.Out = Dict("Force" => zeros(3,Steps))
env.Collect! = function (Params,Out,step)
    Out["Force"][:,step] = sum(env.f[:,Params["left"]],dims=2)
end

quasi_static!([env],Steps,10.0,freq1=1,freq2=1,file_prefix="minimize",start_at=0)


using Plots


dl = 8.1
plot([1:30]*dl/20/30,10*env.Out["Force"][1,:],marker=4,linewidth=6,label=raw"\sigma_x")
plot!([1:30]*dl/20/30,10*env.Out["Force"][2,:],marker=4,linewidth=6,label=raw"\sigma_y")
plot!([1:30]*dl/20/30,10*env.Out["Force"][3,:],marker=4,linewidth=6,label=raw"\sigma_z")
xlabel!("Strain")
ylabel!("Stress")



using DelimitedFiles
mask = y1[1,:].<6
F = zeros(30)
for i in 1:30
    df = readdlm("./output/minimize_env_1_step_$i.data", ',', Float64, '\n', skipstart=2)
    F[i] = sum(df[mask,end-2])
end

x = [i for i in 1:30]/30*8.1/20
y = 10*F
plot(x, y, marker=4, linewidth=6, label=raw"\sigma_x")
xlabel!("Strain")
ylabel!("Stress")

E_ = (y[3]-y[2])/(x[3]-x[2])

#
