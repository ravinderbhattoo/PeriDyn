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

x1 = create_block([1.0,1,1],[40,20,10])
notch_mask = (19 .< x1[1,:] .< 23) .& ((x1[2,:].<6) .| (x1[2,:].>15))
x1 = x1[:,.~notch_mask]

v1 = x1*0
y1 = 1x1

mat_gen1 = GeneralMaterial(y1,v1,x1,x1[1,:]*0 .+ 1, 1000.0, 3.0, 0.5, max_neigh=200)
mat_spec1 = OrdinaryStateBasedSpecific(10.0,10.0,mat_gen1)
block1 = OrdinaryStateBasedMaterial(1,mat_gen1,mat_spec1)

bool1 = y1[1,:].<10
bool2 = y1[1,:].>30
BC1 = Move_BC(1,bool1,[0.0,0,0],1y1[:,bool1])
BC2 = Move_BC(1,bool2,[.5,0,0],1y1[:,bool2])


lazzy_mask = (40 .< x1[1,:] .< 60)
RM1 = SimpleRepulsionModel(2.0,1.0,block1,distanceX=3,max_neighs=200,)

env =  Env(1,[block1],[RM1],[BC1,BC2],1)


quasi_static!([env],50,10.0,freq1=1,freq2=1,file_prefix="minimize",start_at=0)





#
