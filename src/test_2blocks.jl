using Profile


include("./mesh.jl")
include("./simulation.jl")
include("./operators.jl")
include("./materials/material.jl")
include("./materials/bond_based.jl")
include("./materials/ordinary_state_based.jl")
include("./contacts/contacts.jl")
include("./contacts/simple.jl")
include("./util.jl")
include("./solvers.jl")
include("./global_bc.jl")

x1 = create_block([1.0,1,1],[10,50,50])

v1 = x1*0
y1 = 1x1

mat_gen1 = GeneralMaterial(y1,v1,x1,x1[1,:]*0 .+ 1, 1000.0, 3.0, 0.15, max_neigh=200)
mat_spec1 = OrdinaryStateBasedSpecific(1.0,1.0,mat_gen1)
block1 = OrdinaryStateBasedMaterial(1,mat_gen1,mat_spec1)

x2 = create_block([1.0,1,1],[10,10,10])
x2[1,:] .+= 10.5
x2[2,:] .+= 20.5
x2[3,:] .+= 20.5
v2 = x2*0
v2[1,:] .-= 0.01
y2 = 1x2

mat_gen2 = GeneralMaterial(y2,v2,x2,x2[1,:]*0 .+ 1, 10000.0, 3.1, 0.5, max_neigh=200)
mat_spec2 = OrdinaryStateBasedSpecific(10.0,10.0,mat_gen2)
block2 = OrdinaryStateBasedMaterial(2,mat_gen2,mat_spec2)



RM = SimpleRepulsionModel(2.0,1000.0,block1,block2, distanceX=2,max_neighs=200)
RM1 = SimpleRepulsionModel(2.0,1000.0,block1,distanceX=3,max_neighs=200)
RM2 = SimpleRepulsionModel(2.0,1000.0,block2,distanceX=2,max_neighs=200)


env =  Env(1,[block1,block2],[RM,RM1,RM2],[],0.2)


velocity_verlet!([env],100,freq1=20,freq2=20,file_prefix="2blocks",start_at=0)






#
