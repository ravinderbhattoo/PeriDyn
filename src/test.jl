using Profile

include("./mesh.jl")
include("./simulation.jl")
include("./util.jl")
include("./solvers.jl")
include("./operators.jl")
include("./materials/bond_based.jl")
include("./materials/ordinary_state_based.jl")
include("./contacts/contacts.jl")
include("./contacts/simple.jl")

x1 = create_block([1,1,1],[4,10,10])
v1 = x1*0
y1 = 1x1

mat_gen1 = Material(y1,v1,x1,x1[:,1]*0 .+ 1, 1000.0, 3.1, 0.5, max_neigh=100)
mat_spec1 = OrdinaryStateBasedSpecific(10.0,10.0,mat_gen1)
block1 = OrdinaryStateBasedMaterial(1,mat_gen1,mat_spec1)

x2 = create_block([1,1,1],[4,4,4])
x2[:,1] .+= 5
x2[:,2:3] .+= 3
v2 = x2*0
v2[:,1] .+= -0.01
y2 = 1x2

mat_gen2 = Material(y2,v2,x2,x2[:,1]*0 .+ 1, 1000.0, 3.1, 0.5, max_neigh=8)
mat_spec2 = OrdinaryStateBasedSpecific(10.0,10.0,mat_gen2)
block2 = OrdinaryStateBasedMaterial(2,mat_gen2,mat_spec2)

RM = SimpleRepulsionModel(2.0,100.0,block1,block2,distanceX=2)

env = Env([block1,block2],[RM],1500.0,0.2)

pos = velocity_verlet(env)


#
