using Profile
using PeriDyn 
const PD = PeriDyn


x1 = PD.create_block([1.0,1,1],[10,50,50])

v1 = x1*0
y1 = 1x1

mat_gen1 = PD.GeneralMaterial(y1,v1,x1,x1[1,:]*0 .+ 1, 1000.0, 3.0, 0.15, max_neigh=200)
mat_spec1 = PD.OrdinaryStateBasedSpecific(1.0,1.0,mat_gen1)
block1 = PD.OrdinaryStateBasedMaterial(1,mat_gen1,mat_spec1)

x2 = PD.create_block([1.0,1,1],[10,10,10])
x2[1,:] .+= 10.5
x2[2,:] .+= 20.5
x2[3,:] .+= 20.5
v2 = x2*0
v2[1,:] .-= 0.01
y2 = 1x2

mat_gen2 = PD.GeneralMaterial(y2,v2,x2,x2[1,:]*0 .+ 1, 10000.0, 3.1, 0.5, max_neigh=200)
mat_spec2 = PD.OrdinaryStateBasedSpecific(10.0,10.0,mat_gen2)
block2 = PD.OrdinaryStateBasedMaterial(2,mat_gen2,mat_spec2)



RM = PD.SimpleRepulsionModel(2.0,1000.0,block1,block2, distanceX=2,max_neighs=200)
RM1 = PD.SimpleRepulsionModel(2.0,1000.0,block1,distanceX=3,max_neighs=200)
RM2 = PD.SimpleRepulsionModel(2.0,1000.0,block2,distanceX=2,max_neighs=200)


env =  PD.Env(1,[block1,block2],[RM,RM1,RM2],[],0.2)


PD.velocity_verlet!([env],100,freq1=20,freq2=20,file_prefix="2blocks",start_at=0)






#
