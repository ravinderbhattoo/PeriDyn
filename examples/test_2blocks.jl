using Profile
using PeriDyn

x1, v1, y1, vol1 = create_block([1.0,1,1],[100,100,10])
den1 = 1000.0
hor1 = 3.0
cstretch1 = 0.15

mat_gen1 = GeneralMaterial(y1, v1, x1, vol1, den1, hor1, cstretch1, max_neigh=200)

K, G = 1.0, 1.0
mat_spec1 = OrdinaryStateBasedSpecific(K,G,mat_gen1)

id = 1::Int64
block1 = OrdinaryStateBasedMaterial(id, mat_gen1, mat_spec1)

x2, v2, y2, vol2 = create_block([1.0,1,1],[10,10,10])

x2[1,:] .+= 10.5
x2[2,:] .+= 20.5
x2[3,:] .+= 20.5

v2[1,:] .-= 0.01

den2 = 10000.0
hor2 = 3.0
cstretch2 = 0.5

mat_gen2 = GeneralMaterial(y2, v2, x2, vol2, den2, hor2, cstretch2, max_neigh=200)

K2, G2 = 10.0, 10.0

mat_spec2 = OrdinaryStateBasedSpecific(K2, G2, mat_gen2)

id = 2::Int64
block2 = OrdinaryStateBasedMaterial(id, mat_gen2, mat_spec2)


alpha = 2.0
epsilon = 1000.0

RM = SimpleRepulsionModel(alpha, epsilon, block1, block2, distanceX=2, max_neighs=200)
RM1 = SimpleRepulsionModel(alpha, epsilon, block1, distanceX=3, max_neighs=200)
RM2 = SimpleRepulsionModel(alpha, epsilon, block2, distanceX=2, max_neighs=200)
 

velocity_verlet!([env],1000,freq1=20,freq2=20,file_prefix="2blocks",start_at=0)

#
