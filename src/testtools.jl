function virgin_state_1D_2P(k)
    x1, v1, y1, vol1 = create_block([1.0,1,1],[2,2,2])
    y1 = k*y1

    den1 = 1000.0
    vol1 = 1.0
    hor1 = 3.0
    stretch1 = 0.2

    Es = 1.0
    nu = 0.2
    K1 = Es/3/(1-2nu)
    G1 = Es/2/(1+nu)


    mat_gen1 = GeneralMaterial(y1,v1,x1,x1[1,:]*0 .+ vol1, den1, hor1, stretch1, max_neigh=200)
    mat_spec1 = OrdinaryStateBasedSpecific(K1,G1,mat_gen1)
    block1 = OrdinaryStateBasedMaterial(1,mat_gen1,mat_spec1)

    epsilon = 1000.0
    alpha = 2.0

    RM11 = SimpleRepulsionModel(alpha,epsilon,block1, distanceX=2,max_neighs=200)

    dt = 0.2
    env =  Env(1,[block1],[RM11],[],dt)

    return env
end



function virgin_state(k)
    x1, v1, y1, vol1 = create_block([1.0,1,1],[2,2,2])
    y1 = k*y1

    den1 = 1000.0
    hor1 = 3.0
    stretch1 = 0.2

    Es = 1.0
    nu = 0.2
    K1 = Es/3/(1-2nu)
    G1 = Es/2/(1+nu)


    mat_gen1 = GeneralMaterial(y1,v1,x1,x1[1,:]*0 .+ vol1, den1, hor1, stretch1, max_neigh=200)
    mat_spec1 = OrdinaryStateBasedSpecific(K1,G1,mat_gen1)
    block1 = OrdinaryStateBasedMaterial(1,mat_gen1,mat_spec1)


    x2, v2, y2, vol2 = create_block([1.0,1,1],[2,2,2])
    x2[1,:] .+= 1.5
    y2 = k*y2
    y2[1,:] .+= 1.5


    den2 = 10000.0
    hor2 = 3.0
    stretch2 = 0.2
    K2 = K1
    G2 = G1

    mat_gen2 = GeneralMaterial(y2,v2,x2,x2[1,:]*0 .+ vol2, den2, hor2, stretch2, max_neigh=200)
    mat_spec2 = OrdinaryStateBasedSpecific(K2,G2,mat_gen2)
    block2 = OrdinaryStateBasedMaterial(2,mat_gen2,mat_spec2)


    epsilon = 1000.0
    alpha = 2.0

    RM12 = SimpleRepulsionModel(alpha,epsilon,block1,block2, distanceX=2,max_neighs=200)
    RM11 = SimpleRepulsionModel(alpha,epsilon,block1, distanceX=2,max_neighs=200)
    RM22 = SimpleRepulsionModel(alpha,epsilon,block2, distanceX=2,max_neighs=200)

    dt = 0.2
    env =  Env(1,[block1,block2],[RM12,RM11,RM22],[],dt)

    return env
end


function virgin_state2(k)
    epsilon = 1000.0
    alpha = 2.0

    x1, v1, y1, vol1 = create_block([1.0,1,1],[2,2,2])
    y1 = k*y1

    den1 = 1000.0
    hor1 = 3.0
    stretch1 = 0.2
    K1 = 1.0
    G1 = 1.0
    mat_gen1 = GeneralMaterial(y1,v1,x1,x1[1,:]*0 .+ vol1, den1, hor1, stretch1, max_neigh=200)
    mat_spec1 = BondBasedSpecific(epsilon)
    block1 = BondBasedMaterial(1,mat_gen1,mat_spec1)


    x2, v2, y2, vol2 = create_block([1.0,1,1],[2,2,2])
    x2[1,:] .+= 1.5
    y2 = k*y2
    y2[1,:] .+= 1.5


    den2 = 10000.0
    vol2 = 1.0
    hor2 = 3.0
    stretch2 = 0.2
    K2 = 1.0
    G2 = 1.0

    mat_gen2 = GeneralMaterial(y2,v2,x2,x2[1,:]*0 .+ vol2, den2, hor2, stretch2, max_neigh=200)
    mat_spec2 = BondBasedSpecific(epsilon)
    block2 = BondBasedMaterial(2,mat_gen2,mat_spec2)

    RM12 = SimpleRepulsionModel(alpha,epsilon,block1,block2, distanceX=2,max_neighs=200)
    RM11 = SimpleRepulsionModel(alpha,epsilon,block1, distanceX=2,max_neighs=200)
    RM22 = SimpleRepulsionModel(alpha,epsilon,block2, distanceX=2,max_neighs=200)

    dt = 0.2
    env =  Env(1,[block1,block2],[RM12,RM11,RM22],[],dt)

    return env
end
