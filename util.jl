struct AbstractMaterial
    y::Array{Float64,2}
    velocity::Array{Float64,2}
    x::Array{Float64,2}
    particle_size::Float64
    volume::Array{Float64,1}
    density::Float64
    horizon::Float64
    critical_stretch::Float64
    family::Array{Int64,2}
    damage::BitArray{2}
end

abstract type PeridynamicsMaterial end

function Material(y0,v0,x,volume,density,horizon,critical_stretch;particle_size=0,max_neigh=50)
    family = cal_family(x,horizon,max_neigh)
    damage = family.>0
    if particle_size==0
        particle_size = volume[1]^(1/3)
    end
    return AbstractMaterial(y0,v0,x,particle_size,volume,density,horizon,critical_stretch,family,damage)
end



function horizon_correction(dr,r,hr)
    """
    Horizon correction for perdynamic model.
    Args
        dr vector_ij.
    """
    return 1
end

function influence_function(dr)
    return 1/s_magnitude(dr)
end

function dilatation(y,S::AbstractMaterial,m::Array{Float64,1})
    """
    Calculates dilatation for a point.
    Args
        y_i,x_i,info_i,m_i,V,x,y
    """
    theta = 0*S.volume
    for i in 1:size(S.x,1)
        for k in 1:size(S.family,2)
            j = S.family[i,k]
            if (j>0) & (i!=j)
                E = [i,j]
                e = s_magnitude(s_Y(E,y)) - s_magnitude(s_X(E,S.x))
                x = s_X(E,S.x)
                theta[i] += influence_function(x)*s_magnitude(x)*e*horizon_correction(x,S.particle_size,S.horizon)*S.volume[j]
            end
        end
    end
    return 3*theta./m
end


function weighted_volume(S::AbstractMaterial)
    m = 0*S.volume
    for i in 1:size(S.x,1)
        for k in 1:size(S.family,2)
            j = S.family[i,k]
            if (j>0) & (i!=j)
                dr = S.x[j,:]-S.x[i,:]
                m[i] += influence_function(dr)*s_magnitude(dr)^2*horizon_correction(dr,S.particle_size,S.horizon)*S.volume[j]
            end
        end
    end
    return m
end


function cal_family(x::Array{Float64,2},horizon::Float64,max_neigh::Int64)
    family = zeros(Int64,size(x,1),size(x,1))
    for i in 1:size(x,1)
        for j in (i+1):size(x,1)
            if s_magnitude(x[j,:].-x[i,:])<horizon
                family[i,j] = j
                family[j,i] = i
            end
        end
    end
    for i in 1:size(x,1)
        family[i,:] = sort(family[i,:])
    end
    return family[:,end:-1:end-max_neigh+1]
end


function write_data(filename,every,type,pos)
    file = open(filename, "w")
    for i in 1:every:size(pos,3)
        y = pos[:,:,i]
        N = size(y,1)
        write(file, "$N \n\n")
        for j in 1:size(y,1)
            t = type[j]
            a,b,c = y[j,:]
            write(file, "$t $a, $b, $c ")
            write(file,"\n")
        end
    end
    close(file)
end
