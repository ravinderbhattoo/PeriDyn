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

function _influence_function(dr)
    return 1/s_magnitude(dr)
end

function dilatation(y,S::AbstractMaterial,m::Array{Float64,1})
    """
    Calculates dilatation for a point.
    Args
        y_i,x_i,info_i,m_i,V,x,y
    """
    theta = 0*S.volume::Array{Float64}
    j = 1::Int64
    E_n = 0.0::Float64
    E = 0.0::Float64
    x = [1.0,0,0]::Array{Float64}
    e = 0.0::Float64
    for i in 1:size(S.x,2)
        for k in 1:size(S.family,1)
            if S.intact[k,i]
                j = S.family[k,i]
                if (j>0) & (i!=j)
                    x[1],x[2],x[3] = S.x[1,j]-S.x[1,i],S.x[2,j]-S.x[2,i],S.x[3,j]-S.x[3,i]
                    E_n = sqrt((y[1,j]-y[1,i])^2 + (y[2,j]-y[2,i])^2 + (y[3,j]-y[3,i])^2)
                    E = sqrt((S.x[1,j]-S.x[1,i])^2 + (S.x[2,j]-S.x[2,i])^2 + (S.x[3,j]-S.x[3,i])^2)
                    e = E_n - E
                    theta[i] += influence_function(x)*E*e*S.volume[j] *horizon_correction(x,S.particle_size,S.horizon)
                end
            end
        end
    end
    return 3*theta./m
end


function weighted_volume(S::AbstractMaterial)
    m = 0*S.volume::Array{Float64,1}
    dr = zeros(3)::Array{Float64,1}
    j = 0::Int64
    for i in 1:size(S.x,2)
        for k in 1:size(S.family,1)
            j = S.family[k,i]
            if (j>0) & (i!=j)
                dr[1] = S.x[1,j]-S.x[1,i]
                dr[2] = S.x[2,j]-S.x[2,i]
                dr[3] = S.x[3,j]-S.x[3,i]
                m[i] += influence_function(dr)*s_magnitude(dr)^2*horizon_correction(dr,S.particle_size,S.horizon)*S.volume[j]
            end
        end
    end
    return m
end



function cal_family(x::Array{Float64,2},horizon::Float64, max_neigh::Int64)::Array{Float64,2}
    family = zeros(Int64,size(x,2),size(x,2))
    for i in 1:size(x,2)
        a1,b1,c1 = x[1,i],x[2,i],x[3,i]
        for j in (i+1):size(x,2)
            a2,b2,c2 = x[1,j],x[2,j],x[3,j]
            if (a1-a2)*(a1-a2)+(b1-b2)*(b1-b2)+(c1-c2)*(c1-c2)<horizon^2
                family[i,j] = i
                family[j,i] = j
            end
        end
    end
    family = sort(family,dims=1)
    return family[end:-1:end+1-min(size(x,2),max_neigh),1:end]
end


function write_data(filename,every,type,pos)
    file = open(filename, "w+")
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

function write_data(filename::String,type::Array{Int64,1},y::Array{Float64,2})
    file = open(filename, "w+")
    N = size(y,2)
    write(file, "$N \n\n")
    for j in 1:size(y,2)
        t = type[j]
        a,b,c = y[1,j],y[2,j],y[3,j]
        write(file, "$t $a, $b, $c ")
        write(file,"\n")
    end
    close(file)
end

function write_data(filename::String,type::Array{Int64,1},y::Array{Float64,2}, v::Array{Float64,2},f::Array{Float64,2},)
    file = open(filename, "w+")
    N = size(y,2)
    write(file, "$N \n\n")
    for j in 1:size(y,2)
        t = type[j]
        a,b,c = y[1,j],y[2,j],y[3,j]
        a1,b1,c1 = v[1,j],v[2,j],v[3,j]
        a2,b2,c2 = f[1,j],f[2,j],f[3,j]
        write(file, "$t $a, $b, $c $a1, $b1, $c1 $a2, $b2, $c2 ")
        write(file,"\n")
    end
    close(file)
end
