using Pkg

Pkg.activate(".")

using PeriDyn

using PDMesh
using Flux
using Plots
using Random
using Statistics
using Folds

Random.seed!(42)

resolution = 1.0

x, v, y, vol, type = unpack(create(Cuboid([-4 4; -4 4; -4 4]), resolution=resolution, type=1))

hor = 3.0*resolution
mat_gen = GeneralMaterial(y, v, x, vol, type, hor, max_neigh=200)

K, G = 100.0, 100.0
den = 1000.0
cstretch = 0.2

mat_spec = BondBasedSpecific([K], [cstretch], [den])


squareplus(x) = 0.5*(x + sqrt(x^2 + 4))
NNmat_spec = PairwiseNNSpecific([1,100,1], [cstretch], [den]; 
                act=squareplus)

block = PeridynamicsMaterial(mat_gen, mat_spec; name="block 1")
NNblock = PeridynamicsMaterial(mat_gen, NNmat_spec; name="block NN")

MSE(a, b) = mean((a .- b).^2)

M = size(y, 2)
loss = (k, particles, actual) -> begin
    pred = force_density_T(k*y, NNblock; particles)
    MSE(actual, pred)
end


ps = params(NNmat_spec.NNs[1,1])
opt = Flux.ADAM()

ydata = Dict()
for k in 0.85:0.05:1.15
    push!(ydata,  k => force_density_T(k*y, block))
end


using Zygote

function add!(a::Nothing, b::Zygote.Grads)
    b
end


function add!(a::Zygote.Grads, b::Zygote.Grads)
    for p in a.params
        a[p] .= a[p] .+ b[p]
    end
    a
end

function div!(a::Zygote.Grads, i::T) where T <: Real
    for p in a.params
        a[p] ./= i
    end
    a
end

function mul!(a::Zygote.Grads, i::T) where T <: Real
    for p in a.params
        a[p] .*= i
    end
    a
end


ks = [0.85, 0.9, 0.95, 1.05, 1.1, 1.15]
# ks = [0.95]

get_gs = (k) -> begin
    particles = randperm(M)[1:10]    
    actual = ydata[k]
    gs = gradient(() -> loss(k, particles, actual[:, particles]), ps)
end

for i in 1:1000
    
    gs = Folds.mapreduce(get_gs, add!, ks; init=nothing)
    Flux.update!(opt, ps, div!(gs, length(ks)))    
    
    if i%10==0
        l = 0.0
        fig = plot()
        for k in ks
            actual = ydata[k]
            pred = force_density_T(k*y, NNblock)
            l += MSE(actual, pred)
            scatter!(reshape(actual, (:,)), reshape(pred, (:,)))
        end
        println("Loss $i : ", l)
        display(fig)
    end
end

#