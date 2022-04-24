using Pkg

Pkg.activate(".")

using Revise
using Flux
using Statistics

f = (x) -> x^5 - 4*x^4 + 3*x^2 + 4

function data(N)
    N = N // 4    
    x_ = -2.0 : 1/N : 2.0
    (collect(x_), f.(x_))
end


x, y = data(100)

squareplus(x) = 0.5*(x + sqrt(x^2 + 4))


NN = Chain((x)->[x], Dense(1, 5, squareplus), Dense(5, 5, squareplus), Dense(5, 1), x -> first(x) )


loss = () -> begin
    mean((y .- NN.(x)).^2)
end

ps = params(NN)

opt = Flux.ADAM()



for i in 1:10000
    
    gs = gradient(loss, ps)
    Flux.update!(opt, ps, gs)
    println("Loss $i : ", loss())
    
    if i%100 == 0
        fig=plot(x, NN.(x); :linewidth=>3)
        scatter!(x, y)
        display(fig)
    end

end

#