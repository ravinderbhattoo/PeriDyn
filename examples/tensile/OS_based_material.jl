using PeriDyn

function assign_mat(K, G, sigma_y, cstretch, density)
    return PeriDyn.OrdinaryStateBasedSpecific([K], [G], [cstretch], [density])
end

include(joinpath(@__DIR__, "bar.jl"))