using Pkg
Pkg.activate(".")

using Revise
using PeriDyn

out_dir="./output/ttensile_sim_OSB/DSVelocityVerlet"

Threads.@threads for i in 0:100:2000
    PeriDyn.jld2ovito("$(out_dir)/env_1_step_*.jld", i; start=i)
end