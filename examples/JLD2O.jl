using PeriDyn

data_dir="./output/ttensile_sim_OSB/DSVelocityVerlet"

Threads.@threads for i in 0:100:2000
    PeriDyn.jld2ovito("$(data_dir)/env_1_step_*.jld2", i; start=i)
end