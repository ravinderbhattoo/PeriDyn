using Pkg
Pkg.activate(".")

using Revise
using PDBenchmark

PDBenchmark.set_cuda(true)
PDBenchmark.set_multi_threading(true)


const Es = 2000.0
const nu = 0.25 # for bond based peridynamics

K = Es/3/(1-2nu)
G = Es/2/(1+nu)

const rho = 2000.0

wv = sqrt(Es/rho)
println("Wave velocity: ", wv)

const cs = 0.01
const reso = 0.5
const horizon = 3*reso

time_step = 1.0 * reso / wv

C = 18K/(pi*horizon^4)

gen_mat = PDBenchmark.NameParam(:GeneralMaterial, (horizon), Dict(:max_neigh=>200))
spc_mat = [
            PDBenchmark.NameParam(:BondBasedSpecific, ([C], [cs], [rho], ), Dict()),
            PDBenchmark.NameParam(:BondBasedSpecific, ([C], [10*cs], [10*rho], ), Dict())
          ]

solver = :vv

st_st = Dict(
    :qs => Dict(
            :printevery => 1,
            :fwf => 1,
            :Steps => 20,
            ),
    :vv => Dict(
            :printevery => 1,
            :fwf => 10,
            :nuf => 10,
            :Steps => 4000,
            )
    )

out_dir = "ImpactDiskBB_thin4" * string(solver)

try
    foreach(rm, filter(endswith(".data"), readdir("./output/"*out_dir, join=true)))
catch
    nothing
end

test = PDBenchmark.ImpactDisk(;gen_mat=gen_mat, dt=time_step, st_st[solver]..., projectile_velocity=[0.0, 0.0, -0.05],
                                spc_mat=spc_mat, resolution=reso, solver_=solver, out_dir=out_dir, makeplot=true, trueE=Es)

env, solvef! = PDBenchmark.stage!(test)

solvef!(env)


# using JLD
# using Plots
# using Statistics

# damage1 = Float64[]
# damage2 = Float64[]
# time = Float64[]
# for i in 0:10:4000
#     data = load("./output/ImpactDiskBB_vv/env_1_step_$(i).jld")
#     type_ = data["type"]
#     damage_ = data["damage"]
#     mask1 = type_ .== 1
#     mask2 = type_ .== 2

#     push!(damage1, mean(data["damage"][mask1]))
#     push!(damage2, mean(data["damage"][mask2]))
#     push!(time, i)
# end

# fig = plot(time, damage1; label="Disk")
# plot!(time, damage2; label="Projectile")
# xlabel!("Time step")
# ylabel!("Average damage")
# savefig("./output/damage_disk.png")

# using PeriDyn
# PeriDyn.jld2ovito("./output/$(out_dir)/env_1_step_*.jld", st_st[solver][:Steps]; start=0, step=Int64(st_st[solver][:Steps]/100))













