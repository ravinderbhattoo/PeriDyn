################################################
# Import packages
################################################
using Revise
using PeriDyn
using PDMaterialPoints
using Unitful
using Plots

################################################
# Set simulation paramters
################################################
# simulation configuration
DIR = joinpath(homedir(), "Downloads", "PeriDyn")

device = :cpu
PeriDyn.set_device(device)
PDMaterialPoints.set_device(device)

################################################
# Define a material block and its properties [SPECIMEN]
################################################
resolution = 1.0u"μm"
δ = 3 * resolution 
L = 10u"μm"
Es = 700u"GPa"
nu = 0.15
K = Es/3/(1-2nu)
G = Es/2/(1+nu)
den = 2.2 * u"g/cm^3"
cstretch = 0.1
σy = Es * 0.1 / 2

# Define Specimen
out1 = create(Cuboid([-2L 2L; -2L 2L; -L L]); resolution=resolution)  
mat_gen1 = GeneralMaterial(out1, δ; max_neigh=300)

mat_spec1 = BondBasedSpecific([K], [cstretch], [den])
# mat_spec1 = OrdinaryStateBasedSpecific([K], [G], [cstretch], [den])
# mat_spec1 = ElastoPlasticSolidSpecific([K], [G], [cstretch], [den], [σy])
# mat_spec1 = SkipSpecific()

specimen = PeridynamicsMaterial(mat_gen1, mat_spec1; name="Specimen")

################################################
# Define a material block and its properties [Indenter]
################################################
height = 5.0u"μm" 
radius = 5.0u"μm"
Es_indentor = 100*Es
K_indentor = Es_indentor/3/(1-2nu)
G_indentor = Es_indentor/2/(1+nu)

# cone = Cone(radius, height)
# obj = move(PDMaterialPoints.rotate(cone; angle=180.0), by=[0.0u"μm", 0.0u"μm", L + height + resolution])

sphere = Sphere(radius)
obj = move(sphere; by=[0.0u"μm", 0.0u"μm", L + radius + 2*resolution])
out2 = create(obj; resolution=resolution, type=2) 
mat_gen2 = GeneralMaterial(out2, δ; max_neigh=200)
mat_spec2 = BondBasedSpecific([K_indentor], [cstretch], [den])
indentor = PeridynamicsMaterial(mat_gen2, mat_spec2; name="Indentor")

################################################
# Define solver and simulation parameters
################################################
solver = DSVelocityVerlet()
steps, fwf, nuf, apf = 10000, 100, 10, 100
dt = 0.1u"μs" / steps

################################################
# Define contact models
################################################
RMs = []
k = 18*K/(π*δ^4)
RM1 = LinearSpringContactModel(k, specimen; distanceX=3, max_neighs=200)
RM2 = LinearSpringContactModel(k, indentor; distanceX=3, max_neighs=200)
RM3 = LinearSpringContactModel(k, indentor, specimen; distanceD=3.0, distanceX=3, max_neighs=200)
push!(RMs, RM1, RM2, RM3)


################################################
# Create environment
################################################
env = PeriDyn.Env(1, [specimen, indentor], RMs, [], dt; units=true)

# Define boundary conditions
wv = sqrt(Es/den)
println("Wave velocity in specimen: ", wv + 0u"m/s")
speed = 10*resolution / (steps*dt)
println("Speed: ", uconvert(u"m/s", speed))

################################################
# Define boundary conditions
################################################
BC1 = FixBC(env.y[3, :] .< (-L + 3*resolution))
BC2 = ToFroBC(env.y[3, :] .> (L + resolution), [0.0, 0.0, -1]*speed, Int(steps*0.8))
BCs = [BC1, BC2]
env.boundary_conditions = BCs

################################################
# Define default units and strip units from the environment
################################################
PeriDyn.set_default_units(; ulength=u"mm", utime=u"μs", umass=u"mg")
env = PeriDyn.ustrip_to_default(env)

function inunits(x; default=false)
    if default
        x = PeriDyn.DEFAULT_UNITS[dimension(x)]
    end
    val = Float64(uconvert(x, 1*PeriDyn.DEFAULT_UNITS[dimension(x)]))
    if ustrip(val) == 1.0
        return unit(val)
    else
        return val
    end
end

################################################
# Make plot during simulation
################################################
function make_plot(Out)
    x = Out[:Displacement]
    y = Out[:Force][:, 3]
    fig = plot(x, y,
        xlabel="Displacement [$(inunits(u"m"))]",
        ylabel="Force [$(inunits(u"N"))]",
        label="$(typeof(mat_spec1))",
        legend=:topleft)
    savefig(fig, "$(out_dir)/force_displacement.png")
end

indentor_ = (env.type .== 2);
indentor_max = maximum(env.y[3, indentor_])
env.Params = Dict(:indentor => indentor_, :indentor_max=>indentor_max)
env.Out = Dict(
        :Force => zeros(eltype(first(env.f)*first(env.mass)), steps, size(env.f, 1)),
        :Displacement => zeros(eltype(env.y), steps)
    )
env.Collect! = (env, step) -> begin
        mask = env.Params[:indentor]
        env.Out[:Force][step, :] =
            sum(env.f[:, mask]' .* env.mass[mask], dims=1)
        env.Out[:Displacement][step] =
            -(
                maximum(env.y[3, mask])
                - env.Params[:indentor_max]
            )

        if step%100 == 0
            make_plot(env.Out)
        end
end

################################################
# Define output directory
################################################
out_dir="$(DIR)/output/indentor_sim_OSB/$(typeof(solver))"
rm(out_dir; force=true, recursive=true)
mkpath("$(out_dir)")
save_state_ovito_bc!("$(out_dir)/indentor_specimen.data", env)

################################################
# Run simulation
################################################
run!(env, steps, solver;
    filewrite_freq=fwf,
    average_prop_freq=apf,
    neigh_update_freq=nuf,
    out_dir=out_dir,
    start_at=0,
    ext=:data
    )


################################################
# Time benchmark
################################################
# using TimerOutputs
# function dotime(env)
#     enable_timer!(PeriDyn.timings);
#     apply_solver!(env, DSVelocityVerlet());
#     show(PeriDyn.timings)
# end