using Pkg
Pkg.activate("./test")
Pkg.instantiate()

using Revise
using PeriDyn
using PDMesh

# new block
DIR = "."
device = :cpu

PeriDyn.set_device(device)
PDMesh.set_device(device)

# new block
# Define a block of material and its properties
resolution = 0.2
L = 2
Es = 70 # GPa
nu = 0.15
K = Es/3/(1-2nu)
G = Es/2/(1+nu)
den = 2.2 * 1000.0 # Kg/m3
cstretch = 0.02

out1 = create(Cuboid([-2L 2L; -2L 2L; -L L]); resolution=resolution)  # mm
mat_gen1 = GeneralMaterial(out1, 3.0*resolution; max_neigh=200)
# mat_spec1 = BondBasedSpecific([K], [cstretch], [den])
mat_spec1 = OrdinaryStateBasedSpecific([K], [G], [cstretch], [den])
# mat_spec1 = SkipSpecific()
specimen = PeridynamicsMaterial(mat_gen1, mat_spec1; name="Specimen")

# new block
# Define Indentor Conical as of now
height = 2
radius = 2
Es_indentor = 100*Es # GPa
K_indentor = Es_indentor/3/(1-2nu)
G_indentor = Es_indentor/2/(1+nu)

cone = Cone(radius, height)
obj = move(rotate(cone, angle=180), by=[0.0, 0.0, L + height + 3*resolution])
out2 = create(obj; resolution=resolution, type=2)  # mm
mat_gen2 = GeneralMaterial(out2, 3.0*resolution; max_neigh=200)
# mat_spec2 = BondBasedSpecific([K_indentor], [cstretch], [den])
mat_spec2 = OrdinaryStateBasedSpecific([K_indentor], [G_indentor], [cstretch], [den])
indentor = PeridynamicsMaterial(mat_gen2, mat_spec2; name="Indentor")

# new block
solver = DSVelocityVerlet()
Steps, fwf, nuf = 1000, 100, 10
dt = 1.0

# new block
RMs = []
k = 1.0e6
RM1 = LinearRepulsionModel(k, specimen; distanceX=3, max_neighs=200)
RM2 = LinearRepulsionModel(k, indentor; distanceX=3, max_neighs=200)
RM3 = LinearRepulsionModel(k, indentor, specimen; distanceX=3, max_neighs=200)
push!(RMs, RM1, RM2, RM3)

env = PeriDyn.Env(1, [specimen, indentor], RMs, [], dt)

BC1 = FixBC(env.y[3, :] .< -L + 3*resolution)
BC2 = MoveBC(env.y[3, :] .> L + 2*resolution, [0.0, 0.0, -5*resolution*height/Steps/dt])
BCs = [BC1, BC2]

env.boundary_conditions = BCs

out_dir="$(DIR)/output/indentor_sim_OSB/$(typeof(solver))"

mkpath("$(out_dir)/indentor_specimen.data")
save_state_ovito_bc!("$(out_dir)/indentor_specimen.data", env)
run!([env], Steps, solver;
    filewrite_freq=fwf, neigh_update_freq=nuf, out_dir=out_dir, start_at=0, ext=:data)
