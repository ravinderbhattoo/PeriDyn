module PeriDyn

using Base: Bool
using LinearAlgebra

include("./simulation.jl") # define simulation environment
include("./operators.jl") # standars peridynamic operators (redundant)

#material models
include("./materials/material.jl")
include("./materials/bond_based.jl")
include("./materials/ordinary_state_based.jl")
include("./materials/skip_specific.jl")

# contact models
include("./contacts/contacts.jl")

include("./util.jl") # utility functions
include("./boundary_conditions.jl") # boundary condition functions
include("./solvers.jl") # All implemented solvers

end
