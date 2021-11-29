module PeriDyn

using Base: Bool
using LinearAlgebra
using Dates
using Flux
using Zygote: Buffer
using Distributed, ParallelUtilities

VERBOSE = false

include("./simulation.jl") # define simulation environment
include("./operators.jl") # standars peridynamic operators (redundant)

#material models
include("./materials/material.jl")

include("./util.jl") # utility functions

include("./materials/bond_based.jl")
include("./materials/ordinary_state_based.jl")
include("./materials/skip_specific.jl")
include("./materials/pairwiseNN.jl")

# contact models
include("./contacts/contacts.jl")

include("./boundary_conditions.jl") # boundary condition functions
include("./solvers.jl") # All implemented solvers

end
