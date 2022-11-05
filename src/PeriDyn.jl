module PeriDyn

using Base: Bool
using LinearAlgebra
using Dates
using Flux
using Zygote
using StaticArrays
using Folds
using PDMesh
using ProgressBars
using JLD

include("./macros.jl") # macros

SOLVERS =  Dict()
include("./simulation.jl") # define simulation environment

#material models
include("./materials/material.jl")

include("./util.jl") # utility functions

include("./materials/bond_based.jl")
include("./materials/ordinary_state_based.jl")
include("./materials/skip_specific.jl")
include("./materials/pairwiseNN.jl")
include("./materials/EPS.jl")

# contact models
include("./contacts/contacts.jl")

include("./boundary_conditions/boundary_conditions.jl") # boundary condition functions
include("./solvers/solvers.jl") # All implemented solvers
include("./io/io.jl") 

end
