module PeriDyn

using CUDA
using Base: Bool
using Folds
using LinearAlgebra
using Dates
using StaticArrays
using ProgressBars
using JLD2
using Flux
using Zygote
using Optim
using Term
using TimerOutputs
using Unitful

VERSION = "0.1.0"
version() = println(@yellow(@bold "PeriDyn v$VERSION"))

# include all files
include("./misc/logo.jl") # print logo
include("./misc/units.jl") # units
include("./misc/log.jl") # logging
include("./misc/const_macros.jl") # constatnts and macros
include("./misc/repr.jl") # printing and displaying variables
include("./misc/utils.jl") # utility functions
include("./misc/types.jl") # utility functions
include("./peri.jl") # peridynamics functions

# define simulation environment
include("./environment.jl")

# material models
include("./materials/material.jl")

# contact models
include("./contacts/contacts.jl")

# boundary condition functions
include("./boundary_conditions/boundary_conditions.jl")

# All implemented solvers
include("./solvers/solvers.jl")

# All implemented save and load functions
include("./io/io.jl")

# All implemented cuda functions
include("./misc/cuda.jl")

end
