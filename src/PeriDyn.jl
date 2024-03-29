__precompile__()
module PeriDyn

using CUDA
using Base: Bool
using Folds
using LinearAlgebra

using Dates
using StaticArrays
using ProgressBars
using JLD
using Flux
using Zygote

using Optim

# using Flux
# using Zygote

const PDBlockID = Ref{Int64}(1)
SOLVERS =  Dict()
const DEVICE = Ref{Symbol}(:cpu)

function set_device(device)
    device = get_valid_device(device)
    DEVICE[] = device
    println("PeriDyn: DEVICE set to $(DEVICE[])")
end

function get_valid_device(x)
    out = x
    if x==:cuda
        if !CUDA.functional()
            out = :cpu
            println("PeriDyn: CUDA is not available.")
        end
    else
        out = :cpu
    end
    return out
end

function reset_cuda()
    set_device(:cuda)
end


include("./macros.jl") # macros
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
include("./cuda.jl")

end
