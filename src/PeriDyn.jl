module PeriDyn

include("./mesh.jl") # utility fuctions to create mesh
include("./simulation.jl") # define simulation environment
include("./operators.jl") # standars peridynamic operators (redundant)

#material models
include("./materials/material.jl")
include("./materials/bond_based.jl")
include("./materials/ordinary_state_based.jl")

# contact models
include("./contacts/contacts.jl")

include("./util.jl") # utility functions
include("./global_bc.jl") # boundary condition functions
include("./solvers.jl") # All implemented solvers

include("./testtools.jl") # only for testing

end
