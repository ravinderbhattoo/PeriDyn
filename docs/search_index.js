var documenterSearchIndex = {"docs":
[{"location":"mesh.html#Peridynamics-mesh-and-default-shapes","page":"Peridynamics mesh and default shapes","title":"Peridynamics mesh and default shapes","text":"","category":"section"},{"location":"mesh.html","page":"Peridynamics mesh and default shapes","title":"Peridynamics mesh and default shapes","text":"PDMesh.create\r\nPDMesh.delete\r\nPDMesh.move","category":"page"},{"location":"mesh.html#PDMesh.create","page":"Peridynamics mesh and default shapes","title":"PDMesh.create","text":"create(shape::T; resolution=nothing, rand_=0.0, type::Int64=1) where T <: Shape\n\nAbstact function for creating Shape objects.\n\nReturns\n\n- X : Initial reference position\n- V : Initial velocity\n- Y : Initial position\n- volume : Volume per particle point\n- type: Type of particle point\n\n\n\n\n\n","category":"function"},{"location":"mesh.html#PDMesh.delete","page":"Peridynamics mesh and default shapes","title":"PDMesh.delete","text":"delete(obj::T, f::Function) where T\n\nDelete mesh particle for object using boolean array from function f.\n\n\n\n\n\n","category":"function"},{"location":"operatorandutil.html#Operators-and-utilities","page":"Operator and Utility","title":"Operators and utilities","text":"","category":"section"},{"location":"autodocs.html","page":"Autodocs","title":"Autodocs","text":"Modules = [PeriDyn, PDMesh]","category":"page"},{"location":"autodocs.html#PeriDyn.BondBasedMaterial","page":"Autodocs","title":"PeriDyn.BondBasedMaterial","text":"Bond based material type.\n\n\n\n\n\n","category":"type"},{"location":"autodocs.html#PeriDyn.BondBasedSpecific","page":"Autodocs","title":"PeriDyn.BondBasedSpecific","text":"Specific bond based material type.\n\n\n\n\n\n","category":"type"},{"location":"autodocs.html#PeriDyn.BoundaryCondition","page":"Autodocs","title":"PeriDyn.BoundaryCondition","text":"Abstract class for boundary conditions.\n\n\n\n\n\n","category":"type"},{"location":"autodocs.html#PeriDyn.ElastoPlasticSolidMaterial","page":"Autodocs","title":"PeriDyn.ElastoPlasticSolidMaterial","text":"ElastoPlasticSolidMaterial(type::UnitRange{Int64}, general::GeneralMaterial, specific::ElastoPlasticSolidSpecific)\n\n\n\n\n\n","category":"type"},{"location":"autodocs.html#PeriDyn.ElastoPlasticSolidSpecific","page":"Autodocs","title":"PeriDyn.ElastoPlasticSolidSpecific","text":"Specific Elasto Plastic Solid Material type.\n\n\n\n\n\n","category":"type"},{"location":"autodocs.html#PeriDyn.ElastoPlasticSolidSpecific-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Any}","page":"Autodocs","title":"PeriDyn.ElastoPlasticSolidSpecific","text":"ElastoPlasticSolidSpecific(bulk_modulus::Array{Float64, 1}, shear_modulus::Array{Float64,1}, critical_stretch::Array{Float64,1}, density::Array{Float64,1})\n\nSpecific Elasto Plastic Solid Material type.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.ElastoPlasticSolidSpecificFull","page":"Autodocs","title":"PeriDyn.ElastoPlasticSolidSpecificFull","text":"Specific Elasto Plastic Solid Material type.\n\n\n\n\n\n","category":"type"},{"location":"autodocs.html#PeriDyn.FixBC-Tuple{Any}","page":"Autodocs","title":"PeriDyn.FixBC","text":"FixBC(bool; onlyatstart=false)\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.GeneralMaterial-NTuple{6, Any}","page":"Autodocs","title":"PeriDyn.GeneralMaterial","text":"GeneralMaterial(y0,v0,x,volume,type,horizon; particle_size=0,max_neigh=50)\n\nGeneral peridynamics material type.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.GeneralMaterial-Tuple{Dict, Vararg{Any}}","page":"Autodocs","title":"PeriDyn.GeneralMaterial","text":"\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.OrdinaryStateBasedMaterial","page":"Autodocs","title":"PeriDyn.OrdinaryStateBasedMaterial","text":"OrdinaryStateBasedMaterial(type::UnitRange{Int64}, general::GeneralMaterial, specific::OrdinaryStateBasedSpecific)\n\n\n\n\n\n","category":"type"},{"location":"autodocs.html#PeriDyn.OrdinaryStateBasedSpecific","page":"Autodocs","title":"PeriDyn.OrdinaryStateBasedSpecific","text":"Specific ordinary state based materail type.\n\n\n\n\n\n","category":"type"},{"location":"autodocs.html#PeriDyn.OrdinaryStateBasedSpecific-NTuple{4, Vector{Float64}}","page":"Autodocs","title":"PeriDyn.OrdinaryStateBasedSpecific","text":"OrdinaryStateBasedSpecific(bulk_modulus::Array{Float64, 1}, shear_modulus::Array{Float64,1}, critical_stretch::Array{Float64,1}, density::Array{Float64,1})\n\nSpecific ordinary state based materail type.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.PairwiseNNMaterial","page":"Autodocs","title":"PeriDyn.PairwiseNNMaterial","text":"Bond based material type (NN).\n\n\n\n\n\n","category":"type"},{"location":"autodocs.html#PeriDyn.PairwiseNNSpecific","page":"Autodocs","title":"PeriDyn.PairwiseNNSpecific","text":"Specific bond based material type (NN) .\n\n\n\n\n\n","category":"type"},{"location":"autodocs.html#PeriDyn.PeridynamicsMaterial","page":"Autodocs","title":"PeriDyn.PeridynamicsMaterial","text":"Abstract PeridynamicsMaterial type.\n\n\n\n\n\n","category":"type"},{"location":"autodocs.html#PeriDyn.PeridynamicsMaterial-NTuple{4, Any}","page":"Autodocs","title":"PeriDyn.PeridynamicsMaterial","text":"PeridynamicsMaterial(gen, spc)\n\nCreate peridynamics material.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.PeridynamicsMaterial-Tuple{Any, Any, Any, ElastoPlasticSolidSpecific}","page":"Autodocs","title":"PeriDyn.PeridynamicsMaterial","text":"PeridynamicsMaterial(gen, spc::ElastoPlasticSolidSpecific)\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.PeridynamicsMaterial-Tuple{Any, Any, Any, SkipSpecific}","page":"Autodocs","title":"PeriDyn.PeridynamicsMaterial","text":"PeridynamicsMaterial(gen, spc::SkipSpecific)\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.PeridynamicsMaterial-Tuple{Any, Any}","page":"Autodocs","title":"PeriDyn.PeridynamicsMaterial","text":"PeridynamicsMaterial(gen, spc::PeridynamicsMaterial; name=\"PM\")\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.SkipMaterial","page":"Autodocs","title":"PeriDyn.SkipMaterial","text":"SkipMaterial(type::UnitRange{Int64}, general::GeneralMaterial, specific::SkipSpecific)\n\n\n\n\n\n","category":"type"},{"location":"autodocs.html#PeriDyn.SkipSpecific","page":"Autodocs","title":"PeriDyn.SkipSpecific","text":"Specific skip materail type.\n\n\n\n\n\n","category":"type"},{"location":"autodocs.html#PeriDyn.SkipSpecific-Tuple{}","page":"Autodocs","title":"PeriDyn.SkipSpecific","text":"Specific skip materail type.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.Env-Tuple{Int64, Any, Any, Any, Any}","page":"Autodocs","title":"PeriDyn.Env","text":"Env(id::Int64,materials,short_range_repulsion,boundary_conds,dt;state=2)\n\nCreate a GeneralEnv for holding parameters for a simulation.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.LJRepulsionModel-Tuple{Any, Any, PeridynamicsMaterial, PeridynamicsMaterial}","page":"Autodocs","title":"PeriDyn.LJRepulsionModel","text":"LJRepulsionModel(alpha::Float64,epsilon::Float64, mat1::PeridynamicsMaterial,mat2::PeridynamicsMaterial;distanceX=5,max_neighs=50)\n\nLJ repulsive model for 1-2 material blocks.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.LJRepulsionModel-Tuple{Any, Any, PeridynamicsMaterial}","page":"Autodocs","title":"PeriDyn.LJRepulsionModel","text":"LJRepulsionModel(alpha::Float64,epsilon::Float64, mat1::PeridynamicsMaterial;distanceX=5,max_neighs=50)\n\nLJ repulsive model for 1-1 material blocks.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.LinearRepulsionModel-Tuple{Any, PeridynamicsMaterial, PeridynamicsMaterial}","page":"Autodocs","title":"PeriDyn.LinearRepulsionModel","text":"LinearRepulsionModel(stifness::Float64, mat1::PeridynamicsMaterial,mat2::PeridynamicsMaterial;distanceX=5,max_neighs=50)\n\nLinear repulsive model for 1-2 material blocks.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.LinearRepulsionModel-Tuple{Any, PeridynamicsMaterial}","page":"Autodocs","title":"PeriDyn.LinearRepulsionModel","text":"LinearRepulsionModel(stifness::Float64, mat1::PeridynamicsMaterial;distanceX=5,max_neighs=50)\n\nLinear repulsive model for 1-1 material blocks.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.NonLinearRepulsionModel-Tuple{Any, Any, PeridynamicsMaterial, PeridynamicsMaterial}","page":"Autodocs","title":"PeriDyn.NonLinearRepulsionModel","text":"NonLinearRepulsionModel(exponent::Float64,stifness::Float64, mat1::PeridynamicsMaterial,mat2::PeridynamicsMaterial;distanceX=5,max_neighs=50)\n\nNonLinear repulsive model for 1-2 material blocks.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.NonLinearRepulsionModel-Tuple{Any, Any, PeridynamicsMaterial}","page":"Autodocs","title":"PeriDyn.NonLinearRepulsionModel","text":"NonLinearRepulsionModel(exponent::Float64,stifness::Float64, mat1::PeridynamicsMaterial;distanceX=5,max_neighs=50)\n\nNonLinear repulsive model for 1-1 material blocks.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.ShortRangeRepulsionModel-Tuple{Any, PeridynamicsMaterial, PeridynamicsMaterial}","page":"Autodocs","title":"PeriDyn.ShortRangeRepulsionModel","text":"ShortRangeRepulsionModel(spring_const::Float64, mat1::PeridynamicsMaterial,mat2::PeridynamicsMaterial;distanceX=5,max_neighs=50)\n\nShortRange repulsive model for 1-2 material blocks.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.ShortRangeRepulsionModel-Tuple{Any, PeridynamicsMaterial}","page":"Autodocs","title":"PeriDyn.ShortRangeRepulsionModel","text":"ShortRangeRepulsionModel(spring_const::Float64, mat1::PeridynamicsMaterial;distanceX=5,max_neighs=50)\n\nShortRange repulsive model for 1-1 material blocks.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn._tij-NTuple{6, Any}","page":"Autodocs","title":"PeriDyn._tij","text":"_tij(x,mi,thetai,eij,K,G)::Float64\n\ncalculates tij force density magnitude (actually acceleration).\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.apply_bc!-Union{Tuple{T}, Tuple{Any, T, Type{Val{:position}}}} where T<:PeriDyn.BoundaryCondition","page":"Autodocs","title":"PeriDyn.apply_bc!","text":"apply_bc!(env,BC::BoundaryCondition, :position)\n\nApplies general boundary condition to a given material block.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.apply_bc!-Union{Tuple{T}, Tuple{Any, T, Type{Val{:velocity}}}} where T<:PeriDyn.BoundaryCondition","page":"Autodocs","title":"PeriDyn.apply_bc!","text":"apply_bc!(env,BC::BoundaryCondition, :velocity)\n\nApplies general boundary condition to a given material block.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.cal_family!-Tuple{Matrix{Int64}, Matrix{Float64}, Float64}","page":"Autodocs","title":"PeriDyn.cal_family!","text":"cal_family!(family::Array{Int64,2},x::Array{Float64,2}, horizon::Float64)\n\nCalculate family members for each material point (inplace).\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.cal_family-Tuple{Matrix{Float64}, Float64, Int64}","page":"Autodocs","title":"PeriDyn.cal_family","text":"cal_family!(family::Array{Int64,2},x::Array{Float64,2}, horizon::Float64)\n\nCalculate family members for each material point.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.cell_number-NTuple{4, Any}","page":"Autodocs","title":"PeriDyn.cell_number","text":"cell_number(i,j,k,N)\n\nCalculates cell number for given i,j,k cell indices (when in a list).\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.collision_box-Tuple{Matrix{Float64}, Matrix{Float64}, Float64}","page":"Autodocs","title":"PeriDyn.collision_box","text":"collision_box(x1::Array{Float64,2}, x2::Array{Float64,2}, skin::Float64)\n\nCalculates collsion box between two material blocks.\n\nInput Args:\n\nx1 :: Positions of material point (block 1)\nx2 :: Positions of material point (block 2)\nskin :: Extra distance need to consider (usually >= particle size)\n\nOutput Args:\n\nbox_min :: Minimum position limits for overlap\nbox_max :: Miximum position limits for overlap\nifoverlap :: Boolean (true if overlap)\n\nExamples\n\njulia> collision_box(x1, x2, skin)\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.dilatation-Tuple{Matrix{Float64}, GeneralMaterial, Vector{Float64}}","page":"Autodocs","title":"PeriDyn.dilatation","text":"dilatation(y::Array{Float64,2},S::GeneralMaterial,m::Array{Float64,1})\n\nIt gives dilatation as given ordinary state material model.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.dilatation_theta-Tuple{Matrix{Float64}, GeneralMaterial}","page":"Autodocs","title":"PeriDyn.dilatation_theta","text":"dilatation_theta(y::Array{Float64,2},S::GeneralMaterial)\n\nIt gives dilatation as given ordinary state material model.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.force_density_T-Tuple{Matrix{Float64}, BondBasedMaterial}","page":"Autodocs","title":"PeriDyn.force_density_T","text":"force_density_T(mat::BondBasedMaterial)\n\nCalculates force density (actually acceleration) for bond based material type.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.force_density_T-Tuple{Matrix{Float64}, ElastoPlasticSolidMaterial}","page":"Autodocs","title":"PeriDyn.force_density_T","text":"force_density_T(y::Array{Float64,2},mat::ElastoPlasticSolidMaterial)\n\nCalculates force density (actually acceleration) for ordinary state based material type.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.force_density_T-Tuple{Matrix{Float64}, OrdinaryStateBasedMaterial}","page":"Autodocs","title":"PeriDyn.force_density_T","text":"force_density_T(y::Array{Float64,2},mat::OrdinaryStateBasedMaterial)\n\nCalculates force density (actually acceleration) for ordinary state based material type.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.force_density_T-Tuple{Matrix{Float64}, PairwiseNNMaterial}","page":"Autodocs","title":"PeriDyn.force_density_T","text":"force_density_T(mat::PairwiseNNMaterial)\n\nCalculates force density (actually acceleration) for bond based material type.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.force_density_T-Tuple{Matrix{Float64}, SkipMaterial}","page":"Autodocs","title":"PeriDyn.force_density_T","text":"force_density_T(y::Array{Float64,2},mat::SkipMaterial)\n\nCalculates force density (actually acceleration) for skip material type.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.get_cells-Tuple{Matrix{Float64}, Number}","page":"Autodocs","title":"PeriDyn.get_cells","text":"get_cells(x::Array{Float64,2},horizon::Float64)\n\nFill cells with material points.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.horizon_correction-Tuple{Any, Any, Any}","page":"Autodocs","title":"PeriDyn.horizon_correction","text":"horizon_correction(dr,r,hr)\n\nIt gives horizon correction factor (It will give 1 as of now).\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.influence_function-Tuple{Any}","page":"Autodocs","title":"PeriDyn.influence_function","text":"influence_function(dr)\n\nIt gives influence function factor (It will give 1/r as of now).\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.make_NN-Union{Tuple{T}, Tuple{T, Any}} where T","page":"Autodocs","title":"PeriDyn.make_NN","text":"make_NN(layers::Tuple{T}, N) where T\n\nCreate an symmetrical NxN matrix from a vector of length N(N+1)/2.    \n\n================\n\nReturns\n\nM :: Matrix\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.make_matrix-Union{Tuple{Vector{T}}, Tuple{T}} where T","page":"Autodocs","title":"PeriDyn.make_matrix","text":"make_matrix(S::Array{T,1})\n\nCreate an symmetrical NxN matrix from a vector of length N(N+1)/2.    \n\n================\n\nReturns\n\nM :: Matrix\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.make_matrix_gm-Union{Tuple{Vector{T}}, Tuple{T}} where T","page":"Autodocs","title":"PeriDyn.make_matrix_gm","text":"make_matrix(S::Array{T,1})\n\nCreate an symmetrical NxN matrix from a vector of length N(N+1)/2.    \n\n================\n\nReturns\n\nM :: Matrix\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.neigh_cells-NTuple{4, Any}","page":"Autodocs","title":"PeriDyn.neigh_cells","text":"neigh_cells(i,j,k,N)\n\nCalculate all neighboring cells for i,j,k cell.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.plastic_yield_function-Tuple{Any, Any, Any}","page":"Autodocs","title":"PeriDyn.plastic_yield_function","text":"yield function\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.print_data_file!-Tuple{Array{PeriDyn.GeneralEnv}, String, Int64}","page":"Autodocs","title":"PeriDyn.print_data_file!","text":"print_data_file!(envs::Array{GeneralEnv}, file_prefix::String, i::Int64)\n\nWrites data file to disk.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.quasi_static!-Tuple{Any, Int64, Float64}","page":"Autodocs","title":"PeriDyn.quasi_static!","text":"quasi_static!(envs::Any,N::Int64,step_size::Float64; max_iter::Int64=100, min_step_tol_per::Float64=0.5, filewrite_freq::Int64=10, neigh_update_freq::Int64=50, file_prefix::String=\"datafile\",start_at::Int64=0,write_from::Int=0)\n\nImplements quasi static simulation using minimize for each step.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.repulsion_force-Tuple{Any, LJRepulsionModel11}","page":"Autodocs","title":"PeriDyn.repulsion_force","text":"repulsion_force(dr, RepMod::LJRepulsionModel11)\n\nCalculates repulsive acceleration for 1-1 materials block interaction.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.repulsion_force-Tuple{Any, LJRepulsionModel12}","page":"Autodocs","title":"PeriDyn.repulsion_force","text":"repulsion_force(dr, RepMod::LJRepulsionModel12)\n\nCalculates repulsive acceleration for 1-2 materials block interaction.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.repulsion_force-Tuple{Any, NonLinearRepulsionModel11}","page":"Autodocs","title":"PeriDyn.repulsion_force","text":"repulsion_force(dr,RepMod::NonLinearRepulsionModel11)\n\nCalculates repulsive acceleration for 1-1 materials block interaction.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.repulsion_force-Tuple{Any, NonLinearRepulsionModel12}","page":"Autodocs","title":"PeriDyn.repulsion_force","text":"repulsion_force(dr,RepMod::NonLinearRepulsionModel12)\n\nCalculates repulsive acceleration for 1-2 materials block interaction.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.repulsion_force-Tuple{Any, ShortRangeRepulsionModel11}","page":"Autodocs","title":"PeriDyn.repulsion_force","text":"repulsion_force(dr, RepMod::ShortRangeRepulsionModel11)\n\nCalculates repulsive acceleration for 1-1 materials block interaction.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.repulsion_force-Tuple{Any, ShortRangeRepulsionModel12}","page":"Autodocs","title":"PeriDyn.repulsion_force","text":"repulsion_force(dr, RepMod::ShortRangeRepulsionModel12)\n\nCalculates repulsive acceleration for 1-2 materials block interaction.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.save_state!-Tuple{Any, Any}","page":"Autodocs","title":"PeriDyn.save_state!","text":"save_state!(filename, env)\n\nSave env::GeneralEnv to disk.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.set_env_active!-Tuple{Any}","page":"Autodocs","title":"PeriDyn.set_env_active!","text":"set_env_active!(env)\n\nSet environment state active.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.set_env_idel!-Tuple{Any}","page":"Autodocs","title":"PeriDyn.set_env_idel!","text":"set_env_idel!(env)\n\nSet environment state idel.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.set_env_inactive!-Tuple{Any}","page":"Autodocs","title":"PeriDyn.set_env_inactive!","text":"set_env_inactive!(env)\n\nSet environment state inactive.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.set_ghost_atoms!-Tuple{Any, Any}","page":"Autodocs","title":"PeriDyn.set_ghost_atoms!","text":"set_ghost_atoms!(env,ghost_atoms)\n\nSet ghost atoms for an environment.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.short_range_repulsion!-Tuple{Any, Any, Any, Any, PeriDyn.RepulsionModel11}","page":"Autodocs","title":"PeriDyn.short_range_repulsion!","text":"short_range_repulsion!(y,f,type,RepusionModel)\n\nUpdates (inplace) the repulsive acceleration of material points.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.short_range_repulsion!-Tuple{Any, Any, Any, Any, PeriDyn.RepulsionModel12}","page":"Autodocs","title":"PeriDyn.short_range_repulsion!","text":"short_range_repulsion!(y,f,type,RepusionModel)\n\nUpdates (inplace) the repulsive acceleration of material points.\n\nα\n\n1-1 repulsion\n\nInput Args:\n\ny :: Positions of material point\nf :: Acceleration of material points\ntype :: Type of material points\nRepulsionModel :: Repulsion model (see contacts.jl for more details)\n\nOutput Args:\n\nNoting (Inplace updation of f (acceleration))\n\nExamples\n\njulia> short_range_repulsion!(y,f,type,RepusionModel)\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.simulate-Tuple","page":"Autodocs","title":"PeriDyn.simulate","text":"\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.td_norm_fn-Tuple{Any, Any, Any}","page":"Autodocs","title":"PeriDyn.td_norm_fn","text":"td_norm_fn\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.update_acc!-Tuple{PeriDyn.GeneralEnv}","page":"Autodocs","title":"PeriDyn.update_acc!","text":"update_acc!(env::GeneralEnv)\n\nUpdates acceleration of all the material points in a simulation environment.\n\n̈u(xᵢ, t) = from material forces + from contact forces\n\nfrom material force, ̈u(xᵢ, t) = [∑ᵏⱼ₌₁ {T[xᵢ, t]<xⱼ-xᵢ> - T[xⱼ, t]<xᵢ-xⱼ> }*Vⱼ + b(xᵢ, t)] / ρᵢ \n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.update_neighs!-Tuple{Any}","page":"Autodocs","title":"PeriDyn.update_neighs!","text":"function update_neighs!(envs)\n\nUpdates neighbors of each material point for a list of simulation environments.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.update_repulsive_neighs!-Tuple{Any, Any, PeriDyn.RepulsionModel11}","page":"Autodocs","title":"PeriDyn.update_repulsive_neighs!","text":"update_repulsive_neighs!(y,type,RepulsionModel11)\n\nUpdate neighbour list for repulsive force calculation (1-1 interaction).\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.update_repulsive_neighs!-Tuple{Any, Any, PeriDyn.RepulsionModel12}","page":"Autodocs","title":"PeriDyn.update_repulsive_neighs!","text":"update_repulsive_neighs!(y,type,RepulsionModel12)\n\nUpdate neighbour list for repulsive force calculation (1-2 interaction).\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.velocity_verlet!-Tuple{Any, Int64}","page":"Autodocs","title":"PeriDyn.velocity_verlet!","text":"velocityverlet!(envs::Any,N::Int64;filewritefreq=10,neighupdatefreq=50,fileprefix=\"datafile\",startat::Int64=0,write_from::Int=0)\n\nVelocity verlet :).\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.velocity_verlet_step!-Tuple{PeriDyn.GeneralEnv}","page":"Autodocs","title":"PeriDyn.velocity_verlet_step!","text":"velocityverletstep!(env::GeneralEnv)\n\nImplement a single step of velocity verlet algorithm.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.weighted_volume-NTuple{5, Any}","page":"Autodocs","title":"PeriDyn.weighted_volume","text":"weighted_volume(S::GeneralMaterial)\n\nIt gives weighted volume.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.write_data-NTuple{4, Any}","page":"Autodocs","title":"PeriDyn.write_data","text":"write_data(filename,every,type,pos)\n\nIt writes the data file.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.write_data-Tuple{String, Vararg{Any}}","page":"Autodocs","title":"PeriDyn.write_data","text":"write_data(filename::String, args...)\n\nIt writes the data file.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.write_data-Tuple{String, Vector{Int64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}","page":"Autodocs","title":"PeriDyn.write_data","text":"write_data(filename::String,type::Array{Int64,1},y::Array{Float64,2}, v::Array{Float64,2},f::Array{Float64,2},)\n\nIt writes the data file.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.write_data-Tuple{String, Vector{Int64}, Matrix{Float64}}","page":"Autodocs","title":"PeriDyn.write_data","text":"write_data(filename::String,type::Array{Int64,1},y::Array{Float64,2})\n\nIt writes the data file.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PeriDyn.write_data_cell_ids-Tuple{String, Matrix, Any}","page":"Autodocs","title":"PeriDyn.write_data_cell_ids","text":"write_data_cell_ids(filename::String, y::Array{Float64,2}, cells::Any)\n\nIt writes the data file.\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PDMesh.create-Tuple{T} where T<:Shape","page":"Autodocs","title":"PDMesh.create","text":"create(shape::T; resolution=nothing, rand_=0.0, type::Int64=1) where T <: Shape\n\nAbstact function for creating Shape objects.\n\nReturns\n\n- X : Initial reference position\n- V : Initial velocity\n- Y : Initial position\n- volume : Volume per particle point\n- type: Type of particle point\n\n\n\n\n\n","category":"method"},{"location":"autodocs.html#PDMesh.delete-Union{Tuple{T}, Tuple{T, Function}} where T","page":"Autodocs","title":"PDMesh.delete","text":"delete(obj::T, f::Function) where T\n\nDelete mesh particle for object using boolean array from function f.\n\n\n\n\n\n","category":"method"},{"location":"cmodels.html#Contact-models","page":"Contact models","title":"Contact models","text":"","category":"section"},{"location":"cmodels.html","page":"Contact models","title":"Contact models","text":"RepulsionModel11\r\nRepulsionModel11_gf\r\nRepulsionModel11_gcal\r\nRepulsionModel12\r\nRepulsionModel12_gf\r\nRepulsionModel12_gcal\r\n\r\nLJRepulsionModel\r\nLJRepulsionModel11\r\nLJRepulsionModel12\r\n\r\nShortRangeRepulsionModel\r\nShortRangeRepulsionModel11 \r\nShortRangeRepulsionModel12\r\n\r\nNonLinearRepulsionModel\r\nNonLinearRepulsionModel11\r\nNonLinearRepulsionModel12\r\n\r\nLinearRepulsionModel\r\n\r\nrepulsion_force\r\n\r\nshort_range_repulsion! \r\ncollision_box\r\nupdate_repulsive_neighs!","category":"page"},{"location":"cmodels.html#PeriDyn.LJRepulsionModel","page":"Contact models","title":"PeriDyn.LJRepulsionModel","text":"LJRepulsionModel(alpha::Float64,epsilon::Float64, mat1::PeridynamicsMaterial,mat2::PeridynamicsMaterial;distanceX=5,max_neighs=50)\n\nLJ repulsive model for 1-2 material blocks.\n\n\n\n\n\nLJRepulsionModel(alpha::Float64,epsilon::Float64, mat1::PeridynamicsMaterial;distanceX=5,max_neighs=50)\n\nLJ repulsive model for 1-1 material blocks.\n\n\n\n\n\n","category":"function"},{"location":"cmodels.html#PeriDyn.ShortRangeRepulsionModel","page":"Contact models","title":"PeriDyn.ShortRangeRepulsionModel","text":"ShortRangeRepulsionModel(spring_const::Float64, mat1::PeridynamicsMaterial,mat2::PeridynamicsMaterial;distanceX=5,max_neighs=50)\n\nShortRange repulsive model for 1-2 material blocks.\n\n\n\n\n\nShortRangeRepulsionModel(spring_const::Float64, mat1::PeridynamicsMaterial;distanceX=5,max_neighs=50)\n\nShortRange repulsive model for 1-1 material blocks.\n\n\n\n\n\n","category":"function"},{"location":"cmodels.html#PeriDyn.NonLinearRepulsionModel","page":"Contact models","title":"PeriDyn.NonLinearRepulsionModel","text":"NonLinearRepulsionModel(exponent::Float64,stifness::Float64, mat1::PeridynamicsMaterial,mat2::PeridynamicsMaterial;distanceX=5,max_neighs=50)\n\nNonLinear repulsive model for 1-2 material blocks.\n\n\n\n\n\nNonLinearRepulsionModel(exponent::Float64,stifness::Float64, mat1::PeridynamicsMaterial;distanceX=5,max_neighs=50)\n\nNonLinear repulsive model for 1-1 material blocks.\n\n\n\n\n\n","category":"function"},{"location":"cmodels.html#PeriDyn.LinearRepulsionModel","page":"Contact models","title":"PeriDyn.LinearRepulsionModel","text":"LinearRepulsionModel(stifness::Float64, mat1::PeridynamicsMaterial,mat2::PeridynamicsMaterial;distanceX=5,max_neighs=50)\n\nLinear repulsive model for 1-2 material blocks.\n\n\n\n\n\nLinearRepulsionModel(stifness::Float64, mat1::PeridynamicsMaterial;distanceX=5,max_neighs=50)\n\nLinear repulsive model for 1-1 material blocks.\n\n\n\n\n\n","category":"function"},{"location":"cmodels.html#PeriDyn.repulsion_force","page":"Contact models","title":"PeriDyn.repulsion_force","text":"repulsion_force(dr, RepMod::LJRepulsionModel12)\n\nCalculates repulsive acceleration for 1-2 materials block interaction.\n\n\n\n\n\nrepulsion_force(dr, RepMod::LJRepulsionModel11)\n\nCalculates repulsive acceleration for 1-1 materials block interaction.\n\n\n\n\n\nrepulsion_force(dr, RepMod::ShortRangeRepulsionModel12)\n\nCalculates repulsive acceleration for 1-2 materials block interaction.\n\n\n\n\n\nrepulsion_force(dr, RepMod::ShortRangeRepulsionModel11)\n\nCalculates repulsive acceleration for 1-1 materials block interaction.\n\n\n\n\n\nrepulsion_force(dr,RepMod::NonLinearRepulsionModel12)\n\nCalculates repulsive acceleration for 1-2 materials block interaction.\n\n\n\n\n\nrepulsion_force(dr,RepMod::NonLinearRepulsionModel11)\n\nCalculates repulsive acceleration for 1-1 materials block interaction.\n\n\n\n\n\n","category":"function"},{"location":"cmodels.html#PeriDyn.short_range_repulsion!","page":"Contact models","title":"PeriDyn.short_range_repulsion!","text":"short_range_repulsion!(y,f,type,RepusionModel)\n\nUpdates (inplace) the repulsive acceleration of material points.\n\nα\n\n1-1 repulsion\n\nInput Args:\n\ny :: Positions of material point\nf :: Acceleration of material points\ntype :: Type of material points\nRepulsionModel :: Repulsion model (see contacts.jl for more details)\n\nOutput Args:\n\nNoting (Inplace updation of f (acceleration))\n\nExamples\n\njulia> short_range_repulsion!(y,f,type,RepusionModel)\n\n\n\n\n\nshort_range_repulsion!(y,f,type,RepusionModel)\n\nUpdates (inplace) the repulsive acceleration of material points.\n\n\n\n\n\n","category":"function"},{"location":"cmodels.html#PeriDyn.collision_box","page":"Contact models","title":"PeriDyn.collision_box","text":"collision_box(x1::Array{Float64,2}, x2::Array{Float64,2}, skin::Float64)\n\nCalculates collsion box between two material blocks.\n\nInput Args:\n\nx1 :: Positions of material point (block 1)\nx2 :: Positions of material point (block 2)\nskin :: Extra distance need to consider (usually >= particle size)\n\nOutput Args:\n\nbox_min :: Minimum position limits for overlap\nbox_max :: Miximum position limits for overlap\nifoverlap :: Boolean (true if overlap)\n\nExamples\n\njulia> collision_box(x1, x2, skin)\n\n\n\n\n\n","category":"function"},{"location":"cmodels.html#PeriDyn.update_repulsive_neighs!","page":"Contact models","title":"PeriDyn.update_repulsive_neighs!","text":"update_repulsive_neighs!(y,type,RepulsionModel12)\n\nUpdate neighbour list for repulsive force calculation (1-2 interaction).\n\n\n\n\n\nupdate_repulsive_neighs!(y,type,RepulsionModel11)\n\nUpdate neighbour list for repulsive force calculation (1-1 interaction).\n\n\n\n\n\n","category":"function"},{"location":"toc.html#Table-of-contents","page":"Table of contents","title":"Table of contents","text":"","category":"section"},{"location":"toc.html","page":"Table of contents","title":"Table of contents","text":"Pages = [\r\n            \"peridyn.md\",\r\n            \"mmodels.md\",\r\n            \"cmodels.md\",\r\n            \"solvers.md\",\r\n            \"bc.md\",\r\n            \"operatorandutil.md\",\r\n            \"examples.md\",\r\n            \"mesh.md\",\r\n            \"index.md\", \r\n            \"autodocs.md\"\r\n        ]\r\nDepth = 3","category":"page"},{"location":"zexamples.html#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"zexamples.html#Quasi-static-tensile-test","page":"Examples","title":"Quasi-static tensile test","text":"","category":"section"},{"location":"zexamples.html","page":"Examples","title":"Examples","text":"this is description.","category":"page"},{"location":"zexamples.html","page":"Examples","title":"Examples","text":"a = 1","category":"page"},{"location":"solvers.html#Solvers","page":"Solvers","title":"Solvers","text":"","category":"section"},{"location":"mmodels.html#Material-models","page":"Material models","title":"Material models","text":"","category":"section"},{"location":"bc.html#Boundary-conditions","page":"Boundary Conditions","title":"Boundary conditions","text":"","category":"section"},{"location":"bc.html","page":"Boundary Conditions","title":"Boundary Conditions","text":"FixBC\r\nMoveBC\r\nToFroBC\r\nScaleBC","category":"page"},{"location":"bc.html#PeriDyn.FixBC","page":"Boundary Conditions","title":"PeriDyn.FixBC","text":"FixBC(bool; onlyatstart=false)\n\n\n\n\n\n","category":"type"},{"location":"peridyn.html#PeriDyn.jl","page":"Home","title":"PeriDyn.jl","text":"","category":"section"},{"location":"zzfull.html","page":"-","title":"-","text":"Modules = [PeriDyn]","category":"page"},{"location":"examples.html#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples.html#Basic-examples","page":"Examples","title":"Basic examples","text":"","category":"section"},{"location":"examples.html#Contact-examples","page":"Examples","title":"Contact examples","text":"","category":"section"},{"location":"examples.html#composite-material-examples","page":"Examples","title":"composite-material examples","text":"","category":"section"},{"location":"examples.html#Boundary-condition-examples","page":"Examples","title":"Boundary condition examples","text":"","category":"section"},{"location":"list.html","page":"-","title":"-","text":"CurrentModule = PeriDyn","category":"page"},{"location":"list.html","page":"-","title":"-","text":"","category":"page"},{"location":"index.html#Index","page":"Index","title":"Index","text":"","category":"section"},{"location":"index.html","page":"Index","title":"Index","text":"","category":"page"}]
}
