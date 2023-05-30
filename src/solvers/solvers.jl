"""
This module contains solver definitions.
"""

export Solver, QuasiStaticSolver, DynamicSolver
export simulate!, run!, update_acc!, update_neighs!
export apply_bc_at0, check_boundaries!
export filepath_, print_data_file!, save_state!, save_state_ovito_bc!

"""
    Solver

Solver is an abstract type for solvers.
"""
abstract type Solver end

"""
    QuasiStaticSolver <: Solver

QuasiStaticSolver is an abstract type for quasi-static solvers.
"""
abstract type QuasiStaticSolver <: Solver end

"""
    DynamicSolver <: Solver

DynamicSolver is an abstract type for dynamic solvers.
"""
abstract type DynamicSolver <: Solver end

function _cudaconvert(solver::Solver)
    solver
end

@def qssolver_gf begin
    max_iter::Int64
    x_tol::Float64
    f_tol::Float64
end


"""
    simulate!(args...; out_dir="datafile", append_date=true, kwargs...)

Simulate a list of environments. The last argument should be a solver. The `args` and
`kwargs` are passed to `run!` function. The `out_dir` is the directory where the data
files are saved. The `append_date` is a boolean value. If `true`, the date is appended
to the `out_dir` path.
"""
function simulate!(args...; out_dir="datafile", append_date=true, kwargs...)
    solver = last(args)
    log_info("Using solver: $(solver)")
    foldername = filepath_(out_dir; append_date=append_date)
    log_info("Output folder: $foldername")
    return run!(args...; kwargs..., out_dir=foldername)
end

"""
    run!(envs, N::Int64, solver::Solver; filewrite_freq::Int64=10,
        neigh_update_freq::Int64=1, average_prop_freq::Int64=1,
        out_dir::String="datafile",
        start_at::Int64=0, write_from::Int=0, ext::Symbol=:jld,
        max_part=30)

Run simulation for a list of environments.

# Arguments
- `envs`: Array{GeneralEnv}, the list of environments.
- `N`: Int64, the number of steps.
- `solver`: Solver, the solver.

# Keyword Arguments
- `filewrite_freq`: Int64, the frequency of writing data files to disk. Default is 10.
- `neigh_update_freq`: Int64, the frequency of updating neighbors. Default is 1.
- `average_prop_freq`: Int64, the frequency of calculating average properties. Default is 1.
- `out_dir`: String, the directory where the data files are saved. Default is "datafile".
- `start_at`: Int64, the starting step. Default is 0.
- `write_from`: Int, the starting index of the data files. Default is 0.
- `ext`: Symbol, the extension of the data files. Default is :jld.
- `max_part`: Int, the maximum number of particles in a neighborhood. Default is 30.

# See also
- `simulate!`

"""
function run!(envs, N::Int64, solver::Solver; filewrite_freq::Int64=10,
                neigh_update_freq::Int64=1, average_prop_freq::Int64=1,
                out_dir::String="datafile",
                start_at::Int64=0, write_from::Int=0, ext::Symbol=:jld,
                max_part=30)

    mkpath(out_dir)
    _print_data_file!(step) = print_data_file!(envs, out_dir, step; ext=ext)


    # save initial state
    _print_data_file!("start_with")

    update_neighs!(envs; max_part=max_part)

    # apply boundary conditions at t = t0
    for id in 1:size(envs,1)
        if envs[id].state==2
            apply_bc_at0(envs[id], start_at)
        end
    end

    # save initial state
    _print_data_file!(0+write_from)

    # Starting simulation at start point
    N = N + start_at

    for i in (1+start_at):N

        # apply integrator step for all envs
        for id in 1:size(envs,1)
            if envs[id].state==2
                apply_solver!(envs[id], solver)
                if envs[id].Collect! !== nothing
                    envs[id].Collect!(envs[id], i)
                end
            else
                log_impinfo("Env $id is not active. step = $i")
            end
        end

        # calculate average property
        if i%average_prop_freq==0.0 || (i==1 && write_from==0)
            for env in envs
                # println("Momentum: $i", sum(env.p, dims=2))
            end
        end

        # write file to disk
        if i%filewrite_freq==0.0 || (i==1 && write_from==0)
            _print_data_file!(i+write_from)

        end

        # update neighbors
        if i%neigh_update_freq==0.0
            update_neighs!(envs; max_part=max_part)
        end

        percentage = i/N*100
        if percentage%1==0
            log_impinfo("Status: $(round(percentage, digits=3))% completed. step = $i")
        end

        log_detail("Solver iteration : $i")

    end
end


"""
    filepath_(file_prefix::String; append_date=false)

Returns the path to the output folder. If `append_date` is `true`, the date is appended.

# Arguments
- `file_prefix`: String, the prefix of the output folder.

# Keyword Arguments
- `append_date`: Bool, the boolean value to append date to the output folder. Default is
`false`.

"""
function filepath_(file_prefix; append_date=false)
    if append_date
        dtime=replace(string(ceil(Dates.now(), Dates.Second(1))), ":"=>"-")
        return mkpath("./output/"*file_prefix*"_"*dtime)*"/"
    else
        return mkpath("./output/"*file_prefix)*"/"
    end
end


"""
    save_state!(filename, env)

Save `env::GeneralEnv` to disk.

# Arguments
- `filename`: String, the filename.
- `env`: GeneralEnv, the simulation environment.

# See also
- `save_state_ovito_bc!`
"""
function save_state!(filename, env)
    update_acc!(env)
    intact = sum(env.material_blocks[1].general.intact, dims=1)
    for mat in env.material_blocks[2:end]
        intact = hcat(intact, sum(mat.general.intact, dims=1))
    end

    damage = 1 .- reshape(intact, :) ./ env.intact0
    write_data(filename; id=1:length(env.type), type=env.type, position=env.y, velocity=env.v, acceleration=env.f, mass=env.mass, volume=env.volume, damage=damage)
end

"""
    save_state_ovito_bc!(filename, env)

Save `env::GeneralEnv` to disk for ovito visualization. The boundary conditions are saved
as type.

# Arguments
- `filename`: String, the filename.
- `env`: GeneralEnv, the simulation environment.

# See also
- `save_state!`
"""
function save_state_ovito_bc!(filename, env)
    update_acc!(env)
    intact = sum(env.material_blocks[1].general.intact, dims=1)
    for mat in env.material_blocks[2:end]
        intact = hcat(intact, sum(mat.general.intact, dims=1))
    end

    damage = 1 .- reshape(intact, :) ./ env.intact0

    BC_TYPE = 99*ones(Int64, length(env.type))

    i = 1
    for bc in env.boundary_conditions
        BC_TYPE[bc.bool] .= i
        i += 1
    end

    filename = "$(dirname(filename))/BC_$(basename(filename))"

    write_data(filename; id=1:length(env.type), type=BC_TYPE, position=env.y, velocity=env.v, acceleration=env.f, mass=env.mass, volume=env.volume, damage=damage)
end


"""
    apply_bc_at0!(env, start_at)

Apply boundary conditions at t = t0.

# Arguments
- `env`: GeneralEnv, the simulation environment.
- `start_at`: Int64, the starting step.

"""
function apply_bc_at0(env, start_at)
    if start_at==0
        for bc in env.boundary_conditions
            apply_bc_at0!(env, bc)
        end
    end
end

"""
    update_acc!(env::GeneralEnv)

Updates acceleration of all the material points in a simulation environment.

̈u(xᵢ, t) = forces due to material deformation + forces due to contact

forces due to material deformation:
    > ̈u(xᵢ, t) = [∑ᵏⱼ₌₁ {T[xᵢ, t]<xⱼ-xᵢ> - T[xⱼ, t]<xᵢ-xⱼ> }*Vⱼ + b(xᵢ, t)] / ρᵢ

# Arguments
- `env`: GeneralEnv, the simulation environment.

# See also
- `force_density`
- `short_range_repulsion!`
"""
function update_acc!(env::GeneralEnv)
    # fill force vector with zeros
    fill!(env.f, 0.0)

    ###################
    # MATERIAL FORCES #
    ###################
    # calculate force density * volume i.e
    # ̈ρu(xᵢ, t) = [∑ᵏⱼ₌₁ {T[xᵢ, t]<xⱼ-xᵢ> - T[xⱼ, t]<xᵢ-xⱼ> }*Vⱼ + b(xᵢ, t)]
    # for each block (type)
    for mat in env.material_blocks
        log_detail("Force calculation for Material block $(mat.name).")
        mask = env.bid .== mat.blockid
        limits = (argmax(mask), length(mask) + 1 - argmax(reverse(mask)))
        force_density(env.f, env.y, limits, mat)
        # env.f[:, mask] .+= force_density(env.y[:, mask], mat)

        # if NaN occurs, throw error
        @check_nan env.f "env.f from material forces ($(mat.name))."
    end

    # if NaN occurs, throw error
    @check_nan env.f "env.f from material forces."

    ##################
    # CONTACT FORCES #
    ##################
    # calculate short range repulsion (should be force density * volume)
    for i in 1:size(env.short_range_repulsion, 1)
        log_detail("Force calculation for Short range repulsion $(env.short_range_repulsion[i].name).")
        short_range_repulsion!(env.y, env.f, env.type, env.bid, env.volume, env.short_range_repulsion[i])
    end

    # if NaN occurs, throw error
    @check_nan env.f "env.f from repulsive forces."

    # calculate aceeleration and momentum from force density and velocity respectively
    env.p[1, :] = 1*env.volume
    env.p[2, :] = 1*env.volume
    env.p[3, :] = 1*env.volume
    for i in 1:size(env.material_blocks,1)
        for j in env.material_blocks[i].type
            m = (env.type .== j)
            t = j - env.material_blocks[i].type[1] + 1

            density = Array(env.material_blocks[i].specific.density)[t]

            # acceleration = (force_density[force/vol/vol] * volume) / density
            env.f[:, m] ./= density

            # momentum = velocity * volume * density
            env.p[:, m] .*= env.v[:, m] .* density
        end
    end

end


"""
    function update_neighs!(envs)

Updates neighbors of each material point for a list of simulation environments.

# Arguments
- `envs`: Array{GeneralEnv}, the list of simulation environments.

# See also
- `update_repulsive_neighs!`

"""
function update_neighs!(envs; max_part=30)
    log_info("Updating neighbors for collision ...")
    for env in envs
        if env.state[1]==2
            for RM in env.short_range_repulsion
                update_repulsive_neighs!(env.y,env.type,RM; max_part=max_part)
            end
        end
    end
end

"""
    print_data_file!(envs::Array{GeneralEnv}, file_prefix::String, i::Int64)

Prints data files to disk for a list of simulation environments.

# Arguments
- `envs`: Array{GeneralEnv}, the list of simulation environments.
- `file_prefix`: String, the prefix of the output folder.
- `i`: Int64, the step number.

"""
function print_data_file!(envs::Array{GeneralEnv}, file_prefix::String, i; ext::Symbol=:jld)
    log_info("Writting data file ($i) ...")
    for id in 1:size(envs,1)
        env = envs[id]
        filename = string(file_prefix,"/env_",env.id,"_step_",i,".$(String(ext))")
        filename2 = string(file_prefix,"/env_Out.jld")
        save_state!(filename, env)
        write_global_data(filename2; Out=env.Out)
    end
end

"""
    check_boundaries!(env)

Check if any material point is outside the boundaries. If so, set the type to -1 and
velocity to zero.

# Arguments
- `env`: GeneralEnv, the simulation environment.

"""
function check_boundaries!(env)
    # _min, _max = env.boundaries
    # y = env.y

    # or(a, b) = a | b

    # mask_min = reshape(any(y .< _min, dims=1), :)
    # mask_max = reshape(any(y .> _max, dims=1), :)
    # mask = or.(mask_min, mask_max)

    # env.type[mask] .= -1
    # env.v[:, mask] .= 0.0

    return nothing
end


include("./velocity_verlet.jl")
include("./minimize.jl")
