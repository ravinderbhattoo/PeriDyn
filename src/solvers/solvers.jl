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

function Base.show(io::IO, solver::T) where T <: Solver
    txt = ""
    for k in fieldnames(T)
        txt = txt * " $k=$(getfield(solver, k)),"
    end
    txt = @green("$T") * ": [" * txt[2:end-1] * "]"
    print(io, txt)
end

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

# See also
[`run!`](@ref)
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
- `envs`: Array{Environment}, the list of environments.
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
[`simulate!`](@ref)
"""
function run!(envs, N::Int64, solver::Solver; filewrite_freq::Int64=10,
                neigh_update_freq::Int64=1, average_prop_freq::Int64=1,
                out_dir::String="datafile",
                start_at::Int64=0, write_from::Int=0, ext::Symbol=:jld,
                max_part=30)

    if isa(envs, GeneralEnvironment)
        envs = [envs]
    end

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

    # run simulation
    log_impinfo("\n \
                    Simulation started at step $(start_at).\n \
                    Total number of steps are $N.\n \
                    System sizes are $(mapreduce(x->string(size(x.y))*" ", *, envs)).\n \
                    Output directory: $out_dir\n ")

    filename = joinpath("$out_dir", "log.txt")

    # print average property header if GOBBLE_KEYS[] is true
    PeriDyn.LOGFILE_HANDLE = open(filename, "w")
    write_headers(envs)
    close(PeriDyn.LOGFILE_HANDLE)

    for i in (1+start_at):N

        # apply integrator step for all envs
        for id in 1:size(envs,1)
            if envs[id].state==2
                @timeit timings "applying solver env $id" apply_solver!(envs[id], solver)
                if envs[id].Collect! !== nothing
                    @timeit timings "Collect! env %id" envs[id].Collect!(envs[id], i)
                end
            else
                log_impinfo("Env $id is not active. step = $i")
            end
        end

        # print average property
        if i%average_prop_freq==0.0 || (i==1 && write_from==0)
            PeriDyn.LOGFILE_HANDLE = open(filename, "a")
            log_data(Step=i)
            for env in envs
                printdata(env)
            end
            println()
            write(PeriDyn.LOGFILE_HANDLE, "\n")
            close(PeriDyn.LOGFILE_HANDLE)
        end

        # write file to disk
        if i%filewrite_freq==0.0 || (i==1 && write_from==0)
            @timeit timings "File write" _print_data_file!(i+write_from)

        end

        # update neighbors
        if i%neigh_update_freq==0.0
            @timeit timings "Neigh update" update_neighs!(envs; max_part=max_part)
        end

        log_detail("Solver iteration : $i")

    end

    _print_data_file!("end_with")

    log_impinfo("Runs finished.")
end


"""
    apply_bc_at0!(env, start_at)

Apply boundary conditions at t = t0.

# Arguments
- `env`: the simulation environment.
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
    update_mat_acc!(env)

Update the forces (acceleration) due to material deformation.

```math
\\rho \\ddot{u}(x_i, t) = \\left[\\sum_{j=1}^{N_j} \\left\\{T[x_i, t]\\langle x_j-x_i \\rangle - T[x_j, t]\\langle x_i-x_j \\rangle \\right\\}V_j + b(x_i, t)\\right]
```
where ``T`` is the force density, ``V_j`` is the volume of the ``j``th particle,
``\\rho`` is the density, ``u`` is the displacement, and ``b`` is the body force.
The summation is over all the particles in the neighborhood given by ``N_j``.

# Arguments
- `env`: the simulation environment.

# See also
[`force_density!`](@ref)

"""
function update_mat_acc!(env)
    ###################
    # MATERIAL FORCES #
    ###################
    # calculate force density * volume given by formula mentioned above
    # for each block (type)
    for mat in env.material_blocks
        log_detail("Force calculation for Material block $(mat.name).")
        mask = env.bid .== mat.blockid
        limits = (argmax(mask), length(mask) + 1 - argmax(reverse(mask)))
        @timeit timings mat.name force_density!(env.f, env.y, limits, mat) # (force_density[force/vol/vol] * volume)

        # if NaN occurs, throw error. can take ms on large array so do once in an while
        # @check_nan env.f "env.f from material forces ($(mat.name))."
    end
end


"""
    update_contact_acc!(env)

Update the forces (acceleration) due to contact.

# Arguments
- `env`: the simulation environment.

# See also
[`short_range_contact!`](@ref)
"""
function update_contact_acc!(env)
    ##################
    # CONTACT FORCES #
    ##################
    # calculate short range contact (should be force density * volume)
    for contact_model in env.short_range_contact
        log_detail("Force calculation for Short range contact $(contact_model.name).")
        @timeit timings contact_model.name short_range_contact!(env.y, env.f, env.type, env.bid, env.volume, env.density, contact_model)
    end

    # if NaN occurs, throw error
    # @check_nan env.f "env.f from contact forces."
end

"""
    update_misc!(env)

Update the misc items such as momentum etc.

# Arguments
- `env`: the simulation environment.

"""
function update_misc!(env)
    # calculate aceeleration and momentum from force density and velocity respectively

    Threads.@threads for i in 1:PeriDyn.SPATIAL_DIMENSIONS_REF[]
        env.p[i, :] .= env.mass .* @view env.v[i, :]
    end

    # for i in 1:size(env.material_blocks,1)
        # for j in env.material_blocks[i].type
            # m = (env.type .== j)
            # f = @view env.f[:, m]
            # p = @view env.p[:, m]
            # v = @view env.v[:, m]

            # t = j - env.material_blocks[i].type[1] + 1
            # density = env.material_blocks[i].specific.density[t]

            # acceleration = (force_density[force/vol/vol] * volume) / density
            # f .*= (1/density)

            # momentum = velocity * volume * density
            # p .*= v .* density
    #     end
    # end
end


"""
    update_acc!(env)

Updates acceleration of all the material points in a simulation environment.

# Arguments
- `env`: the simulation environment.

# See also
[`update_mat_acc!`](@ref)
[`update_contact_acc!`](@ref)
[`update_misc!`](@ref)
"""
function update_acc!(env)
    # fill force vector with zeros
    fill!(env.f, eltype(env.f)(0.0))

    ###################
    # MATERIAL FORCES #
    ###################
    @timeit timings "Material force cal" update_mat_acc!(env)


    ##################
    # CONTACT FORCES #
    ##################
    @timeit timings "Contact force cal" update_contact_acc!(env)

    @timeit timings "Misc cal" update_misc!(env)
end


"""
    function update_neighs!(envs)

Updates neighbors of each material point for a list of simulation environments.

# Arguments
- `envs`: Array{Environment}, the list of simulation environments.

# See also
- `update_contact_neighs!`

"""
function update_neighs!(envs; max_part=30)
    log_info("Updating neighbors for collision ...")
    for env in envs
        if env.state[1]==2
            for RM in env.short_range_contact
                log_detail("Updating neighbors for $(RM.name) ...")
                @timeit timings "NU $(RM.name)" update_contact_neighs!(env.y,env.type,RM; max_part=max_part)
            end
        end
    end
end

"""
    check_boundaries!(env)

Check if any material point is outside the defined boundaries. Do nothing.

# Arguments
- `env`: the simulation environment.

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
