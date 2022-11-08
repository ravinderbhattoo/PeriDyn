export save_state!, update_acc!, update_neighs!, run!, simulate!, Solver

abstract type QuasiStaticSolver end
abstract type DynamicSolver end
Solver = Union{QuasiStaticSolver,DynamicSolver}

@def qssolver_gf begin
    max_iter::Int64
    x_tol::Float64
    f_tol::Float64
end


"""
"""
function simulate!(args...; out_dir="datafile", append_date=true, kwargs...)
    solver = last(args)
    println("Using solver: $(solver)")
    foldername = filepath_(out_dir; append_date=append_date)
    println("Output folder: $foldername")
    return run!(args...; kwargs..., out_dir=foldername)
end


function run!(envs, N::Int64, solver; filewrite_freq::Int64=10,
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
                    envs[id].Collect!(envs[id],i)
                end
            else
                println("Env $id is not active. step = $i")
            end
        end

        # calculate average property
        if i%average_prop_freq==0.0 || (i==1 && write_from==0)
            for env in envs
                println("Momentum: $i", sum(env.p, dims=2))
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

        println(round(i/N*100, digits=3),"%")
    end
end


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
"""
function save_state!(filename, env)
    update_acc!(env)
    intact = sum(env.material_blocks[1].general.intact, dims=1)
    for mat in env.material_blocks[2:end]
        intact = hcat(intact, sum(mat.general.intact, dims=1))
    end

    damage = 1.0 .- vec(intact) ./ env.intact0
    write_data(filename; id=1:length(env.type), type=env.type, position=env.y, velocity=env.v, acceleration=env.f, mass=env.mass, volume=env.volume, damage=damage)
end

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

̈u(xᵢ, t) = from material forces + from contact forces

from material force, ̈u(xᵢ, t) = [∑ᵏⱼ₌₁ {T[xᵢ, t]<xⱼ-xᵢ> - T[xⱼ, t]<xᵢ-xⱼ> }*Vⱼ + b(xᵢ, t)] / ρᵢ

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
    for i in 1:size(env.material_blocks,1)
        # mask = false
        # for j in env.material_blocks[i].type
        #     m = (env.type .== j)
        #     mask = mask .| m
        # end
        mask = env.bid .== env.material_blocks[i].blockid
        env.f[1:end, mask] .+= @timeit force_density(env.y[:, mask], env.material_blocks[i]) "material force block $i"

        # if NaN occurs, throw error
        @check_nan env.f "env.f from material forces ($(env.material_blocks[i].name))."
    end

    # if NaN occurs, throw error
    @check_nan env.f "env.f from material forces."

    ##################
    # CONTACT FORCES #
    ##################
    # calculate short range repulsion (should be force density * volume)
    for i in 1:size(env.short_range_repulsion,1)
        @timeit short_range_repulsion!(env.y, env.f, env.type, env.bid, env.volume, env.short_range_repulsion[i]) "Contact forces for RM $i"
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
            t = j - env.material_blocks[i].type.start + 1

            # acceleration = (force_density[force/vol/vol] * volume) / density
            env.f[1:end, m] ./= env.material_blocks[i].specific.density[t]

            # momentum = velocity * volume * density
            env.p[1:end, m] .*= env.v[:, m] .* env.material_blocks[i].specific.density[t]
        end
    end
end


"""
    function update_neighs!(envs)

Updates neighbors of each material point for a list of simulation environments.
"""
function update_neighs!(envs; max_part=30)
    println("\nUpdating neighbors for collision..................")
    for env in envs
        if env.state[1]==2
            for rm in 1:size(env.short_range_repulsion,1)
                update_repulsive_neighs!(env.y,env.type,env.short_range_repulsion[rm]; max_part=max_part)
            end
        end
    end
end

"""
    print_data_file!(envs::Array{GeneralEnv}, file_prefix::String, i::Int64)

Writes data file to disk.
"""
function print_data_file!(envs::Array{GeneralEnv}, file_prefix::String, i; ext::Symbol=:jld)
    println("\nWritting data file................................$i\n")
    for id in 1:size(envs,1)
        env = envs[id]
        filename = string(file_prefix,"/env_",env.id,"_step_",i,".$(String(ext))")
        filename2 = string(file_prefix,"/env_Out.jld")
        save_state!(filename, env)
        write_data(filename2; Out=env.Out)
    end
    println("Done")
end

function check_boundaries!(env)
    _min, _max = env.boundaries
    y = env.y

    or(a, b) = a | b

    mask_min = reshape(any(y .< _min, dims=1), :)
    mask_max = reshape(any(y .> _max, dims=1), :)
    mask = or.(mask_min, mask_max)

    env.type[mask] .= -1
    env.v[:, mask] .= 0.0
#     # env.volume[mask] .= 9999999

#     # env.y[1, mask] .= _min[1]
#     # env.y[2, mask] .= _min[2]
#     # env.y[3, mask] .= _min[3]

#     # env.f[1, mask] .= _min[1]
#     # env.f[2, mask] .= _min[2]
#     # env.f[3, mask] .= _min[3]
end


include("./velocity_verlet.jl")
include("./minimize.jl")
