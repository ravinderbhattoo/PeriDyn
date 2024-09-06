export filepath_, save_state!, save_state_ovito_bc!, print_data_file!
export write_data, write_global_data, load_from_file!, save_to_file



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
        return mkpath(file_prefix*"_"*dtime)*"/"
    else
        return mkpath(file_prefix)*"/"
    end
end


"""
    save_state!(filename, env)

Save `env` to disk.

# Arguments
- `filename`: String, the filename.
- `env`: GeneralEnvironment, the simulation environment.

# Keyword Arguments
- `force`: Bool, the boolean value to force saving to data file format. Default is `false`.

# See also
- `save_state_ovito_bc!`
"""
function save_state!(filename, env; force=false)
    update_acc!(env)
    damage = cal_damage(env)
    if filename[end-4:end] == ".data"
        if (size(env.y, 2) > 500000) && !force
            filename = filename[1:end-5] * ".jld2"
            log_warning("Saving to JLD2 file instead of data file. Size of the data is too large.")
        end
    end
    write_data(filename; id=1:length(env.type), type=env.type, position=env.y, velocity=env.v, acceleration=env.f, mass=env.mass, volume=env.volume, damage=damage)
end

"""
    save_state_ovito_bc!(filename, env)

Save `env` to disk for ovito visualization. The boundary conditions are saved
as type.

# Arguments
- `filename`: String, the filename.
- `env`: GeneralEnvironment, the simulation environment.

# Keyword Arguments
- `force`: Bool, the boolean value to force saving to data file format. Default is `false`.

# See also
- `save_state!`
"""
function save_state_ovito_bc!(filename, env; force=false)
    update_acc!(env)
    damage = cal_damage(env)

    BC_TYPE = 99*ones(Int64, length(env.type))

    i = 1
    for bc in env.boundary_conditions
        BC_TYPE[bc.bool] .= i
        i += 1
    end

    filename = "$(dirname(filename))/BC_$(basename(filename))"

    if filename[end-4:end] == ".data"
        if (size(env.y, 2) > 500000) && !force
            filename = filename[1:end-5] * ".jld2"
            log_warning("Saving to JLD2 file instead of data file. Size of the data is too large.")
        end
    end
    write_data(filename; id=1:length(env.type), type=BC_TYPE, position=env.y, velocity=env.v, acceleration=env.f, mass=env.mass, volume=env.volume, damage=damage)
end


"""
    print_data_file!(envs::Array{GeneralEnvironment}, file_prefix::String, i::Int64)

Prints data files to disk for a list of simulation environments.

# Arguments
- `envs`: Array{GeneralEnvironment}, the list of simulation environments.
- `file_prefix`: String, the prefix of the output folder.
- `i`: Int64, the step number.

"""
function print_data_file!(envs::Array{GeneralEnvironment}, file_prefix::String, i; ext::Symbol=:jld)
    log_info("Writting data file ($i) ...")
    for id in 1:size(envs,1)
        env = envs[id]
        filename = joinpath(file_prefix, "env_$(env.id)_step_$i.$(String(ext))")
        filename2 = joinpath(file_prefix,"env_Out.jld2")
        save_state!(filename, env)
        @assert isa(env.Out, Dict) @red @bold "Cannot save global data. \
                env.Out is not a Dict. Please collect data in env.Out as Dict."
        write_global_data(filename2; env.Out...)
    end
end



"""
    write_data(filename; kwargs...)

Writes the data file.

# Arguments
- `filename`: String, the name of the file to be written.
- `kwargs`: Keyword arguments, additional options for writing the file.

Note: This function supports writing files in the `.data` and `.jld2` formats.
"""
function write_data(filename; kwargs...)
    if split(filename, ".")[end] == "data"
        write_ovito(filename; kwargs...)
    elseif split(filename, ".")[end] == "jld2"
        jldsave(filename; kwargs...)
    elseif split(filename, ".")[end] == "jld"
        jldsave(filename; kwargs...)
    end
end

"""
    write_global_data(filename; kwargs...)

Writes the global data file.

# Arguments
- `filename`: String, the name of the file to be written.
- `kwargs`: Keyword arguments, additional options for writing the file.

Note: This function writes a data file in the `.jld2` format.
"""
function write_global_data(filename; kwargs...)
    jldsave(filename; kwargs...)
end


"""
    load_from_file!(filename)

Loads data from a file into the specified environment.

# Arguments
- `filename`: String, the name of the file to load data from.

# Returns
- `env`: the simulation environment.
"""
function load_from_file(filename)
    return load(filename)["env"]
end


"""
    save_to_file(env, filename)

Saves data from the specified environment to a file.

# Arguments
- `env`: GeneralEnvironment, the environment to save.
- `filename`: String, the name of the file to save data to.

"""
function save_to_file(env, filename)
    jldsave(filename; env=env)
end


include("./ff_ovito.jl")