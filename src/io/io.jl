export write_data

"""
    write_data(filename; kwargs...)

Writes the data file.

Arguments:
- `filename`: String, the name of the file to be written.
- `kwargs`: Keyword arguments, additional options for writing the file.

Note: This function supports writing files in the `.data` and `.jld` formats.
"""
function write_data(filename; kwargs...)
    if split(filename, ".")[end] == "data"
        write_ovito(filename; kwargs...)
    elseif split(filename, ".")[end] == "jld"
        items = Any[]
        for name in keys(kwargs)
            push!(items, String(name))
            push!(items, Array(kwargs[name]))
        end
        save(filename, items...)
    end
end


"""
    write_global_data(filename; kwargs...)

Writes the global data file.

Arguments:
- `filename`: String, the name of the file to be written.
- `kwargs`: Keyword arguments, additional options for writing the file.

Note: This function writes a data file in the `.jld` format.
"""
function write_global_data(filename; kwargs...)
    items = Any[]
    for name in keys(kwargs)
        push!(items, String(name))
        push!(items, kwargs[name])
    end
    save(filename, items...)
end

"""
    load_from_file!(env::GeneralEnv, filename)

Loads data from a file into the specified environment.

Arguments:
- `env`: GeneralEnv, the environment to load the data into.
- `filename`: String, the name of the file to load data from.

Returns:
- Nothing

Note: This function assumes the presence of specific data fields in the file, such as "position", "velocity", and "acceleration".
"""
function load_from_file!(env::GeneralEnv, filename)
    y, v, f = load("output/ttensile_sim_OSB/DSVelocityVerlet/env_1_step_100.jld", "position", "velocity", "acceleration")
    org_type = typeof(env.y)
    env.y .= org_type(y)
    env.v .= org_type(v)
    env.f .= org_type_(f)
    return nothing
end

include("./ff_ovito.jl")