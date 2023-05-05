export write_data

"""
    write_data(filename; kwargs...)

It writes the data file.
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

It writes the data file.
"""
function write_global_data(filename; kwargs...)
    items = Any[]
    for name in keys(kwargs)
        push!(items, String(name))
        push!(items, kwargs[name])
    end
    save(filename, items...)
end

function load_from_file!(env::GeneralEnv, filename)
    y, v, f = load("output/ttensile_sim_OSB/DSVelocityVerlet/env_1_step_100.jld", "position", "velocity", "acceleration")
    org_type = typeof(env.y)
    env.y .= org_type(y)
    env.v .= org_type(v)
    env.f .= org_type_(f)
    return nothing
end

include("./ff_ovito.jl")