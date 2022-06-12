export write_data

"""
    write_data(filename; kwargs...)

It writes the data file.
"""
function write_data(filename; kwargs...)
    col_names = keys(kwargs)
    col_values = values(kwargs)
    if split(filename, ".")[end] == "data"
        write_ovito(filename, col_values...; col_names=col_names)
    elseif split(filename, ".")[end] == "jld"
        items = Any[]
        for i in 1:length(col_names)
            push!(items, String(col_names[i]))
            push!(items, col_values[i]) 
        end
        save(filename, items...)
    end
end

include("./ff_ovito.jl")