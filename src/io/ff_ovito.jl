export write_ovito, write_ovito_cell_ids, jld2ovito

"""
    jld2ovito(file, N; step=100)

It writes the data file.
"""
function jld2ovito(file, N; start=0, step=100)
    for i in start:step:N
        x = replace(file, "*"=>"$i")
        println(x)
        out = jldread(x)
        write_ovito(x*".data"; out...)
    end
end

"""
    jld2array(file, N; step=100)

It loads the data file.
"""
function jld2array(file, N; start=0, step=100)
    out = []
    for i in start:step:N
        x = replace(file, "*"=>"$i")
        println(x)
        push!(out, jldread(x))
    end
    return out
end


"""
    jldread(filename::String, args...)

It reads the jld data file.
"""
function jldread(filename::String)
    inp = Dict()
    out = load(filename)
    for pair in out
        key, value = pair
        merge!(inp, Dict(Symbol(key) => value))
    end
    return inp
end

"""
    write_ovito(filename::String; kwargs...)

It writes the data file.
"""
function write_ovito(filename::String; kwargs...)
    col_names = sort(collect(keys(kwargs)))
    args = [Array(kwargs[k]) for k in col_names]
    item = kwargs[:id]

    if length(size(item)) > 1
        N = size(item, 2)
    else
        N = length(item)
    end

    if ~isa(col_names, Nothing)
        comment="# " *mapfoldl(String, (a, b)->"$a | $b", col_names)
    else
        comment=""
    end
    file = open(filename, "w+")
    function getitems(i, arg)
        if length(size(arg)) == 1
            return arg[i]
        else
            return arg[:, i]
        end
    end

    write(file, "$N \n$(comment)\n")

    CUDA.allowscalar() do
        for j in 1:N
            row = foldr((a, b) -> "$a $b", foldr((a, b) -> (a..., b...), [getitems(j, arg) for arg in args]))
            write(file, "$row \n")
        end
    end
    close(file)
end


"""
    write_ovito_cell_ids(filename::String, y::Array{Float64,2}, cells::Any)

It writes the data file.
"""
function write_ovito_cell_ids(filename::String, y::Matrix, horizon)
    cells, _ = get_cells(y, horizon)
    file = open(filename, "w+")
    N = size(y, 2)
    write(file, "$N \n\n")
    ind = 1
    for i1 in 1:length(cells)
        for i2 in 1:length(cells[i1])
            for i3 in 1:length(cells[i1][i2])
                ids = cells[i1][i2][i3]
                for j in 1:length(ids)
                    a, b, c = y[1, ids[j]], y[2, ids[j]], y[3, ids[j]]
                    write(file, "$ind, $a, $b, $c \n")
                end
                ind += 1
            end
        end
    end
    close(file)
end

