export write_ovito, write_ovito_cell_ids, jld2ovito

"""
    jld2ovito(file, N; start=0, step=100)

Converts JLD files to Ovito data files.

Arguments:
- `file`: String, the base name of the JLD files.
- `N`: Int, the number of files to convert.
- `start`: Int, the index of the first file to convert. Default is 0.
- `step`: Int, the step size between files to convert. Default is 100.

Note: This function iterates over a range of files and converts each JLD file to an Ovito data file using the `write_ovito` function.
"""
function jld2ovito(file, N; start=0, step=100)
    for i in start:step:N
        x = replace(file, "*"=>"$i")
        log_info("$(x)")
        out = jldread(x)
        write_ovito(x*".data"; out...)
    end
end

"""
    jld2array(file, N; start=0, step=100)

Loads data from JLD files into an array.

Arguments:
- `file`: String, the base name of the JLD files.
- `N`: Int, the number of files to load.
- `start`: Int, the index of the first file to load. Default is 0.
- `step`: Int, the step size between files to load. Default is 100.

Returns:
- Array, an array containing the data loaded from the JLD files.

Note: This function iterates over a range of files and loads data from each JLD file using the `jldread` function.
"""
function jld2array(file, N; start=0, step=100)
    out = []
    for i in start:step:N
        x = replace(file, "*"=>"$i")
        log_info("$(x)")
        push!(out, jldread(x))
    end
    return out
end


"""
    jldread(filename::String)

Reads data from a JLD file.

Arguments:
- `filename`: String, the name of the JLD file to read.

Returns:
- Dict, a dictionary containing the data read from the JLD file.

Note: This function uses the `load` function from the JLD package to read the data from the file, and converts it to a dictionary format.
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

Writes data to an Ovito data file.

Arguments:
- `filename`: String, the name of the file to write.
- `kwargs`: Keyword arguments, the data to write to the file.

Note: This function writes data in Ovito data file format, with each column specified by a keyword argument.
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
    write_ovito_cell_ids(filename::String, y::Matrix, horizon)

Writes cell IDs to an Ovito data file.

Arguments:
- `filename`: String, the name of the file to write.
- `y`: Matrix, the coordinates of the particles.
- `horizon`: The horizon value.

Note: This function writes the cell IDs to an Ovito data file, based on the particle coordinates and horizon value.
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

