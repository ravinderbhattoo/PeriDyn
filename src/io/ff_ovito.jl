export write_ovito, write_ovito_cell_ids

"""
    write_ovito(filename::String, args...)

It writes the data file.
"""
function write_ovito(filename::String)
    inp = Dict()
    out = Dict(load(filename))
    for item in out
        merge!(inp, Dict(Symbol(item[1]) => item[2]))
    end    
    write_data(filename*".data"; inp...)
end

"""
    write_ovito(filename::String, args...)

It writes the data file.
"""
function write_ovito(filename::String, args...; col_names=nothing)
    if ~isa(col_names, Nothing)
        comment="# " *mapfoldl(String, (a, b)->"$a | $b", col_names)
    else
        comment=""
    end
    file = open(filename, "w+")
    if length(size(args[1])) == 1
        N = length(args[1])
    else
        N = size(args[1], 2)
    end
    function getitems(i, arg)
        if length(size(arg)) == 1
            return arg[i]
        else
            return arg[:, i]
        end
    end

    write(file, "$N \n$(comment)\n")
    for j in 1:N
        row = foldr((a, b) -> "$a $b", foldr((a, b) -> (a..., b...), [getitems(j, arg) for arg in args]))
        write(file, "$row \n")
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

