"""
    Functions for printing and displaying variables
"""


"""
    array_repr(item)

Represent an array in nice format.

# Arguments
- `item`: Any, the array.

"""
function array_repr(item)
    txt = ""
    et = eltype(item)
    s = size(item)
    if eltype(item) <: Quantity
        txt = txt * "$(unit(et)) ($(dimension(et))) "
    end
    if prod(size(item)) > 16
        txt = txt * "size$s\n"
    elseif eltype(item) <: Quantity
        str_ = split(repr("text/plain", item), ":")[2]
        txt = txt * "$(str_) \n"
    else
        try
            str_ = split(repr("text/plain", item), ":")[2]
            txt = txt * "size$s $(str_)\n"
        catch e
            txt = txt * "\n $(repr("text/plain", item))\n"
        end
    end
    return txt
end


"""
    variable_txt(item)

Represent a variable in nice format.

# Arguments
- `item`: Any, the variable.

"""
function variable_txt(item)
    if isa(item, AbstractArray)
        txt = array_repr(item)
    else
        txt = "$(item)\n"
    end
    return txt
end


"""
    variable_color(x; kalar="#773399")

Color a variable.

# Arguments
- `x`: Any, the variable to be colored.

# Keyword Arguments
- `kalar`: String, the color to use.

"""
function variable_color(x; kalar="#773399") #883355
    name = @bold "{$kalar}$(string(x)){/$kalar}"
    return name
end

"""
    DPanel(x; title_style="bold blue", title_justify=:center, kwargs...)

Display a variable in a panel.

# Arguments
- `x`: Any, the text to be displayed.

# Keyword Arguments
- `title_style`: String, the style of the title.
- `title_justify`: Symbol, the justification of the title.
- `kwargs...`: Any, keyword arguments to be passed to `Panel`.

"""
function DPanel(x; title_style="bold blue",
    title_justify=:center,
    kwargs...)
Panel(x; title_style=title_style,
title_justify=title_justify,
kwargs...)
end

"""
    SPanel(x; title_style="bold green", title_justify=:center, kwargs...)

Display a variable in a panel.

# Arguments
- `x`: Any, the text to be displayed.

# Keyword Arguments
- `title_style`: String, the style of the title.
- `title_justify`: Symbol, the justification of the title.
- `kwargs...`: Any, keyword arguments to be passed to `Panel`.

"""
function SPanel(x; title_style="bold green",
    title_justify=:center,
    style="blue bold",
    justify=:center,
    kwargs...)
Panel(x; title_style=title_style,
title_justify=title_justify,
style=style,
justify=justify,
kwargs...)
end