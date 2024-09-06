using Pkg
Pkg.activate(".")
Pkg.instantiate()

function load_packages()
    Pkg.add("Documenter")
    Pkg.add("DocThemeIndigo")
    Pkg.add("Revise")
    Pkg.add("Dates")

    Pkg.develop(path="../../PDMaterialPoints/")
    Pkg.develop(path="../")
    Pkg.instantiate()
end

try
    load_packages()
catch
    @warn "Could not load packages. Re-doing..."
finally
    rm("Project.toml")
    rm("Manifest.toml")
    load_packages()
end

using Revise
using Documenter, DocThemeIndigo
using PeriDyn, PDMaterialPoints
using Dates

indigo = DocThemeIndigo.install(PeriDyn)

format = Documenter.HTML(;
    footer="Updated: $(mapreduce(x-> " "*x*" ", *, split(string(now()), "T")))",
    #. "*string(Documenter.HTML().footer),
    prettyurls = get(ENV, "CI", nothing) == "true",
    # assets = [indigo],
    assets = ["assets/font.css", "assets/color.css"],
    size_threshold = 1024 * 1024 * 1 # MB

)

# format = Documenter.LaTeX(;platform="none")

function doit()
    makedocs(sitename="PeriDyn",
    build = "../docs",
    format = format,
    # remotes = nothing,
    modules=[PeriDyn],
    pages = [
        "Home" => "index.md",
        "Quick start" => "examples.md",
        "Table of contents" => "toc.md",
        "Material geometry" => "materialgeometry.md",
        "Material models" => "mmodels.md",
        "Contact models" => "cmodels.md",
        "Boundary conditions" => "bc.md",
        "Solvers" => "solvers.md",
        # "Operator and Utility" => "operatorandutil.md",
        "I/O" => "io.md",
        "Index" => "list.md",
        "Autodocs" => "autodocs.md"
    ]
    )
end

doit()

