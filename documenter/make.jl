push!(LOAD_PATH,"../src/")

using Documenter #DocThemeIndigo
using PeriDyn, PDMesh

# indigo = DocThemeIndigo.install(PeriDyn)

makedocs(sitename="PeriDyn",
        build = "../docs",
        format = Documenter.HTML(;
            prettyurls = get(ENV, "CI", nothing) == "true",
            assets=["indigo.css"]
            ),
        # format = Documenter.LaTeX(),
        modules=[PeriDyn, PDMesh],
        pages = [
            "Home" => "peridyn.md",
            "Table of contents" => "toc.md",
            "Material models" => "mmodels.md",
            "Contact models" => "cmodels.md",
            "Solvers" => "solvers.md",
            "Boundary Conditions" => "bc.md",
            "Operator and Utility" => "operatorandutil.md",
            "Examples" => "examples.md",
            "Index" => "index.md",
            "Autodocs" => "autodocs.md"
        ]
        )