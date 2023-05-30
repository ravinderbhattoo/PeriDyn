push!(LOAD_PATH,"../src/")

using Documenter, DocThemeIndigo
using PeriDyn, PDMaterialPoints
using Dates

indigo = DocThemeIndigo.install(PeriDyn)

makedocs(sitename="PeriDyn",
        build = "../docs",
        format = Documenter.HTML(;
            footer="Updated: $(now()). "*string(Documenter.HTML().footer),
            prettyurls = get(ENV, "CI", nothing) == "true",
            # assets = [indigo],
            assets = ["assets/font.css", "assets/purple.css"]
            ),
        # format = Documenter.LaTeX(),
        modules=[PeriDyn],
        pages = [
            "Home" => "index.md",
            # "Table of contents" => "toc.md",
            # "Material models" => "mmodels.md",
            "Contact models" => "cmodels.md",
            "Solvers" => "solvers.md",
            "Boundary Conditions" => "bc.md",
            # "Operator and Utility" => "operatorandutil.md",
            "Examples" => "examples.md",
            "Material Geometry" => "materialgeometry.md",
            # "Index" => "list.md",
            # "Autodocs" => "autodocs.md"
        ]
        )