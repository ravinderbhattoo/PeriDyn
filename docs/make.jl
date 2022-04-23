using Documenter, DocThemeIndigo
using PeriDyn, PDMesh

indigo = DocThemeIndigo.install(PeriDyn)

makedocs(sitename="PeriDyn", 
        modules=[PeriDyn, PDMesh],
        format = Documenter.HTML(;
            prettyurls = false,
            assets=String[indigo]),
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