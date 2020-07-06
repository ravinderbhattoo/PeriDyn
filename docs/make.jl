using Documenter, PeriDyn

makedocs(sitename = "PeriDyn",
format = Documenter.HTML(prettyurls = true),
pages = [
    "Table of Contents" => "toc.md",
    "PeriDyn" => "index.md",
    "Examples" => "zexamples.md",
    "Index" => "list.md",
    "Autodocs" => "zzfull.md",
]
)
