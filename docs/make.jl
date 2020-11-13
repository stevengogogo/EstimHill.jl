using EstimHill
using Documenter

makedocs(;
    modules=[EstimHill],
    authors="Steven Chiu <stevengogogo4321@gmail.com> and contributors",
    repo="https://github.com/stevengogogo/EstimHill.jl/blob/{commit}{path}#L{line}",
    sitename="EstimHill.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://stevengogogo.github.io/EstimHill.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/stevengogogo/EstimHill.jl",
)
