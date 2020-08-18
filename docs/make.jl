using Solar
using Documenter

makedocs(;
    modules=[Solar],
    authors="Alejandro Morales Sierra <morales.s.alejandro@gmail.com> and contributors",
    repo="https://github.com/AleMorales/Solar.jl/blob/{commit}{path}#L{line}",
    sitename="Solar.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://AleMorales.github.io/Solar.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/AleMorales/Solar.jl",
)
