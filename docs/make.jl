using Documenter, VPL, DocumenterMarkdown

makedocs(;
    modules = [Sky],
    format = Markdown(),#Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages=[
        "API" => "Sky.md",
    ],
    repo="https://github.com/AleMorales/Sky.jl/blob/{commit}{path}#L{line}",
    sitename="Sky.jl",
    authors="Alejandro Morales, Wageningen University & Research"
)

# deploydocs(;
#     repo="github.com/AleMorales/VPL.jl",
# )
