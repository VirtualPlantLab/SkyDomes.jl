using Sky
using Documenter

makedocs(;
    doctest = false,
    modules = [Sky],
    authors = "Alejandro Morales Sierra <alejandro.moralessierra@wur.nl> and contributors",
    repo = "https://github.com/VirtualPlantLab/Sky.jl/blob/{commit}{path}#{line}",
    sitename = "Sky.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        edit_link = "master",
        assets = String[]),
    pages = [
        "API" => "index.md",
    ])

deploydocs(;
    repo = "github.com/VirtualPlantLab/Sky.jl.git",
    devbranch = "master")
