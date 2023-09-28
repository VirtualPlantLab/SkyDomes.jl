using SkyDomes
using Documenter

makedocs(;
    doctest = false,
    modules = [SkyDomes],
    authors = "Alejandro Morales Sierra <alejandro.moralessierra@wur.nl> and contributors",
    repo = "https://github.com/VirtualPlantLab/SkyDomes.jl/blob/{commit}{path}#{line}",
    sitename = "SkyDomes.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        edit_link = "master",
        assets = String[]),
    pages = [
        "API" => "index.md",
    ])

deploydocs(;
    repo = "github.com/VirtualPlantLab/SkyDomes.jl.git",
    devbranch = "master")
