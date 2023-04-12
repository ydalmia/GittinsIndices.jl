using GittinsIndices
using Documenter

DocMeta.setdocmeta!(GittinsIndices, :DocTestSetup, :(using GittinsIndices); recursive=true)

makedocs(;
    modules=[GittinsIndices],
    authors="Yash Dalmia, Santino Ramos, Amelie Buc",
    repo="https://github.com/ydalmia/GittinsIndices.jl/blob/{commit}{path}#{line}",
    sitename="GittinsIndices.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ydalmia.github.io/GittinsIndices.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ydalmia/GittinsIndices.jl",
    devbranch="main",
)
