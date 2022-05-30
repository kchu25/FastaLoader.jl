using FastaLoader
using Documenter

DocMeta.setdocmeta!(FastaLoader, :DocTestSetup, :(using FastaLoader); recursive=true)

makedocs(;
    modules=[FastaLoader],
    authors="Shane Kuei Hsien Chu (skchu@wustl.edu)",
    repo="https://github.com/kchu25/FastaLoader.jl/blob/{commit}{path}#{line}",
    sitename="FastaLoader.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kchu25.github.io/FastaLoader.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kchu25/FastaLoader.jl",
    devbranch="main",
)
