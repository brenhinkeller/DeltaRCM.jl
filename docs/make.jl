using DeltaRCM
using Documenter

DocMeta.setdocmeta!(DeltaRCM, :DocTestSetup, :(using DeltaRCM); recursive=true)

makedocs(;
    modules=[DeltaRCM],
    authors="C. Brenhin Keller <cbkeller@dartmouth.edu> and contributors",
    sitename="DeltaRCM.jl",
    checkdocs=:exports,
    format=Documenter.HTML(;
        canonical="https://brenhinkeller.github.io/DeltaRCM.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Guide" => "guide.md",
        "Output and checkpoints" => "output.md",
        "Customization" => "customization.md",
        "Differences from pyDeltaRCM" => "differences.md",
        "API" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/brenhinkeller/DeltaRCM.jl",
    devbranch="main",
)
