using Documenter
using SingleModeApproximation

makedocs(;
    modules=[SingleModeApproximation],
    authors="zongyy",
    repo="https://github.com/zongyy/SingleModeApproximation.jl/blob/{commit}{path}#{line}",
    sitename="SingleModeApproximation.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://zongyy.github.io/SingleModeApproximation.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Quantum System" => "quantumsystem.md",
        "Hartree-Fock Approximation" => "hartreefock.md",
        "Single-Mode Approximation" => "singlemode.md",
    ],
)

deploydocs(;
    repo="github.com/ZongYongyue/SingleModeApproximation.jl",
    devbranch="main",
    push_preview=true,
)
