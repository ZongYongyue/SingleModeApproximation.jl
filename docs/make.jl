using Documenter
using MeanFieldTheories

makedocs(;
    modules=[MeanFieldTheories],
    authors="zongyy",
    repo="https://github.com/zongyy/MeanFieldTheories.jl/blob/{commit}{path}#{line}",
    sitename="MeanFieldTheories.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://zongyy.github.io/MeanFieldTheories.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Quantum System" => "quantumsystem.md",
        "Hartree-Fock Approximation" => [
            "Real Space"      => "hartreefock_real.md",
            "Momentum Space"  => "hartreefock_momentum.md",
        ],
        "Single-Mode Approximation" => "singlemode.md",
    ],
    checkdocs=:none,
    warnonly=true,
)

deploydocs(;
    repo="github.com/ZongYongyue/MeanFieldTheories.jl",
    devbranch="main",
    push_preview=true,
)
