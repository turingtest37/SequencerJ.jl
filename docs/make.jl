using Documenter
using SequencerJulia

makedocs(
    sitename = "SequencerJulia",
    format = Documenter.HTML(),
    modules = [SequencerJulia],
    pages = [
        "Quick Start" => "index.md",
        "Metrics" => "algorithms.md",
        "API" => "api.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
