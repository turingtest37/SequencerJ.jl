using Documenter
using SequencerJ

makedocs(
    sitename = "SequencerJ",
    format = Documenter.HTML(),
    modules = [SequencerJ],
    pages = [
        "ABCs of SequencerJ" => "index.md",
        "Metrics" => "algorithms.md",
        "Examples" => "examples.md",
        "API" => "api.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
