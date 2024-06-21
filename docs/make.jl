using Documenter
using GasChromatographySimulator
using ForwardDiff, Measurements

makedocs(
            sitename = "GasChromatographySimulator.jl",
            pages = Any[
                "Home" => "index.md",
                "Installation" => "installation.md",
                "Usage" => "usage.md",
                "Examples" => "examples.md",
                "Additional features" => "additional.md",
                "Functions" => "functions.md",
                "References" => "references.md"
            ]
        )

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/GasChromatographyToolbox/GasChromatographySimulator.jl",
    devbranch = "main"
)