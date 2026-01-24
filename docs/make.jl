using Documenter
using GasChromatographySimulator
using ForwardDiff, Measurements

makedocs(
            sitename = "GasChromatographySimulator.jl",
            format = Documenter.HTML(
                size_threshold_warn = 500_000,
                size_threshold = 2_000_000,
                example_size_threshold = 20_000,
            ),
            pages = Any[
                "Home" => "index.md",
                "Installation" => "installation.md",
                "Usage" => "usage.md",
                "ODE Solvers" => "ode_solvers.md",
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