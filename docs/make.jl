using Documenter
using GasChromatographySimulator

makedocs(
    sitename = "GasChromatographySimulator",
    format = Documenter.HTML(),
    modules = [GasChromatographySimulator]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
