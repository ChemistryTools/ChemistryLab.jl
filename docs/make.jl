using ChemistryLab
using Documenter
using PrettyTables

include("pages.jl")

DocMeta.setdocmeta!(
    ChemistryLab,
    :DocTestSetup,
    :(using ChemistryLab, DynamicQuantities, OrderedCollections, ModelingToolkit);
    recursive=true,
)

ENV["FORCE_COLOR"] = "true"
ENV["COLUMNS"] = "200"
ENV["LINES"] = "100"
ENV["GKSwstype"] = "100"   # headless GR backend — prevents Plots from hanging in doc builds

makedocs(;
    clean=false,
    modules=[ChemistryLab],
    authors="Jean-François Barthélémy and Anthony Soive",
    sitename="ChemistryLab.jl",
    format=Documenter.HTML(;
        canonical="https://jfbarthelemy.github.io/ChemistryLab.jl",
        edit_link="main",
        assets=["assets/favicon.ico", "assets/custom.css"],
        prettyurls=(get(ENV, "CI", nothing) == "true"),
        collapselevel=1,
    ),
    pages=pages,
    warnonly=[:missing_docs, :docs_block],
    draft=false,
)

deploydocs(; repo="github.com/jfbarthelemy/ChemistryLab.jl", devbranch="main")
