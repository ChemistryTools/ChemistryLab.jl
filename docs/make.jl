using ChemistryLab
using Documenter
using PrettyTables

include("pages.jl")

DocMeta.setdocmeta!(
    ChemistryLab,
    :DocTestSetup,
    :(using ChemistryLab, DynamicQuantities, OrderedCollections, ModelingToolkit);
    recursive=true
)

makedocs(;
    modules=[ChemistryLab],
    authors="Jean-François Barthélémy and Anthony Soive",
    sitename="ChemistryLab.jl",
    format=Documenter.HTML(;
        canonical="https://jfbarthelemy.github.io/ChemistryLab.jl",
        edit_link="main",
        assets=["assets/favicon.ico"],
    ),
    pages=pages,
    warnonly=[:missing_docs, :docs_block],
)

deploydocs(; repo="github.com/jfbarthelemy/ChemistryLab.jl", devbranch="main")
