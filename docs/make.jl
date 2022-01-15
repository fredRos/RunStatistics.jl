# Use
#
#     DOCUMENTER_DEBUG=true julia --color=yes make.jl local [nonstrict] [fixdoctests]
#
# for local builds.

using Documenter
using RunStatistics

# Doctest setup
DocMeta.setdocmeta!(
    RunStatistics,
    :DocTestSetup,
    :(using RunStatistics);
    recursive=true,
)

makedocs(
    sitename = "RunStatistics",
    modules = [RunStatistics],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://bat.github.io/RunStatistics.jl/stable/"
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => [
        "Introduction" => "introduction.md",
        "Guide" => "guide.md",
        ],
        "API" => "api.md",
        "LICENSE" => "LICENSE.md",
    ],
    doctest = ("fixdoctests" in ARGS) ? :fix : true,
    linkcheck = !("nonstrict" in ARGS),
    strict = !("nonstrict" in ARGS),
)

deploydocs(
    repo = "github.com/bat/RunStatistics.jl.git",
    devbranch="main",
    forcepush = true,
    push_preview = true,
)
