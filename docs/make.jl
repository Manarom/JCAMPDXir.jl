push!(LOAD_PATH,"../src/")
#include("../test/tests data/TestingData.jl")
using Documenter,JCAMPDXir
mathengine = Documenter.MathJax3()
makedocs(
        sitename = "JCAMPDXir.jl",
        repo="https://github.com/Manarom/JCAMPDXir.jl/blob/{commit}{path}#{line}",
        highlightsig = false,
        checkdocs = :none,
        format=Documenter.HTML(size_threshold = 2000 * 2^10),
        pages=[
                "Home" => "index.md",
                "JCAMPDXir" =>"JCAMPDXir.md"
                #"Examples"=>["Examples" =>"pluto_tests_git.md"
                ]
                #="Modules" => [
                    "PlanckFunctions" =>"PlanckFunctions.md"
                    "TestingData"=>"TestingData.md"
                 ] =#
)
deploydocs(;
                repo="https://github.com/Manarom/JCAMPDXir.jl/blob/{commit}{path}#{line}", 
                devbranch = "main",
                devurl="dev",
                target = "build",
                branch = "gh-pages",
                versions = ["stable" => "v^", "v#.#" ]
        )