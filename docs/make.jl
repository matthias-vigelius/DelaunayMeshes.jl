using Documenter

push!(LOAD_PATH, "..\\src\\")

using DelaunayMeshes

makedocs(
  format = :html,
  sitename = "DelaunayMeshes",
  pages = ["User Guide" => "UserGuide.md", "API" => "index.md"]
)

deploydocs(
    repo = "github.com/matthias-vigelius/DelaunayMeshes.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
