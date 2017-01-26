using Documenter

push!(LOAD_PATH, "..\\src\\")

using QuadEdge
using Triangulation


makedocs(
  format = :html,
  sitename = "DelaunayTriangulation",
  pages = ["index.md"]
)
