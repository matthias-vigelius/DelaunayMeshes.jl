if (!(pwd() in LOAD_PATH))
  push!(LOAD_PATH, pwd())
end

include("QuadEdgeTests.jl")
include("TriangulationTests.jl")
include("DelaunayMeshTests.jl")
