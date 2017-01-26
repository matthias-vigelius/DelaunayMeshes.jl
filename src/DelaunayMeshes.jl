"""
This module provides high-level functionality to generate, manage and refine
meshes based on constrained Delaunay triangulations.
"""

module DelaunayMeshes

import DiffEqPDEBase
import GeometricalPredicates

import Base.convert
import Base.push!


include("QuadEdge.jl")
include("Triangulation.jl")
include("DelaunayMesh.jl")

export getvoronoivertices, setboundingbox
export scalepoints, unscalepoints, addconstraint!
export getdelaunaycoordinates
end # module
