typealias Constraint Vector{Int}

"""
    computescaleandshiftfromboundingbox(boundingBox::Vector{Float64})

Computes the scale and horizontal/vertical shift from a given bounding box.
This establishes an aspect-ratio-conserving transformation from the given
bounding box into a square that fits comfortably into the outer bounding triangle
of the DelaunayMeshes. The square (blue) fills ten per cent of the maximally
inscribed square (dashed).

![inscribed](inscribedsquare.svg)

Specifically, if ``(c_x, c_y)`` denote the center of the bounding box and
``\Delta = \max((x_\mathrm{max} - x_\mathrm{min}), (x_\mathrm{max} - x_\mathrm{min}))``
is the maximum extent of the bounding box, then the coordinate transformation to
the new primed coordinates is given by
``(x', y') = 0.1 (x - c_x, y - c_y)/\Delta + (c_x', c_y')``,
where ``(c_x', c_y') = (1.5, 1.25)`` denotes the center of the inscribed square.
"""
function computescaleandshiftfromboundingbox(boundingBox::Vector{Float64})
  local scale = max(boundingBox[2]-boundingBox[1], boundingBox[4]-boundingBox[3])
  local cx = 0.5*(boundingBox[2]+boundingBox[1])
  local cy = 0.5*(boundingBox[4]+boundingBox[3])
  return scale, cx, cy
end

"""
Defines a mesh. A mesh consists of a bounding box providing bounds to the
vertex coordinates, constraints consisting of closed polygons and the
corresponding Delaunay tesselation.
"""
type Mesh
   boundingBox::Nullable{Vector{Float64}}
   scale::Nullable{Float64}
   shiftX::Nullable{Float64}
   shiftY::Nullable{Float64}

   tesselation::DelaunayMeshes.DelaunayTesselation
   constraints::Vector{Constraint}
   faceLocation::Vector{Bool}
   voronoiVertices::Nullable{Array{Float64, 2}}

   Mesh() = new(
     Nullable{Vector{Float64}}(),
     Nullable{Vector{Float64}}(),
     Nullable{Vector{Float64}}(),
     Nullable{Vector{Float64}}(),
     DelaunayMeshes.DelaunayTesselation(),
     Vector{Constraint}(),
     [true, false],
     Nullable{Array{Float64, 2}}()
   )
end

"""
    getvoronoivertices(mesh::Mesh)

Computes the Voronoi vertices of the given DelaunayMeshes.

# Remarks
* The results are cached but need to be re-computed whenever a new vertex or
  a new constraint is pushed into the DelaunayMeshes.
"""
function getvoronoivertices(mesh::Mesh)
  if (!isnull(mesh.voronoiVertices))
    return mesh.voronoiVertices
  end

  mesh.voronoiVertices = fill(Inf, size(mesh.tesselation.faces, 1), 2)

  # go through all dual edges and compute barycenter of corresponding face
  for ei in 2:2:length(mesh.tesselation.edges)
    local curFace = mesh.tesselation.edges[ei].origin
    if get(mesh.voronoiVertices)[curFace, :] == [Inf, Inf]
      local curFaceEdge = rot(ei)
      local curFaceVertexIndices =
       lface_vertices(mesh.tesselation, curFaceEdge)
      local curFaceVertexCoordinates =
       map((i::Int) -> mesh.tesselation.vertices[i], curFaceVertexIndices)
      local baryCenter = sum(curFaceVertexCoordinates)/3.
      get(mesh.voronoiVertices)[curFace, :] = baryCenter
    end
  end

  return get(mesh.voronoiVertices)
end

function setboundingbox(mesh::Mesh, boundingBox::Vector{Float64})
  @assert isnull(mesh.boundingBox)
  mesh.boundingBox = boundingBox
  mesh.scale, mesh.shiftX, mesh.shiftY =
   computescaleandshiftfromboundingbox(get(mesh.boundingBox))
end

function scalepoints(mesh::Mesh, points::Array{Float64, 2})
  local scaledX = 0.05*(points[:, 1] - get(mesh.shiftX))/get(mesh.scale) + 1.5
  local scaledY = 0.05*(points[:, 2] - get(mesh.shiftY))/get(mesh.scale) + 1.25
  local scaledPoints = [scaledX scaledY]

  return scaledPoints
end

function unscalepoints(mesh::Mesh, points::Array{Float64, 2})
  local unScaledX = (points[:, 1]  - 1.5) * get(mesh.scale)/0.05 + get(mesh.shiftX)
  local unScaledY = (points[:, 2] - 1.25) * get(mesh.scale)/0.05 + get(mesh.shiftY)
  local unScaledPoints = [unScaledX unScaledY]

  return unScaledPoints
end

"""
    push!(mesh::Mesh, points::Array{Float64, 2})

Inserts a number of vertices into the mesh.

# Remarks
* If the mesh was initialized without a bounding box, it is determined by
  the first set of vertices that are pushed into the triangualation. All
  consecutive push operations cannot insert vertices that exceed the bounding
  box limits
* New faces that are created by a push operation are marked as interior by
  default, i.e. it is forbidden to push vertices into an exterior region.
"""
function push!(mesh::Mesh, points::Array{Float64, 2})
  local minX = minimum(points[:,1])
  local maxX = maximum(points[:,1])
  local minY = minimum(points[:,2])
  local maxY = maximum(points[:,2])

  if isnull(mesh.boundingBox)
    setboundingbox(mesh, [minX, maxX, minY, maxY])
  else
    @assert (minX >= get(mesh.boundingBox)[1])
    @assert (maxX <= get(mesh.boundingBox)[2])
    @assert (minY >= get(mesh.boundingBox)[3])
    @assert (maxY <= get(mesh.boundingBox)[4])
  end

  local scaledPoints = scalepoints(mesh, points)
  push!(mesh.tesselation, scaledPoints)

  # mark faces as interior
  local noMissingFaceLocations =
   length(mesh.tesselation.faces) - length(mesh.faceLocation)
  local truth = fill(true, noMissingFaceLocations)
  mesh.faceLocation = [mesh.faceLocation; truth ]

  # invalidate voronoi vertex cache
  mesh.voronoiVertices = Nullable{Array{Float64, 2}}()
end

"""
    getfacesinsideregion(mesh::Mesh, ei::Int)

Finds all faces that are located inside the same region as the face identified
by `org(ei)`.

# Remarks
* The algorithm performs a breadth-first search over all adjacent triangles
  while omitting triangles that requires crossing a boundary.
"""
function getfacesinsideregion(mesh::Mesh, ei::Int)
  local foundFaces = Vector{Int}()
  local queuedEdges = [ei]
  while (length(queuedEdges) > 0)
    local currentEdge = shift!(queuedEdges)
    local currentFace = org(mesh.tesselation, currentEdge)

    # do nothing if it was already discovered
    if (currentFace ∉ foundFaces)
      push!(foundFaces, currentFace)

      # insert adjacent faces if they do not require boundary crossing
      local ec = currentEdge
      if (!isconstraint(mesh.tesselation, rot(ec)))
        push!(queuedEdges, sym(ec))
      end
      ec = onext(mesh.tesselation, ec)
      if (!isconstraint(mesh.tesselation, rot(ec)))
        push!(queuedEdges, sym(ec))
      end
      ec = onext(mesh.tesselation, ec)
      if (!isconstraint(mesh.tesselation, rot(ec)))
        push!(queuedEdges, sym(ec))
      end
    end
  end

  return foundFaces
end

"""
    addconstraint!(mesh::Mesh, vertices::Vector{Int})

Adds a constraint, given as a list of existing vertices, into the DelaunayMeshes.

# Remarks
* Constraints are closed polygons. A constraint edge from the last vertex to
  the first vertex is automatically inserted to ensure closure.
* All faces interior of the constrained region are marked as "exterior" to the
  DelaunayMeshes. A face is interior to the constrained region if it is located
  to the right of the given constraint polygon (i.e. the constraint vertices
  circle the interior region counter-clockwise).
* Constraint edges must not intersect. For performance reasons, this is not
  checked.
* Any vertex of the triangulation can intersect either zero or two constraint edges.
"""
function addconstraint!(mesh::Mesh, vertices::Vector{Int})
  # insert constraints between vertices and last one to close
  for i=1:length(vertices)-1
    insertconstraint!(mesh.tesselation, vertices[i], vertices[i+1])
  end
  local ci = insertconstraint!(mesh.tesselation, vertices[end], vertices[1])

  # find vertices to the right of the constraint
  local facesInside = getfacesinsideregion(mesh, rot(ci))
  mesh.faceLocation[facesInside] = false

  # invalidate voronoi vertex cache
  mesh.voronoiVertices = Nullable{Array{Float64, 2}}()
end

function istrianglevalid(mesh::Mesh, x::Tuple{Int, Int, Int})
  local rei = rot(x[1])
  local faceI = dest(mesh.tesselation, rei)
  return mesh.faceLocation[faceI]
end

function verticesfromtriangle(mesh::Mesh, x::Tuple{Int, Int, Int})
  local v1 = org(mesh.tesselation, x[1])
  local v2 = org(mesh.tesselation, x[2])
  local v3 = org(mesh.tesselation, x[3])

  return [v1 v2 v3]
end

"""
    convert(::Type{DiffEqPDEBase.SimpleFEMMesh}, mesh::Mesh)

Converts the given mesh to a format that can be used by the FEM-solvers.
"""
function convert(::Type{DiffEqPDEBase.SimpleFEMMesh}, mesh::Mesh)
  local allTriangles = DelaunayMeshes.getalltriangles(mesh.tesselation)
  local validTriangles = filter((x)->istrianglevalid(mesh,x), allTriangles)
  local triangles = map((x)->verticesfromtriangle(mesh,x), validTriangles)

  local redTriangles = reduce(vcat, triangles)
  local redVertices = reduce(hcat, mesh.tesselation.vertices)'

  return DiffEqPDEBase.SimpleFEMMesh(mesh.tesselation.vertices, redTriangles)
end
