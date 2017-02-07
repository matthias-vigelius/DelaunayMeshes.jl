typealias Point Vector{Float64}
typealias DelaunayTesselation SubDivision{Point}

minX = 1.0
maxX = 2.0 - eps(Float64)
deltaX = maxX - minX

"""
    DelaunayTesselation()

Creates an initial triangulation consisting of a maximal triangle and the
corresponding faces.

# Remarks
* The inner face coordinates are the centroid of the triangle.
* The outer face coordinates are located at ``(∞, ∞)``.
"""
function DelaunayTesselation()
  local v1 = convert(Point, [minX, minX])
  local v2 = convert(Point, [maxX, minX])
  local v3 = convert(Point, [minX + 0.5 * deltaX, maxX])
  local innerFace = convert(Point, (1./3.) * (v1 + v2 + v3))
  local outerFace = convert(Point, [Inf, Inf])
  local sd = create_subdivision(v1, v2, v3, innerFace, outerFace)

  return sd
end

"""
    getalltriangles(sd::SubDivision)

Returns all triangles as edge lists.

# Remarks
* The boundary edges are included.
"""
function getalltriangles(sd::DelaunayTesselation)
  local triangles = Vector{Tuple{Int, Int, Int}}()
  local nedges = length(sd.edges)
  for ei in 1:2:nedges
    local e2 = lnext(sd, ei)
    local e3 = lnext(sd, e2)
    if (e2 > ei && e3 > ei)
      push!(triangles, (ei, e2, e3))
    end
  end

  return triangles
end

function getdelaunaycoordinates(sd::DelaunayTesselation)
  local x = Vector{Float64}()
  local y = Vector{Float64}()
  for i in 1:4:length(sd.edges)
    local vo = org(sd, i)
    local vd = dest(sd, i)
    if (vo > 3 && vd > 3)
      push!(x, sd.vertices[vo][1])
      push!(x, sd.vertices[vd][1])
      push!(x, NaN)
      push!(y, sd.vertices[vo][2])
      push!(y, sd.vertices[vd][2])
      push!(y, NaN)
    end
  end

  return x, y
end

"""
    swap!(sd::DelaunayTesselation, ei::Int)

Swaps edge `ei` such that the new edge connects the apexes of the triangles
adjacent to the old edge.

# Remarks
* Specifially, if `a = oprev(ei)` and `b = oprev(sym(ei))` before the swap, then
  the edge `ei` will connect `dest(a)` with `dest(b)` after the swap.
"""
function swap!(sd::DelaunayTesselation, ei::Int)
  local es = sym(ei)
  local a = oprev(sd, ei)
  local b = oprev(sd, es)
  local leftFaceIndex = lface(sd, ei)
  local rightFaceIndex = rface(sd, ei)

  # update vertex cache
  sd.vertexCache[org(sd, ei)] = a
  sd.vertexCache[org(sd, es)] = b

  splice!(sd, ei, a)
  splice!(sd, es, b)
  splice!(sd, ei, (lnext(sd, a)))
  splice!(sd, es, (lnext(sd, b)))
  endpoints!(sd, ei, dest(sd, a), dest(sd, b), leftFaceIndex, rightFaceIndex)

  #update faces of adjacent triangles
  local ln = lnext(sd, ei)
  while (ln != ei)
    sd.edges[sym(rot(ln))].origin = leftFaceIndex
    ln = lnext(sd, ln)
  end
  local lnr = lnext(sd, es)
  while (lnr != es)
    sd.edges[sym(rot(lnr))].origin = rightFaceIndex
    lnr = lnext(sd, lnr)
  end
end

"""
    rightof(sd::DelaunayTesselation, vi::Int, ei::Int)
    rightof(sd::DelaunayTesselation, x::Point, ei::Int)

Computes the location of vertex `vi` relative to line defined by `ei`.

# Remarks
* Internally, the computation is forwarded to
  `GeometricalPredicates.orientation`.

# Returns
* True, if `vi` is located to the right of the edge
"""
function rightof(sd::DelaunayTesselation, x::Point, ei::Int)
  local op = sd.vertices[org(sd, ei)]
  local dp = sd.vertices[dest(sd, ei)]
  local o = GeometricalPredicates.Point2D(op[1], op[2])
  local d = GeometricalPredicates.Point2D(dp[1], dp[2])
  local edge = GeometricalPredicates.Line(o,d)
  local test = GeometricalPredicates.Point2D(x[1], x[2])

  return GeometricalPredicates.orientation(edge, test) == -1
end

function rightof(sd::DelaunayTesselation, vi::Int, ei::Int)
  local x = sd.vertices[vi]
  return rightof(sd, x, ei)
end

function closeto(sd::DelaunayTesselation, vi::Int, vj::Int)
  local vip = sd.vertices[vi]
  local vjp = sd.vertices[vj]
  return vip ≈ vjp
end

function edgevect(sd::DelaunayTesselation, ei::Int)
  local v1 = org(sd, ei)
  local v2 = dest(sd, ei)
  return (sd.vertices[v2])-(sd.vertices[v1])
end

function onlinesegment(
  xi::Point,
  x1::Point,
  x2::Point)
    if (xi ≈ x1) || (xi ≈ x2)
      return true
    end
    local d1i = normalize(xi - x1)
    local di2 = normalize(x2 - xi)
    return d1i ≈ di2
end

function onlinesegment(sd::DelaunayTesselation, vi::Int, v1::Int, v2::Int)
  local xi = sd.vertices[vi]
  local x1 = sd.vertices[v1]
  local x2 = sd.vertices[v2]
  if (vi == 284)
    local d1i = normalize(xi - x1)
    local di2 = normalize(x2 - xi)
  end
  return onlinesegment(xi, x1, x2)
end

function orientation(v1::Point, v2::Point, v3::Point)
  local xa = v1[1]
  local ya = v1[2]
  local xb = v2[1]
  local yb = v2[2]
  local xc = v3[1]
  local yc = v3[2]

  local det = (xa*yb - xa*yc - xb*ya + xb*yc + xc*ya - xc*yb)
  if abs(det) <= 1e-10
    return 0
  elseif det > 0.
    return 1
  else
    return 2
  end
end

function orientation(sd::DelaunayTesselation, v1::Int, v2::Int, v3::Int)
  return orientation(sd.vertices[v1], sd.vertices[v2], sd.vertices[v3])
end

"""
    pointsintersect(p1::Point, q1::Point, p2::Point, q2::Point)

Checks if line segments given by their end points intersect.

# Remarks
* Taken from geeksforgeeks.org/check-if-two-given-line-segments-intersect
"""
function pointsintersect(p1::Point, q1::Point, p2::Point, q2::Point)
  local o1 = orientation(p1,q1,p2)
  local o2 = orientation(p1,q1,q2)
  local o3 = orientation(p2,q2,p1)
  local o4 = orientation(p2,q2,q1)
  if (o1 != o2 && o3 != o4)
    return true
  elseif (o1 == 0 && onlinesegment(p2, p1, q1))
    return true
  elseif (o2 == 0 && onlinesegment(q2, p1, q1))
    return true
  elseif (o3 == 0 && onlinesegment(p1, p2, q2))
    return true
  elseif (o4 == 0 && onlinesegment(q1, p2, q2))
    return true
  else
    return false
  end
end

function pointsintersect(sd::DelaunayTesselation, p1::Int, q1::Int, p2::Int, q2::Int)
  return pointsintersect(sd.vertices[p1],sd.vertices[q1],sd.vertices[p2],sd.vertices[q2])
end

function edgeintersects(sd::DelaunayTesselation, e0::Int, v1::Int, v2::Int)
  return pointsintersect(sd, org(sd, e0), dest(sd, e0), v1, v2)
end

function edgeintersects(sd::DelaunayTesselation, e0::Int, e1::Int)
  local vo = org(sd, e1)
  local vd = dest(sd, e1)
  return edgeintersects( sd, e0, vo, vd)
end

"""
    locatevertex(sd::DelaunayTesselation, x::Point, ei::Int)

Locates point `x` in existing subdivision by traversing the triangulation
starting at edge `ei`.

# Remarks
* The point coordinates must be exactly inside the enclosing triangle.
* The algorithm used here [Brown]_ is more stable than the original version

.. [Brown] Brown, Faigle: A robust efficient algorithm for point location in triangulations
"""
function locatevertex(sd::DelaunayTesselation, x::Point, ei::Int)
  cei = ei
  if rightof(sd, x, ei)
    cei = sym(ei)
  end

  local count = 0

  while true
    count += 1
    if (count > length(sd.edges))
      throw(DomainError())
    end

    if x ≈ sd.vertices[org(sd, cei)] || x ≈ sd.vertices[dest(sd, cei)]
      return cei
    end

    local oni = onext(sd, cei)
    local odp = dprev(sd, cei)
    local isrightofOnext = rightof(sd, x, oni)
    local isrightofDprev = rightof(sd, x, odp)
    local op = convert(Int, !isrightofOnext) + 2 * convert(Int, !isrightofDprev)

    if op == 0
      return cei
    elseif op == 1
      cei = oni
    elseif op == 2
      cei = odp
    else
      local e0 = edgevect(sd, sym(oni))
      local xe = x - (sd.vertices[dest(sd, oni)])
      if (dot(e0,xe)) < 0.
        cei = odp
      else
        cei = oni
      end
    end
  end
end

function insertnewvertex!(sd::DelaunayTesselation, e0::Int, vi::Int)
  # add new edge
  local bi = makeedge!(sd)
  local se = bi
  local oldFace = lface(sd, e0)
  endpoints!(sd, bi, org(sd, e0), vi, oldFace, oldFace)

  # connect new edge to enclosing edges
  splice!(sd, bi, e0)

  local nei = e0
  local bni = bi

  while true
    local bsi = sym(bni)
    bni = connect!(sd, nei, bsi, [0.,0.])
    nei = oprev(sd, bni)
    if (lnext(sd, nei) == se)
      break
    end
  end

  return se, nei
end

"""
    insertpoint!(sd::DelaunayTesselation, x::Point, ei::Int)

Inserts a given point into an existing DelaunayMeshes.

# Remarks
* The algorithm identifies the triangle in the existing triangulation that contains
  the vertex to be inserted, taking <paramref name ="ei"/> as start point for the search.
  It then connects the newly-inserted point to all vertices of the enclosing triangle.
* If an existing vertex is already located at the point coordinates an error
  is thrown.

#Returns
* new edge, last edge of target triangle (lnext(lastEdge) = newEdge)
"""
function insertpoint!(sd::DelaunayTesselation, x::Point, ei::Int)
  # find triangle
  local e0 = locatevertex(sd, x, ei)

  # check if the vertex is already in the triangulation
  local eorg = org(sd, e0)
  local edest = dest(sd, e0)
  if (x ≈ sd.vertices[eorg] || x ≈ sd.vertices[edest])
    error("Vertex $x is already in DelaunayMeshes.")
  else
    # insert point as new vertex
    local vi = length(sd.vertices) + 1
    push!(sd.vertices, x)
    push!(sd.vertexCache, -1)

    # check if the point is on any edge of the triangle
    local onEdgeIndex = Nullable{EdgeIndex}()
    if (onlinesegment(sd, vi, eorg, edest))
      onEdgeIndex = Nullable{EdgeIndex}(e0)
    end
    if (isnull(onEdgeIndex))
      local e1 = lnext(sd, e0)
      if (onlinesegment(sd, vi, org(sd, e1), dest(sd, e1)))
        onEdgeIndex = Nullable{EdgeIndex}(e1)
      else
        local e2 = lnext(sd, e1)
        if (onlinesegment(sd, vi, org(sd, e2), dest(sd, e2)))
          onEdgeIndex = Nullable{EdgeIndex}(e2)
        end
      end
    end
    if (!isnull(onEdgeIndex))
      # if the edge is a segment edge we need to repair segment markers
      local ct = sd.edges[get(onEdgeIndex)].constraintType
      if (ct == NoConstraint)
        local ep = oprev(sd, get(onEdgeIndex))
        deleteedge!(sd, get(onEdgeIndex))
        newEdge = insertnewvertex!(sd, ep, vi)
        return newEdge
      else
        local lp = lprev(sd, get(onEdgeIndex))
        local ln = lnext(sd, get(onEdgeIndex))
        deleteedge!(sd, get(onEdgeIndex))
        newEdge = insertnewvertex!(sd, lp, vi)
        local aln = lnext(sd, lp)
        local apn = lprev(sd, ln)
        setconstraint!(sd, aln, ct)
        setconstraint!(sd, apn, ct)
        return newEdge
      end
    else
      newEdge = insertnewvertex!(sd, e0, vi)
      return newEdge
    end
  end
end

"""
    incircle(sd::DelaunayTesselation, v1::Int, v2::Int, v3::Int, v4::Int)

Tests if the vertex `v4` is strictly inside the circumcircle of the triangle
`v1`-`v2`-`v3`.
"""
function incircle(sd::DelaunayTesselation, v1::Int, v2::Int, v3::Int, v4::Int)
  local x1 = GeometricalPredicates.Point2D(sd.vertices[v1][1], sd.vertices[v1][2])
  local x2 = GeometricalPredicates.Point2D(sd.vertices[v2][1], sd.vertices[v2][2])
  local x3 = GeometricalPredicates.Point2D(sd.vertices[v3][1], sd.vertices[v3][2])
  local x4 = GeometricalPredicates.Point2D(sd.vertices[v4][1], sd.vertices[v4][2])
  local triangle = GeometricalPredicates.Primitive(x1,x2,x3)

  return GeometricalPredicates.incircle(triangle, x4) > 0
end

"""
    restoredelaunay!(sd::DelaunayTesselation, ei::Int, si::Int, xi::Int)

Restores Delaunay constraint for edge structure by swapping edges of adjacent triangles.

# Remarks
* It examines all triangles to check if the new vertex `xi` is inside its
  circumcircle.
* The algorithm cycles through all adjacent edges starting at `startEdge`
  and stopping at `stopEdge`.
* If an edge is swapped, its neighbours are examined as well.
"""
function restoredelaunay!(sd::DelaunayTesselation, startEdge::Int, stopEdge::Int, xi::Int)
  local curEdge = startEdge

  while true
    local t = oprev(sd, curEdge)
    local eo = org(sd, curEdge)
    local ed = dest(sd, curEdge)
    local td = dest(sd, t)

    if (rightof(sd, td, curEdge)
       && incircle(sd, eo, td, ed, xi)
       && !isconstraint(sd, curEdge))
      swap!(sd, curEdge)
      curEdge = oprev(sd, curEdge)
    else
      local on = onext(sd, curEdge)
      if (on == stopEdge)
        break
      end
      curEdge = lprev(sd, on)
    end
  end
end

function push!(sd::DelaunayTesselation, x::Array{Float64}, ei::Int)
    local newVertexIndex = length(sd.vertices) + 1
    stopEi, startEi = insertpoint!(sd, x, ei)
    restoredelaunay!(sd, startEi, stopEi, newVertexIndex)
    #=
    local faultyTriangles = TestHelpers.TestDelaunayness(sd)
    @assert(length(faultyTriangles) == 0)
    =#
    return startEi
end

function push!(sd::DelaunayTesselation, points::Array{Float64, 2})
  local startEdge = 1
  for i in 1:Int(length(points)/2)
    startEdge = push!(sd, [points[i, 1], points[i, 2]], startEdge)
  end
end

# for a given triangle, find edge that intersects with v1-v2
function findintersectionedgeoftriangle(
  sd::DelaunayTesselation, e0::Int, v1::Int, v2::Int)

  local ecur = e0
  while true
    if edgeintersects(sd, ecur, v1, v2)
      return ecur
    end

    ecur = lnext(sd, ecur)

    if (ecur == e0)
      throw(ErrorException("Cannot find intersecting edge for triangle $e0."))
    end
  end
end

"""
    findintersectingedges(sd::DelaunayTesselation, v1::Int, v2::Int, eg::Int)

Finds all edges of the given tesselation that intersect with a virtual edge
from `v1` to `v2`.
"""
function findintersectingedges(sd::DelaunayTesselation, v1::Int, v2::Int, eg::Int)
  # find starting edge
  local e0 = locatevertex(sd, sd.vertices[v1], eg)

  # org should be either e0 or lnext e0
  local es = org(sd, e0) == v1 ? e0 : lnext(sd, e0)

  # we need to find the existing edges that enclose the virtual edge
  local ef = es
  local rightOfTest = rightof(sd, v2, es)
  while true
    local nextEdge = onext(sd, ef)
    local rightOfNext = rightof(sd, v2, nextEdge)
    if rightOfNext && !rightOfTest
      break
    end
    ef = nextEdge
    rightOfTest = rightOfNext
  end

  # if enclosing edge has same dest, virtual edge is already in tesselation
  local intersecting = Vector{Int}()
  local alreadyIn = dest(sd, ef) == v2

  if (alreadyIn)
    push!(intersecting, ef)
  else
    while true
      local en = lnext(sd, ef)
      if dest(sd, en) == v2
        break
      end
      local eintersect = findintersectionedgeoftriangle(sd, en, v1, v2)
      ef = sym(eintersect)
      push!(intersecting, eintersect)
    end
  end

  return intersecting, alreadyIn
end

"""
    checkconvexityquadriliteral(sd::DelaunayTesselation, vertices::Vector{Int})

Checks if the given quadriliteral is strictly convex.

# Remarks
* The end points of the quadriliteral must be given in CCW order
"""
function checkconvexityquadriliteral(
  sd::DelaunayTesselation, vertices::Vector{Int})

  return pointsintersect(sd, vertices[1], vertices[3], vertices[2], vertices[4])
end

"""
    getadjacentquadriliteral(sd::DelaunayTesselation, ei::Int)

Gets vertices of quadriliteral formed by the two triangles adjacent to given edge in ccw order.
"""
function getadjacentquadriliteral(sd::DelaunayTesselation, ei::Int)
  local e1 = lnext(sd, ei)
  local e2 = lnext(sd, e1)
  local v2 = org(sd, ei)
  local v0 = org(sd, e1)
  local v1 = org(sd, e2)
  local es = sym(ei)
  local es1 = lnext(sd, es)
  local v3 = dest(sd, es1)
  [v0, v1, v2, v3]
end

"""
    removeintersectingedges!(
      sd::DelaunayTesselation, v1::Int, v2::Int, intersecting::Vector{Int})

Removes all intersecting edges by iteratively swapping diagonals of enclosing quadriliterals.

# Remarks
* For details see Sloan: A fast algorithm for generating constrained Delaunay triangulations.
* Returns the index of the constraint edge
"""
function removeintersectingedges!(
  sd::DelaunayTesselation,
  v1::Int, v2::Int,
  intersecting0::Vector{Int})

  local intersecting = copy(intersecting0)
  local lastSwapped = Nullable{Int}()

  while (length(intersecting) > 0)
    local i = shift!(intersecting)
    local adjacQuadriVert = getadjacentquadriliteral(sd, i)
    if !checkconvexityquadriliteral(sd, adjacQuadriVert)
      # put it at the back of list
      push!(intersecting, i)
    else
      # swap diagonals
      swap!(sd, i)

      # check if swapped edge intersects with virtual edge and put it on list
      # don't put it on the list if org/dest is either of v1/v2
      local io = org(sd, i)
      local id = dest(sd, i)

      if (io == v1) || (io == v2) || (id == v1) || (id == v2)
        lastSwapped = Nullable{Int}(i)
      else
        if (edgeintersects(sd, i, v1, v2))
          push!(intersecting, i)
        else
          lastSwapped = Nullable{Int}(i)
        end
      end
    end
  end

  return get(lastSwapped)
end

function restoredelaunayOfIntersecting!(
  sd::DelaunayTesselation, intersecting::Vector{Int}, v1::Int, v2::Int)

  while (length(intersecting) > 0)
    local i = shift!(intersecting)

    if ((org(sd, i) != v1) && (dest(sd, i) != v2)) && ((org(sd, i) != v2) && (dest(sd, i) != v1))

      local adj = getadjacentquadriliteral(sd, i)

      # check if one of each vertex is contained in other
      local t1in2 = incircle(sd, adj[1], adj[3], adj[4], adj[2])
      local t2in1 = incircle(sd, adj[1], adj[2], adj[3], adj[4])

      # if so swap and put new edge back on list of intersecting
      if t1in2 || t2in1
        swap!(sd, i)
        push!(intersecting, i)
      end
    end
  end
end

"""
    insertconstraint!(sd::DelaunayTesselation, v1::Int, v2::Int)

Inserts constraint going from `v1` to `v2`.

# Remarks
* The constrained region is assumed to be on the right of the directed
  edge.
"""
function insertconstraint!(sd::DelaunayTesselation, v1::Int, v2::Int)
    intersecting, alreadyIn = findintersectingedges(sd, v1, v2, 1)

    if alreadyIn
      local ci = intersecting[1]
      local co = org(sd, ci)
      local ct = v1 == co ? OnRightSide : OnLeftSide
      setconstraint!(sd, ci, ct)
      return ct == OnRightSide ? ci : sym(ci)
    else
      # edge is not already in - we need to swap
      ci = removeintersectingedges!(sd, v1, v2, intersecting)
      co = org(sd, ci)
      ct = v1 == co ? OnRightSide : OnLeftSide
      setconstraint!(sd, ci, ct)
      restoredelaunayOfIntersecting!(sd, intersecting, v1, v2)
      return ct == OnRightSide ? ci : sym(ci)
    end
end

"""
Represents a triangle. It can be constructed from three vertices or a starting edge.
It finds the edges enclosing the triangle and its inner face.
"""
type Triangle
  vertices::Tuple{VertexIndex, VertexIndex, VertexIndex}
  edges::Tuple{Int, Int, Int}
  face::VertexIndex

  function Triangle(sd::DelaunayTesselation, v1::VertexIndex, v2::VertexIndex, v3::VertexIndex)
    # get starting edge from vertex cache
    local ei = sd.vertexCache[v1]
    local e1 = onext(sd, ei)
    while (e1 != ei && dest(sd, e1) != v2)
      e1 = onext(sd, e1)
    end
    if (dest(sd, e1) != v2)
      error("Vertices $v1 - $v2 - $v3 do not belong to the same triangle")
    end
    local e2 = lnext(sd, e1)
    if (dest(sd, e2) != v3)
      error("Vertices $v1 - $v2 - $v3 do not belong to the same triangle")
    end
    local e3 = lnext(sd, e2)
    local face = dest(sd, rot(e1))
    new((v1, v2, v3), (e1, e2, e3), face)
  end

  function Triangle(sd::DelaunayTesselation, e1::EdgeIndex)
    local e2 = lnext(sd, e1)
    local e3 = lnext(sd, e2)
    local v1 = org(sd, e1)
    local v2 = org(sd, e2)
    local v3 = org(sd, e3)
    local face = dest(sd, rot(e1))
    new((v1, v2, v3), (e1, e2, e3), face)
  end
end

# Iterates all triangles
type TrianglesIteratorState
  curIndex::VertexIndex
  visited::Array{VertexIndex}
  curElement::Nullable{VertexIndex}
end

function findNextIteratorState!(sd::DelaunayTesselation, state::TrianglesIteratorState)
  # are there more triangles?
  state.curElement = Nullable{VertexIndex}()
  state.curIndex = state.curIndex + 2
  while (state.curIndex < length(sd.edges))
    cf = org(sd, state.curIndex)
    if (cf > 0 && cf ∉ state.visited)
      state.curElement = cf
      push!(state.visited, cf)
      break
    end
    state.curIndex = state.curIndex + 2
  end

  return state
end

function Base.start(sd::DelaunayTesselation)
  startState = TrianglesIteratorState(0, Array{VertexIndex, 1}(), Nullable{VertexIndex}())
  state = findNextIteratorState!(sd, startState)

  return state
end

function Base.done(sd::DelaunayTesselation, state::TrianglesIteratorState)
  return isnull(state.curElement)
end

function Base.next(sd::DelaunayTesselation, state::TrianglesIteratorState)
  ci = rot(state.curIndex)
  tri = Triangle(sd, ci)

  findNextIteratorState!(sd, state)

  return tri, state
end
