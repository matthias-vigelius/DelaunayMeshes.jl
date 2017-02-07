"""
Maximal ratio of circumcircle radius and length of shortest edge before  a triangle
is considered bad.
"""
beta = sqrt(2.)

minArea = 1e-5

"""
    compute_circumcenter_radius(mesh::Mesh, p::VertexIndex, r::VertexIndex, q::VertexIndex)

Computes the position and radius of the triangle circumcenter.

# Returns
* cx, cy, rs with ``c_x``, ``c_y`` the position and ``r_s`` the radius square of the circum center.
"""
function compute_circumcenter_radius(mesh::Mesh, p::VertexIndex, r::VertexIndex, q::VertexIndex)
  local pr = mesh.tesselation.vertices[p]
  local qr = mesh.tesselation.vertices[q]
  local rr = mesh.tesselation.vertices[r]
  local ax = pr[1]
  local ay = pr[2]
  local bx = qr[1]
  local by = qr[2]
  local cx = rr[1]
  local cy = rr[2]
  local sx = (ax^2.0*(by-cy) + ay^2.0*(by-cy) + bx^2.0*(cy-ay) + by^2.0*(cy-ay) + cx^2.0*(ay-by) + cy^2.0*(ay-by))/2.0
  local sy = (ax^2.0*(cx-bx) + ay^2.0*(cx-bx) + bx^2.0*(ax-cx) + by^2.0*(ax-cx) + cx^2.0*(bx-ax) + cy^2.0*(bx-ax))/2.0
  local a = ax*(by - cy) + bx * (cy-ay) + cx*(ay-by)
  local b = ax^2.0 * (bx*cy - by * cx) + ay^2.0 * (bx*cy - by*cx) + bx^2.0*(ay*cx-ax*cy) + by^2.0*(ay*cx-ax*cy) + cx^2.0*(ax*by-ay*bx) + cy^2.0*(ax*by-ay*bx)
  local rs = ((b/a + (sx*sx + sy*sy)/(a*a)))
  return sx/a, sy/a, rs
end

"""
    compute_offcenter(mesh:Mesh, p::VertexIndex, q::VertexIndex, r::VertexIndex)

Computes position of offcenter[^1] for given triangle ``p - q - r``.

# Remarks
* ``p - q`` must be the shortest edge of the triangle, i.e. the smallest angle is located at ``r``.

# References
[^1] Alper Üngör (2009). Off-centers: A new type of Steiner points for computing size-optimal quality-guaranteed Delaunay triangulations, Computational Geometry, *42* (2)

![SteinerPoint](SteinerPoint.svg)
"""
function compute_offcenter(mesh::Mesh, p::VertexIndex, q::VertexIndex, r::VertexIndex)
  local pr = mesh.tesselation.vertices[p]
  local qr = mesh.tesselation.vertices[q]
  local rr = mesh.tesselation.vertices[r]

  # position of bisection start of pq and circum center
  local bsr = 0.5 * (pr + qr)
  ccx, ccy, ccrs = compute_circumcenter_radius(mesh, p, q, r)

  # length of pq and radius of offcenter
  local lpq = norm(pr - qr)
  local ocr = beta * lpq * 0.999999

  # distance edge midpoint to center of offcircle ``c_2``
  local dobsoc = sqrt(ocr*ocr - (lpq/2.)*(lpq/2.) )

  # distance edge midpoint to center of outcircle
  local dbscc = norm(bsr - [ccx, ccy])
  if (dbscc <= dobsoc)
    # we take the out circle center
    return [ccx, ccy]
  else
    # we take off center

    # bisection: pq-midpoint
    local c2 = bsr + dobsoc * normalize([ccx, ccy] - bsr)
    local offcenter = c2 + ocr * normalize([ccx, ccy] - bsr)
    return offcenter
  end
end

"""
    vertex_encroaches_segment(mesh::Mesh, ei::Int, vertex::Point)

Checks if a potential vertex encroaches the segment denoted by `ei`.

# Remarks
* If edge `ei` is not a segment, i.e. not marked as a boundary, it always
  returns `false`
"""
function vertex_encroaches_segment(mesh::Mesh, ei::Int, vertex::Point)
    if (!isconstraint(mesh.tesselation, ei))
      return false
    end

    vio = org(mesh.tesselation, ei)
    vid = dest(mesh.tesselation, ei)
    v1 = mesh.tesselation.vertices[vio]
    v2 = mesh.tesselation.vertices[vid]
    mp = 0.5*(v1 + v2)
    dr = vertex - mp # from midpoint to new vertex
    dx = v1 - mp # from midpoint to end point

    return (dr[1]^2 + dr[2]^2) <= (dx[1]^2 + dx[2]^2)
  end


""""
    get_shortest_edge(mesh::Mesh, ai::VertexIndex, bi::VertexIndex, ci::VertexIndex)

Given a triangle with CCW-ordered vertices `ai` - `bi` - `ci`, it computes the
edge length and returns a sorted tuple `p` - `q` - `r` such that the triangle
``p - q - r`` has its shortest edge between `p` and `q`.
"""
function get_shortest_edge(mesh::Mesh, ai::VertexIndex, bi::VertexIndex, ci::VertexIndex)
  # nomenclature as in Bronstein et al..
  xa = mesh.tesselation.vertices[ai]
  xb = mesh.tesselation.vertices[bi]
  xc = mesh.tesselation.vertices[ci]
  vab = (xb - xa) # edge a
  vac = (xc - xa) # edge b
  vbc = (xb - xc) # edge a
  c = norm(vab)
  b = norm(vac)
  a = norm(vbc)
  #cosa = (b'*b + c'*c - a'*a)/(2. * b' * c)
  #cosb = (a'*a + c'*c - b'*b)/(2. * c' * a)
  #cosc = (a'*a + b'*b - c'*c)/(2. * a' * b)
  ordered = [(a, ai, bi, ci), (b, bi, ci, ai), (c, ci, ai, bi)]

  # order in place
  if (ordered[1][1] > ordered[2][1])
    ordered[1], ordered[2] = ordered[2], ordered[1]
  end
  if (ordered[2][1] > ordered[3][1])
    ordered[2], ordered[3] = ordered[3], ordered[2]
  end
  if (ordered[1][1] > ordered[2][1])
    ordered[1], ordered[2] = ordered[2], ordered[1]
  end

  return (ordered[1][2], ordered[1][3], ordered[1][4])
end

 """
    function refine_triangle(mesh::Mesh, tri::Triangle)

Refines the triangle located to the left of edge `ei` by adding a new vertex
at the offcenter or the midpoint of encroached edges.

# Remarks
* If a new vertex at the off-center of the triangle does not encroach any segments
  it is inserted at the off-center.
* Otherwise, a new vertex is inserted at the midpoint of each edge that the
  off-center would encroach on.
  try
"""
function refine_triangle(mesh::Mesh, tri::Triangle)

  # get ordered vertices
  p, q, r = get_shortest_edge(mesh, tri.vertices[1], tri.vertices[2], tri.vertices[3])

  # get offcenter
  offcenter = compute_offcenter(mesh, p, q, r)

  # get indices of all encroached segments
  indices = [i for i in 1:4:length(mesh.tesselation.edges)]
  encroachedIndices = filter(i->vertex_encroaches_segment(mesh, i, offcenter), indices)

  if isempty(encroachedIndices)
    # insert off-cente  r
    # we need to push the *unscaled point* so we push it straight to tesselation
    push_scaled!(mesh, offcenter')
  else
    # insert mid points
    for ei in encroachedIndices
      # get midpoint and insert it
      oi = org(mesh.tesselation, ei)
      di = dest(mesh.tesselation, ei)
      mp = 0.5 * (mesh.tesselation.vertices[oi] + mesh.tesselation.vertices[di])

      # we need to push the *unscaled point* so we push it straight to tesselation
      newVertexIndex = length(mesh.tesselation.vertices) + 1
      push_scaled!(mesh, mp')

      # we need to fix the face locations for the new faces.. we do that by
      # iterating over the star of the new vertex, starting at the segment
      startEdge = mesh.tesselation.vertexCache[newVertexIndex]
      while (!isconstraint(mesh.tesselation, startEdge))
        startEdge = onext(mesh.tesselation, startEdge)
      end
      faceIsInside = mesh.tesselation.edges[startEdge].constraintType == OnRightSide
      curEdge = startEdge
      while (true)
        # mark left face as inside/outside
        mesh.faceLocation[lface(mesh.tesselation, curEdge)] = faceIsInside
        curEdge = onext(mesh.tesselation, curEdge)
        if (curEdge == startEdge)
          break
        elseif (isconstraint(mesh.tesselation, curEdge))
          faceIsInside = !faceIsInside
        end
      end
    end
  end
end

"""
    check_triangle(mesh::Mesh, tri::Triangle)

Checks if the triangle satisfies the quality measure ``r/\bar{pq} \le \beta``,
with ``r`` the radius of the circumcircle and ``\bar{pq}`` the length of its
shortest edge.

# Returns
* `true`, if the triangle satisfies the quality measure.
"""
function check_triangle(mesh::Mesh, tri::Triangle)
  # TODO: we really need to get rid of the double computation of offcenter etc.
  #       maybe put a nullable in Triangle?

  # check if the area is below limit
  ra = mesh.tesselation.vertices[tri.vertices[1]]
  rb = mesh.tesselation.vertices[tri.vertices[2]]
  rc = mesh.tesselation.vertices[tri.vertices[3]]

  xa, ya = ra[1], ra[2]
  xb, yb = rb[1], rb[2]
  xc, yc = rc[1], rc[2]

  area = 0.5*abs((xa-xc)*(yb-ya)-(xa-xb)*(yc-ya))
  if area < minArea
    return true
  end

  # get ordered vertices and circumcircle
  p, q, r = get_shortest_edge(mesh, tri.vertices[1], tri.vertices[2], tri.vertices[3])
  sx, sy, rs = DelaunayMeshes.compute_circumcenter_radius(mesh, 1, 2, 3)

  lpq = norm(mesh.tesselation.vertices[p] - mesh.tesselation.vertices[q])
  return sqrt(rs)/lpq <= DelaunayMeshes.beta
end

"""
    refine_valid_triangles!(mesh::

Gets the first triangle that does not satisfy the quality measure or `null`.
"""
function refine_valid_triangles!(mesh::Mesh)
  foundBadTriangle = true
  counter = 0
  while foundBadTriangle && counter < 10000
    counter = counter + 1
    foundBadTriangle = false

    # get the first valid but bad triangle
    itState = start(mesh.tesselation)
    while (!done(mesh.tesselation, itState))
      tri, itState = next(mesh.tesselation, itState)
      if (mesh.faceLocation[tri.face] && !check_triangle(mesh, tri))

        # the triangle is bad and valid - we need to refine it
        refine_triangle(mesh, tri)
        foundBadTriangle = true

        # break out to outer loop
        break
      end
    end
  end
end
