  """
Maximal ratio of circumcircle radius and length of shortest edge before  a triangle
is considered bad.
"""
beta = sqrt(2.)

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
  local c2 = bsr + dobsoc * normalize([ccx, ccy] - bsr)

  # bisection: pq-midpoint
  local offcenter = c2 + ocr * normalize([ccx, ccy] - bsr)
  return offcenter
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
