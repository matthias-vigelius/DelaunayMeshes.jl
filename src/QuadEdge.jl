"
Implementation of the QuadEdge-Datastructure. [Guibas85]_

.. [Guibas85] Leonidas Guibas and Jorge Stolfi, Primitives for the manipulation of general subdivisions and the computation of Voronoi diagrams, *ACM Transactions on Graphics*, 4(2), 1985, 75-123.
"
@enum ConstraintType NoConstraint OnLeftSide OnRightSide

type Edge
  origin::Int
  next::Int
  constraintType::ConstraintType
end

type SubDivision{T}
  edges::Vector{Edge}
  vertexCache::Vector{Int}
  vertices::Vector{T}
  faces::Vector{T}
end

"""
    makeEdge!(sd::SubDivision)

Creates an empty edge and adds it to given SubDivision.

#returns
* Number of newly-added edge.
"""
function makeedge!(sd::SubDivision)
  local newEdgeNo = length(sd.edges) + 0x01
  push!(sd.edges, Edge(0, newEdgeNo, NoConstraint))
  push!(sd.edges, Edge(0, newEdgeNo+3, NoConstraint))
  push!(sd.edges, Edge(0, newEdgeNo+2, NoConstraint))
  push!(sd.edges, Edge(0, newEdgeNo+1, NoConstraint))
  return newEdgeNo
end

##
# basic edge index operations
"""
    base(ei::Int)

Computes base index of index `ei`.
"""
function base(ei::Int)
  return ei - ((ei-1) % 4)
end

"""
    sym(ei::Int)

Computes index of symmetric index `ei`.
"""
function sym(ei::Int)
  local qi = (ei-1) % 4
  if (qi < 2)
     return ei + 2
  else
     return ei - 2
   end
 end

"""
    rot(ei::Int)

Computes the index of the dual edge to index `ei`.

# Remarks
* The dual edge is obtained by rotating `ei` counter-clockwise by 90 degrees.
* The start(end) points of the dual edge can be interpreted as representing
  the right(left) face.
"""
function rot(ei::Int)
  local qi = (ei-1) % 4
  if (qi <3)
    return ei + 1
  else
    return ei - 3
  end
end

"""
    invrot(ei::Int)

Inverse operation of `rot(ei)`.
"""
function invrot(ei::Int)
  local qi = (ei-1) % 4
  if (qi == 0)
    return ei + 3
  else
    return ei - 1
  end
end

##
# Basic edge operations
"""
    org(sd::SubDivision, ei::Int)

Gets the index of the origin vertex of edge index `ei`.
"""
function org(sd::SubDivision, ei::Int)
  return sd.edges[ei].origin
end

"""
    dest(sd::SubDivision, ei::Int)

Gets the index of the destination vertex of edge index `ei`.
"""
function dest(sd::SubDivision, ei::Int)
  local di = sym(ei)
  return sd.edges[di].origin
end

"""
    onext(sd::SubDivision, ei::Int)

Gets the next edge starting at the origin of `ei` in counter-clockwise
direction.
"""
function onext(sd::SubDivision, ei::Int)
  return sd.edges[ei].next
end

"""
    oprev(sd::SubDivision, ei::Int)

Gets the next edge starting at the origin of `ei` in clockwise direction.
"""
function oprev(sd::SubDivision, ei::Int)
  return rot(onext(sd, rot(ei)))
end

"""
    dprev(sd::SubDivision, ei::Int)

Gets the next edge ending at the destination of `ei` in clockwise direction.
"""
function dprev(sd::SubDivision, ei::Int)
  return invrot(onext(sd, invrot(ei)))
end

"""
    lnext(sd::SubDivision, ei::Int)

Gets the next edge around the left face of `ei` and starting at destination
of `ei`.
"""
function lnext(sd::SubDivision, ei::Int)
  return rot(onext(sd, invrot(ei)))
end

"""
    lprev(sd::SubDivision, ei::Int)

Gets the previous edge around the left face of `ei` and ending at source
of `ei`.
"""
function lprev(sd::SubDivision, ei::Int)
  return sym(onext(sd, ei))
end

"""
    lface(sd::SubDivision, ei::Int)

Returns the index to face left of given edge `ei`.
"""
function lface(sd::SubDivision, ei::Int)
  return org(sd, sym(rot(ei)))
end


"""
    lface(sd::SubDivision, ei::Int)

Returns the index to face left of given edge `ei`.
"""
function rface(sd::SubDivision, ei::Int)
  return org(sd, rot(ei))
end

function lface_edges(sd::SubDivision, ei::Int)
  faceEdges = Vector{Int}()
  push!(faceEdges, ei)
  local currentEdge = lnext(sd, ei)
  while (currentEdge != ei)
    push!(faceEdges, currentEdge)
    currentEdge = lnext(sd, currentEdge)
  end

  return faceEdges
end

function lface_vertices(sd::SubDivision, ei::Int)
  local faceVertices = Vector{Int}()
  push!(faceVertices, org(sd, ei))

  local currentEdge = lnext(sd, ei)
  while (currentEdge != ei)
    push!(faceVertices, org(sd, currentEdge))
    currentEdge = lnext(sd, currentEdge)
  end

  return faceVertices
end

function ostar(sd::SubDivision, vi::Int)
  local ostarEdges = Vector{Int}()
  local startEdge = sd.vertexCache[vi]
  push!(ostarEdges, startEdge)
  local curEdge = onext(sd, startEdge)
  while (curEdge != startEdge)
    push!(ostarEdges, curEdge)
    curEdge = onext(sd, curEdge)
  end

  return ostarEdges
end

function isconstraint(sd::SubDivision, ei::Int)
  return !(sd.edges[ei].constraintType == NoConstraint)
end

function setconstraint!(sd::SubDivision, ei::Int, ct::ConstraintType)
  sd.edges[ei].constraintType = ct
  local sct = ct == OnRightSide ? OnLeftSide : OnRightSide
  sd.edges[sym(ei)].constraintType = sct
end

"""
    endpoints!(sd::SubDivision, ei::Int, vio::Int, vid::Int, lfi::Int, rfi::Int)

Assigns the start(end) point `vio`(`vid`) and the left(right) faces `lfi`(`rfi`)
to given edge `ei`.

# Remarks
* The start(end) points of the symmetric and symmetric dual edges are assigned
  correctly.
"""
function endpoints!(sd::SubDivision, ei::Int, vio::Int, vid::Int, lfi::Int, rfi::Int)
  local sei = sym(ei)

  # update end points
  sd.edges[ei].origin = vio
  sd.edges[sei].origin = vid

  # update end points of dual edges
  local re = rot(ei)
  local sre = sym(re)
  sd.edges[re].origin = rfi
  sd.edges[sre].origin = lfi

  # update vertex cache
  sd.vertexCache[vio] = ei
  sd.vertexCache[vid] = sei
end

"""
    splice!(sd::SubDivision, ai::Int, bi::Int)

Combines or splices two given edge rings `ai` and `bi`.

# Remarks
* If the origins of `ai` and `bi` are identical, it cuts the two edge rings
* If the origins differ, it joins the two edge rings
* After the operation, the new edge rings satisfy the identities
  `onext(ai) = onext(bi)` and `onext(bi) = onext(ai)` and the corresponding
  identities for the dual rings.
* The origins of the edges are neither evaluated nor are theyaffected by this
  operation and need to be updated manually.
"""
function splice!(sd::SubDivision, ai::Int, bi::Int)
  local ain = sd.edges[ai].next
  local bin = sd.edges[bi].next
  local alphai = rot(ain)
  local betai  = rot(bin)
  local alphain = sd.edges[alphai].next
  local betain  = sd.edges[betai].next

  sd.edges[ai].next = bin
  sd.edges[bi].next = ain
  sd.edges[alphai].next = betain
  sd.edges[betai].next = alphain
end

function create_subdivision{T}(v1::T, v2::T, v3::T, innerFace::T, outerFace::T)
  local sd = SubDivision(Array{Edge}(0), [-1,-1,-1], [v1, v2, v3],[innerFace, outerFace])

  # assign edges
  local ea = makeedge!(sd)
  endpoints!(sd, ea, 1, 2, 1, 2)

  local eb = makeedge!(sd)
  endpoints!(sd, eb, 2, 3, 1, 2)
  splice!(sd, sym(ea), eb)

  local ec = makeedge!(sd)
  endpoints!(sd, ec, 3, 1, 1, 2)
  splice!(sd, sym(eb), ec)
  splice!(sd, sym(ec), ea)

  return sd
end

###
# High level functions
"""
    connect!{T}(sd::SubDivision, ai::Int, bi::Int, newFace::T)

Connects the end vertex of edge `ai` to the start vertex `bi`.

# Remarks
* This operation adds another face to the subdivision. The new face will be
  to the right of the new edge and will be correctly assigned to all adjacent
  edges.
"""
function connect!{T}(sd::SubDivision{T}, ai::Int, bi::Int, newFace::T)
  local newEdge = makeedge!(sd)
  local newOrigin = dest(sd, ai)
  local newdestination = org(sd, bi)
  local oldFaceIndex = dest(sd, rot(ai))
  local newFaceIndex = length(sd.faces) + 1
  endpoints!(sd, newEdge, newOrigin, newdestination, oldFaceIndex, newFaceIndex)

  splice!(sd, newEdge, lnext(sd, ai))
  splice!(sd, sym(newEdge), bi)

  # update vertex cache
  sd.vertexCache[newOrigin] = newEdge
  sd.vertexCache[newdestination] = sym(newEdge)

  # add new face
  push!(sd.faces, newFace)

  local symNewEdge = sym(newEdge)
  local newFaceEdgeIt = lnext(sd, symNewEdge)
  while (newFaceEdgeIt != symNewEdge)
    sd.edges[sym(rot(newFaceEdgeIt))].origin = newFaceIndex
    newFaceEdgeIt = lnext(sd, newFaceEdgeIt)
  end

  return newEdge
end

function deleteedge!(sd::SubDivision, ei::Int)
  local op = oprev(sd, ei)
  local es = sym(ei)
  local esp = oprev(sd, es)

  local ov = org(sd, ei)
  local dv = org(sd, es)

  # remember faces
  local keepFaceIndex = lface(sd, ei)
  local deleteFaceIndex = rface(sd, ei)

  splice!(sd, ei, op)
  splice!(sd, es, esp)

  # mark edge as invalid
  local baseIndex = base(ei)
  sd.edges[baseIndex] = Edge(-1, -1, NoConstraint)
  sd.edges[baseIndex + 1] = Edge(-1, -1, NoConstraint)
  sd.edges[baseIndex + 2] = Edge(-1, -1, NoConstraint)
  sd.edges[baseIndex + 3] = Edge(-1, -1, NoConstraint)

  # update vertex cache
  sd.vertexCache[ov] = op
  sd.vertexCache[dv] = esp

  # update faces
  local keepFaceEdges = lface_edges(sd, op)
  for curEdge in keepFaceEdges
    sd.edges[sym(rot(curEdge))].origin = keepFaceIndex
  end
end

"""
    removedeletededges!(sd::SubDivision)

Removes all edges that were marked as deleted in DelaunayMeshes.
"""
function removedeletededges!(sd::SubDivision)

  # compute partial sum of offsets
  local delEdgeIndex = find((x::Edge) -> x.origin == -1, sd.edges[1:4:end])
  local offsets = zeros(Int, length(sd.edges[1:4:end]))
  offsets[delEdgeIndex] = 4
  cumsum!(offsets, offsets)

  # substract offset of edges
  for (i, e) in enumerate(sd.edges)
    local baseOffset = Int((base(e.next) - 1)/4 + 1)
    e.next -= offsets[baseOffset]
  end

  # remove all offending edges
  filter!((x::Edge) -> x.origin != -1, sd.edges)

  # rebuild vertex cache
  fill!(sd.vertexCache, -1)
  for i in 1:2:length(sd.edges)
    local org = sd.edges[i].origin
    if (sd.vertexCache[org] == -1)
      sd.vertexCache[org] = i
    end
  end
end
