function testvertexcache(sd::DelaunayMeshes.SubDivision)
  for (vi, v) in enumerate(sd.vertexCache)
    @assert v >= 1
    @assert sd.edges[v].origin == vi
  end
end

function setuptwopointtesselation(v1::DelaunayMeshes.Point, v2::DelaunayMeshes.Point)
    local tess = DelaunayMeshes.DelaunayTesselation()

    # push two more points and connect them
    push!(tess.vertices, v1)
    push!(tess.vertices, v2)
    push!(tess.vertexCache, -1)
    push!(tess.vertexCache, -1)

    local ea = 1
    local eb = 5
    local ec = 9
    local ed = DelaunayMeshes.makeedge!(tess)
    DelaunayMeshes.endpoints!(tess, ed, 1, 4, 1, 1)
    DelaunayMeshes.splice!(tess, ed, ea)
    local ee = DelaunayMeshes.connect!(tess, ea, DelaunayMeshes.Sym(ed), [1.0,1.0])
    local ef = DelaunayMeshes.connect!(tess, eb, DelaunayMeshes.Sym(ee), [2.0,2.0])

    local eg = DelaunayMeshes.makeedge!(tess)
    DelaunayMeshes.endpoints!(tess, eg, 4, 5, 3, 3)
    DelaunayMeshes.splice!(tess, eg, DelaunayMeshes.Sym(ee))
    local eh = DelaunayMeshes.connect!(tess, DelaunayMeshes.sym(ee), DelaunayMeshes.sym(eg), [3.0,3.0])
    local ei = DelaunayMeshes.connect!(tess, eb, DelaunayMeshes.sym(eh), [4.0,4.0])

    #=
    for (i,e) in enumerate(tess.edges)
      println("$i, $e")
    end
    =#
    return tess
end

function locateandassertintriangle(
  sd::DelaunayMeshes.DelaunayTesselation, x::DelaunayMeshes.Point, ei::Int)
  local e1 = DelaunayMeshes.locatevertex(sd, x, ei)
  local e2 = DelaunayMeshes.lnext(sd, e1)
  local e3 = DelaunayMeshes.lnext(sd, e2)
  #println("$e1")

  return  !DelaunayMeshes.rightof(sd, x, e1) && !DelaunayMeshes.RightOf(sd, x, e2) && !DelaunayMeshes.RightOf(sd, x, e3)
end

function testdelaunayness(sd::DelaunayMeshes.DelaunayTesselation)
  local triangles = DelaunayMeshes.getalltriangles(sd)
  local faultyTriangles = Vector{Tuple{Int, Int, Int, Int}}()

  for curTriangle in triangles
    local v1 = DelaunayMeshes.org(sd, curTriangle[1])
    local v2 = DelaunayMeshes.org(sd, curTriangle[2])
    local v3 = DelaunayMeshes.org(sd, curTriangle[3])

    if (v1 > 3 && v2 > 3 && v3 > 3)
      for (v, vv) in enumerate(sd.vertices)
        if (v != v1 && v != v2 && v != v3)
          if (DelaunayMeshes.incircle(sd, v1, v2, v3, v))
            push!(faultyTriangles, (v1, v2, v3, v))
          end
        end
      end
    end
  end

  if (length(faultyTriangles) > 0)
    for f in faultyTriangles
      println("$f: $(sd.vertices[f[1]]), $(sd.vertices[f[2]]), $(sd.vertices[f[3]]), $(sd.vertices[f[4]])")
    end
  end

  return faultyTriangles
end

function insertsingleedge!(sd::DelaunayMeshes.SubDivision, v0::Int, v1::Int, lfi::Int, rfi::Int)
  local newEdge = DelaunayMeshes.makeedge!(sd)
  DelaunayMeshes.endpoints!(sd, newEdge, v0, v1, lfi, rfi)
  return newEdge
end

function edgeconnectsvertices(sd::DelaunayMeshes.SubDivision, e0::Int, v0::Int, v1::Int)
  local vo = DelaunayMeshes.org(sd, e0)
  local vd = DelaunayMeshes.dest(sd, e0)
  return (vo == v0 && vd == v1) || (vo == v1 && vd == v0)
end

function findedgeconnectingvertices(sd::DelaunayMeshes.SubDivision, v0::Int, v1::Int)
  local estart = sd.vertexCache[v0]
  local ecur = DelaunayMeshes.onext(sd, estart)

  while (DelaunayMeshes.dest(sd, ecur) != v1 && ecur != estart)
    ecur = DelaunayMeshes.onext(sd, ecur)
  end

  @assert DelaunayMeshes.dest(sd, ecur) == v1
  return ecur
end

function validateconstraints(sd::DelaunayMeshes.SubDivision, constraints::Array{Int})

  # collect all ordered constraints from tesselation
  local foundConstraints = Array{Int, 2}()
  local baseEdges = Array{Int}()
  for i=1:4:length(sd.edges)
    if (DelaunayMeshes.isconstraint(sd, i))
      local baseInd = sd.edges[i].constraintType == DelaunayMeshes.OnRightSide ? i : i + 2
      local ov = DelaunayMeshes.org(sd, baseInd)
      local dv = DelaunayMeshes.dest(sd, baseInd)
      if (size(foundConstraints, 1) == 0)
        foundConstraints = [ov dv]
        baseEdges = [baseInd]
      else
        foundConstraints = [foundConstraints; ov dv]
        baseEdges = [baseEdges; baseInd]
      end
    end
  end

  local sortedFoundConstraints = sortrows(foundConstraints)
  @test sortedFoundConstraints == sortrows(constraints)

  # check that constraints are not intersected by any other edges
  for bi in baseEdges
    local bov = DelaunayMeshes.org(sd, bi)
    local bdv = DelaunayMeshes.dest(sd, bi)
    for i=1:4:length(sd.edges)
      local tov = DelaunayMeshes.org(sd, i)
      local tdv = DelaunayMeshes.dest(sd, i)
      # check for intersection if they do not have common vertices
      if (bov != tov && bov != tdv && bdv != tov && bdv != tdv)
         @test !(DelaunayMeshes.edgeintersects(sd, bi, i))
      end
    end
  end
end

function getrandomnumbersinsideboundingbox(bb::Vector{Float64}, n::Int)
  local sx = (bb[2] - bb[1])
  local sy = (bb[4] - bb[3])
  local pointsX = rand(Float64, n, 1) * sx + bb[1]
  local pointsY = rand(Float64, n, 1) * sy + bb[3]
  return [pointsX pointsY]
end

function pointsinconstrainedregion(
  mesh::DelaunayMeshes.Mesh, points::Array{Float64, 2})

  local sb = DelaunayMesh.scalepoints(mesh,
    [
      -1.0 -1.0;
      -0.5 -1.0;
      -0.5 0.0;
      -1.0 0.0;
      -1.0 0.0;
       1.0 0.0;
       1.0 1.0;
      -1.0 1.0;
       0.5 -1.0;
       1.0 -1.0;
       1.0 1.0;
       0.5 1.0;
    ]
  )

  local npoints = size(points, 1)
  pointsInRegion = fill(false, npoints)

  for i in 1:npoints
     local cp = points[i,:]
     local inSquare1 =
     [
       cp[1] >= sb[1,1] cp[2] >= sb[1,2];
       cp[1] <= sb[2,1] cp[2] >= sb[2,2];
       cp[1] <= sb[3,1] cp[2] <= sb[3,2];
       cp[1] >= sb[4,1] cp[2] <= sb[4,2];
     ]
    if (inSquare1 == trues(4,2))
      pointsInRegion[i] = true
    end
    local inSquare2 =
     [
       cp[1] >= sb[5,1] cp[2] >= sb[5,2];
       cp[1] <= sb[6,1] cp[2] >= sb[6,2];
       cp[1] <= sb[7,1] cp[2] <= sb[7,2];
       cp[1] >= sb[8,1] cp[2] <= sb[8,2];
     ]
    if (inSquare2 == trues(4,2))
      pointsInRegion[i] = true
    end
    local inSquare3 =
     [
       cp[1] >= sb[9,1] cp[2] >= sb[9,2];
       cp[1] <= sb[10,1] cp[2] >= sb[10,2];
       cp[1] <= sb[11,1] cp[2] <= sb[11,2];
       cp[1] >= sb[12,1] cp[2] <= sb[12,2];
     ]
    if (inSquare3 == trues(4,2))
      pointsInRegion[i] = true
    end
  end

  return pointsInRegion
end

function getplotedges(sd::DelaunayMeshes.DelaunayTesselation)
  local x = Vector{Float64}()
  local y = Vector{Float64}()
  for i in 1:4:length(sd.edges)
    local vo = DelaunayMeshes.org(sd, i)
    local vd = DelaunayMeshes.dest(sd, i)
    if (vo > -1)
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
